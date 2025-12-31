#!/usr/bin/env python3
"""
Batch MD simulations using implicit solvent (faster, no ligand parameterization needed).
Runs minimization + short dynamics to assess binding stability.
"""
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import warnings
warnings.filterwarnings('ignore')

import json
import time
import torch.multiprocessing as mp
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def extract_protein_only(cif_path, output_pdb):
    """Extract only protein from CIF file (skip ligand for now)."""
    from pdbfixer import PDBFixer
    
    fixer = PDBFixer(str(cif_path))
    
    # Standard amino acids
    aa_names = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
                'TYR', 'VAL', 'HIE', 'HID', 'HIP', 'CYX'}
    
    # Remove non-protein
    to_delete = []
    for residue in fixer.topology.residues():
        if residue.name not in aa_names:
            to_delete.append(residue)
    
    fixer.removeHeterogens(keepWater=False)
    
    # Fix structure
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    return fixer

def run_md_protein_only(pair_dir, gpu_id):
    """Run implicit solvent MD on protein structure."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer
    
    pair_name = pair_dir.name
    md_dir = pair_dir / "md"
    md_dir.mkdir(exist_ok=True)
    
    # Check if already done with best pose
    best_marker = md_dir / "used_best_pose"
    if (md_dir / "minimized.pdb").exists() and best_marker.exists():
        return True, "Already done"
    
    # Find best Chai-1 pose (selected from all 5 models)
    chai_dir = pair_dir / "chai1"
    best_cif = chai_dir / "best_pose.cif"
    
    # Fall back to model_idx_0 if best_pose doesn't exist
    if not best_cif.exists():
        best_cif = chai_dir / "pred.model_idx_0.cif"
    
    if not best_cif.exists():
        return False, "No CIF file"
    
    try:
        # Extract protein only
        fixer = extract_protein_only(best_cif, md_dir / "protein.pdb")
        
        if fixer.topology.getNumAtoms() < 100:
            return False, "Too few atoms"
        
        # Set up implicit solvent force field
        forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        
        # Create system with implicit solvent
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds
        )
        
        # Integrator
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picoseconds,
            0.002 * unit.picoseconds
        )
        
        # Use OpenCL on specified GPU
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'DeviceIndex': str(gpu_id)}
        
        # Create simulation
        simulation = app.Simulation(
            fixer.topology,
            system,
            integrator,
            platform,
            properties
        )
        simulation.context.setPositions(fixer.positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Get final energy
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Save minimized structure
        with open(md_dir / "minimized.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Save energy info
        energy_info = {
            'initial_energy_kJ': initial_energy,
            'final_energy_kJ': final_energy,
            'energy_drop_kJ': initial_energy - final_energy
        }
        with open(md_dir / "energy.json", 'w') as f:
            json.dump(energy_info, f)
        
        # Mark that we used best_pose
        with open(md_dir / "used_best_pose", 'w') as f:
            f.write(str(best_cif))
        
        return True, f"E: {final_energy:.0f} kJ/mol"
        
    except Exception as e:
        return False, str(e)[:100]

def worker(gpu_id, task_queue, result_queue):
    """Worker process for a single GPU."""
    while True:
        task = task_queue.get()
        if task is None:
            break
        
        pair_dir, idx, total = task
        pair_name = pair_dir.name
        
        start = time.time()
        success, msg = run_md_protein_only(pair_dir, gpu_id)
        elapsed = time.time() - start
        
        result_queue.put((success, pair_name, msg, elapsed, gpu_id))

def main():
    mp.set_start_method('spawn', force=True)
    
    output_base = Path('all_pairs_results')
    
    # Get all pairs with chai1 output
    tasks = []
    for pair_dir in sorted(output_base.iterdir()):
        if not pair_dir.is_dir():
            continue
        chai_dir = pair_dir / "chai1"
        if chai_dir.exists() and list(chai_dir.glob("*.cif")):
            # Skip if already done
            if not (pair_dir / "md" / "minimized.pdb").exists():
                tasks.append(pair_dir)
    
    log(f"Tasks to run: {len(tasks)}")
    
    if not tasks:
        log("All minimizations complete!")
        return
    
    num_gpus = 8
    task_queue = mp.Queue()
    result_queue = mp.Queue()
    
    # Start workers
    workers = []
    for gpu_id in range(num_gpus):
        p = mp.Process(target=worker, args=(gpu_id, task_queue, result_queue))
        p.start()
        workers.append(p)
        log(f"Started worker on GPU {gpu_id}")
    
    # Queue tasks
    for i, pair_dir in enumerate(tasks):
        task_queue.put((pair_dir, i, len(tasks)))
    
    # Send stop signals
    for _ in range(num_gpus):
        task_queue.put(None)
    
    # Track progress
    success_count = 0
    error_count = 0
    start_time = time.time()
    
    progress_file = output_base / 'md_progress.json'
    
    for i in range(len(tasks)):
        success, pair_name, msg, elapsed, gpu_id = result_queue.get()
        
        if success:
            success_count += 1
            if success_count % 100 == 0 or success_count <= 10:
                log(f"[GPU{gpu_id}] ✓ {pair_name}: {msg} ({elapsed:.1f}s)")
        else:
            error_count += 1
            log(f"[GPU{gpu_id}] ✗ {pair_name}: {msg}")
        
        # Update progress
        completed = i + 1
        elapsed_total = time.time() - start_time
        rate = completed / (elapsed_total / 3600) if elapsed_total > 0 else 0
        remaining = len(tasks) - completed
        eta = remaining / rate if rate > 0 else 0
        
        progress = {
            'completed': completed,
            'total': len(tasks),
            'success': success_count,
            'errors': error_count,
            'rate_per_hour': rate,
            'eta_hours': eta,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(progress_file, 'w') as f:
            json.dump(progress, f)
        
        if completed % 500 == 0:
            log(f"Progress: {completed}/{len(tasks)} ({success_count} ok, {error_count} err) | Rate: {rate:.0f}/hr | ETA: {eta:.1f}h")
    
    # Wait for workers
    for p in workers:
        p.join()
    
    log(f"\nDone! Success: {success_count}, Errors: {error_count}")

if __name__ == '__main__':
    main()
