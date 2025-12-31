#!/usr/bin/env python3
"""
Batch MD simulations for all valid receptor-ligand pairs.
Runs on 8 GPUs in parallel using OpenCL.
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

def run_md_simulation(pair_dir, gpu_id, short_run=True):
    """Run MD simulation for a single pair on specified GPU."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer
    
    pair_name = pair_dir.name
    md_dir = pair_dir / "md"
    md_dir.mkdir(exist_ok=True)
    
    # Find best Chai-1 pose
    chai_dir = pair_dir / "chai1"
    best_cif = chai_dir / "pred.model_idx_0.cif"
    
    if not best_cif.exists():
        return False, "No CIF file"
    
    # Check if already done
    if (md_dir / "trajectory.dcd").exists():
        return True, "Already done"
    
    try:
        # Load structure with PDBFixer
        fixer = PDBFixer(str(best_cif))
        
        # Fix structure
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        
        # Set up force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller
        modeller = app.Modeller(fixer.topology, fixer.positions)
        
        # Add solvent (smaller box for speed)
        modeller.addSolvent(
            forcefield,
            model='tip3p',
            padding=0.6 * unit.nanometers,
            ionicStrength=0.15 * unit.molar
        )
        
        # Create system
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
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
            modeller.topology,
            system,
            integrator,
            platform,
            properties
        )
        simulation.context.setPositions(modeller.positions)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Save minimized structure
        state = simulation.context.getState(getPositions=True)
        with open(md_dir / "minimized.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Short equilibration run (10ps for testing, increase for production)
        if short_run:
            steps = 5000  # 10ps
        else:
            steps = 500000  # 1ns
        
        # Add reporters
        simulation.reporters.append(
            app.DCDReporter(str(md_dir / "trajectory.dcd"), 1000)
        )
        simulation.reporters.append(
            app.StateDataReporter(
                str(md_dir / "energy.csv"),
                1000,
                step=True,
                potentialEnergy=True,
                temperature=True
            )
        )
        
        # Run dynamics
        simulation.step(steps)
        
        # Save final state
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        with open(md_dir / "final.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        return True, f"{steps} steps"
        
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
        success, msg = run_md_simulation(pair_dir, gpu_id)
        elapsed = time.time() - start
        
        result_queue.put((success, pair_name, msg, elapsed, gpu_id))

def main():
    mp.set_start_method('spawn', force=True)
    
    output_base = Path('all_pairs_results')
    
    # Load valid pairs
    with open(output_base / 'valid_for_md.json') as f:
        valid_pairs = json.load(f)
    
    log(f"Valid pairs for MD: {len(valid_pairs)}")
    
    # Get pair directories
    tasks = []
    for pair_name in valid_pairs:
        pair_dir = output_base / pair_name
        if pair_dir.exists() and (pair_dir / "chai1").exists():
            # Skip if already done
            if not (pair_dir / "md" / "trajectory.dcd").exists():
                tasks.append(pair_dir)
    
    log(f"Tasks to run: {len(tasks)}")
    
    if not tasks:
        log("All MD simulations already complete!")
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
    
    # Progress file
    progress_file = output_base / 'md_progress.json'
    
    for i in range(len(tasks)):
        success, pair_name, msg, elapsed, gpu_id = result_queue.get()
        
        if success:
            success_count += 1
            log(f"[GPU{gpu_id}] ✓ {pair_name}: {msg} ({elapsed:.0f}s)")
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
        
        if completed % 100 == 0:
            log(f"Progress: {completed}/{len(tasks)} ({success_count} ok, {error_count} err) | Rate: {rate:.0f}/hr | ETA: {eta:.1f}h")
    
    # Wait for workers
    for p in workers:
        p.join()
    
    log(f"\nDone! Success: {success_count}, Errors: {error_count}")

if __name__ == '__main__':
    main()
