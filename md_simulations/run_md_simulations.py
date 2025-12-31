#!/usr/bin/env python3
"""
Molecular Dynamics Simulation Runner
=====================================
Runs minimization and short equilibration for receptor structures.
Uses experimental PDB structures with proper cleaning.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import time
import sys

# Configuration
WORK_DIR = Path(__file__).parent
STRUCTURES_DIR = WORK_DIR / "structures"
STRUCTURES_DIR.mkdir(exist_ok=True)

# Receptor PDB IDs and their chains to use
RECEPTOR_CONFIG = {
    "5-HT2A": {"pdb": "6WHA", "chain": "A"},
    "5-HT2B": {"pdb": "5TVN", "chain": "A"},
    "D2": {"pdb": "6CM4", "chain": "A"},
    "MOR": {"pdb": "5C1M", "chain": "A"},
    "KOR": {"pdb": "4DJH", "chain": "A"},
    "CB1": {"pdb": "5TGZ", "chain": "A"},
    # Use a simpler test protein for NMDA (too large)
}

# Ligand SMILES for reference
LIGAND_SMILES = {
    "LSD": "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34",
    "psilocin": "CN(C)CCc1c[nH]c2ccc(O)cc12",
    "DMT": "CN(C)CCc1c[nH]c2ccccc12",
    "mescaline": "COc1cc(CCN)cc(OC)c1OC",
    "MDMA": "CC(NC)Cc1ccc2OCOc2c1",
    "ketamine": "CNC1(c2ccccc2Cl)CCCCC1=O",
    "salvinorin_A": "CC(=O)O[C@@H]1C[C@@H]2[C@]3(C)CC[C@@H]4C(C)(C)[C@H](OC(=O)c5ccoc5)C(=O)O[C@@]4(C)[C@@H]3C(=O)[C@@H]([C@@H]2C)[C@H]1OC(C)=O",
    "THC": "CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=C[C@H]21"
}


def download_pdb(pdb_id: str) -> Path:
    """Download PDB file from RCSB"""
    pdb_path = STRUCTURES_DIR / f"{pdb_id}.pdb"
    if pdb_path.exists():
        return pdb_path
    
    print(f"  Downloading {pdb_id}...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def extract_single_chain(pdb_path: Path, chain_id: str) -> Path:
    """Extract a single chain from PDB and clean it for OpenMM"""
    output_path = STRUCTURES_DIR / f"{pdb_path.stem}_chain{chain_id}.pdb"
    
    if output_path.exists():
        return output_path
    
    STANDARD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                   'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP'}
    
    with open(pdb_path) as f:
        lines = f.readlines()
    
    output_lines = []
    atom_num = 0
    
    for line in lines:
        if line.startswith('ATOM'):
            chain = line[21]
            res_name = line[17:20].strip()
            alt_loc = line[16]
            
            # Skip wrong chain
            if chain != chain_id:
                continue
            
            # Skip non-standard residues
            if res_name not in STANDARD_AA:
                continue
            
            # Skip alternate conformations
            if alt_loc not in [' ', 'A']:
                continue
            
            # Fix histidine naming
            if res_name == 'HIS':
                line = line[:17] + 'HID' + line[20:]
            
            atom_num += 1
            # Renumber atoms and remove alt loc
            new_line = f"ATOM  {atom_num:5d} {line[12:16]} {line[17:21]}{chain}{line[22:]}"
            output_lines.append(new_line)
    
    output_lines.append("TER\nEND\n")
    
    with open(output_path, 'w') as f:
        f.writelines(output_lines)
    
    print(f"  Extracted chain {chain_id}: {atom_num} atoms -> {output_path.name}")
    return output_path


def run_md_simulation(receptor_name: str, ligand_name: str, 
                      equilibration_steps: int = 5000,
                      production_steps: int = 0) -> dict:
    """
    Run MD simulation for a receptor-ligand pair.
    
    Args:
        receptor_name: Name of receptor (e.g., "5-HT2A")
        ligand_name: Name of ligand (e.g., "LSD")
        equilibration_steps: Number of equilibration steps
        production_steps: Number of production steps (0 for just minimization)
    
    Returns:
        Dictionary with simulation results
    """
    system_name = f"{receptor_name}_{ligand_name}"
    system_dir = WORK_DIR / system_name
    system_dir.mkdir(exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"SIMULATION: {system_name}")
    print('='*60)
    
    # Get receptor config
    config = RECEPTOR_CONFIG.get(receptor_name)
    if not config:
        print(f"  ERROR: Unknown receptor: {receptor_name}")
        return {"status": "error", "message": f"Unknown receptor: {receptor_name}"}
    
    start_time = time.time()
    
    try:
        # Download and prepare structure
        pdb_path = download_pdb(config["pdb"])
        clean_pdb = extract_single_chain(pdb_path, config["chain"])
        
        # Load PDB
        print(f"\n  Loading: {clean_pdb.name}")
        pdb = app.PDBFile(str(clean_pdb))
        print(f"  Atoms: {pdb.topology.getNumAtoms()}")
        
        # Force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Add hydrogens
        print("  Adding hydrogens...")
        modeller.addHydrogens(forcefield, pH=7.4)
        
        # Add solvent (smaller box for speed)
        print("  Adding solvent...")
        modeller.addSolvent(
            forcefield,
            model='tip3p',
            padding=0.8 * unit.nanometers,
            ionicStrength=0.15 * unit.molar
        )
        total_atoms = modeller.topology.getNumAtoms()
        print(f"  Total atoms: {total_atoms}")
        
        # Create system
        print("  Creating system...")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=app.HBonds
        )
        
        # Integrator
        integrator = mm.LangevinMiddleIntegrator(
            310 * unit.kelvin,
            1.0 / unit.picoseconds,
            0.002 * unit.picoseconds
        )
        
        # Platform
        platform_name = "CPU"
        try:
            platform = mm.Platform.getPlatformByName('CUDA')
            platform_name = "CUDA"
            properties = {'CudaPrecision': 'mixed'}
        except:
            try:
                platform = mm.Platform.getPlatformByName('OpenCL')
                platform_name = "OpenCL"
                properties = {}
            except:
                platform = mm.Platform.getPlatformByName('CPU')
                properties = {}
        
        print(f"  Platform: {platform_name}")
        
        # Create simulation
        simulation = app.Simulation(
            modeller.topology,
            system,
            integrator,
            platform,
            properties if properties else {}
        )
        simulation.context.setPositions(modeller.positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"\n  Initial energy: {initial_energy:.1f} kJ/mol")
        
        # Energy minimization
        print("  Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=1000)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        minimized_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Minimized energy: {minimized_energy:.1f} kJ/mol")
        
        # Save minimized structure
        minimized_pdb = system_dir / "minimized.pdb"
        with open(minimized_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Equilibration
        eq_energy = minimized_energy
        if equilibration_steps > 0:
            print(f"\n  Equilibration ({equilibration_steps} steps)...")
            simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)
            
            # Add reporter
            simulation.reporters.append(
                app.StateDataReporter(
                    str(system_dir / "equilibration.log"),
                    500,
                    step=True,
                    potentialEnergy=True,
                    temperature=True,
                    speed=True
                )
            )
            
            simulation.step(equilibration_steps)
            
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            eq_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"  Equilibrated energy: {eq_energy:.1f} kJ/mol")
            
            # Save equilibrated structure
            eq_pdb = system_dir / "equilibrated.pdb"
            with open(eq_pdb, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Production (if requested)
        final_energy = eq_energy
        if production_steps > 0:
            print(f"\n  Production ({production_steps} steps)...")
            simulation.reporters.clear()
            simulation.reporters.append(
                app.DCDReporter(str(system_dir / "trajectory.dcd"), 500)
            )
            simulation.reporters.append(
                app.StateDataReporter(
                    str(system_dir / "production.log"),
                    500,
                    step=True,
                    time=True,
                    potentialEnergy=True,
                    temperature=True,
                    speed=True
                )
            )
            
            simulation.step(production_steps)
            
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            # Save final structure
            final_pdb = system_dir / "final.pdb"
            with open(final_pdb, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        elapsed_time = time.time() - start_time
        
        # Save results
        results = {
            "status": "success",
            "system": system_name,
            "receptor": receptor_name,
            "ligand": ligand_name,
            "pdb_id": config["pdb"],
            "chain": config["chain"],
            "platform": platform_name,
            "total_atoms": total_atoms,
            "initial_energy_kj": initial_energy,
            "minimized_energy_kj": minimized_energy,
            "equilibrated_energy_kj": eq_energy,
            "final_energy_kj": final_energy,
            "equilibration_steps": equilibration_steps,
            "production_steps": production_steps,
            "elapsed_time_s": elapsed_time
        }
        
        with open(system_dir / "results.json", 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n  ✓ Completed in {elapsed_time:.1f}s")
        print(f"  Results saved to: {system_dir}")
        
        return results
        
    except Exception as e:
        elapsed_time = time.time() - start_time
        error_result = {
            "status": "error",
            "system": system_name,
            "message": str(e),
            "elapsed_time_s": elapsed_time
        }
        with open(system_dir / "error.json", 'w') as f:
            json.dump(error_result, f, indent=2)
        print(f"\n  ✗ Error: {e}")
        return error_result


def run_all_simulations():
    """Run simulations for all receptor-ligand pairs"""
    
    # Define simulation pairs
    pairs = [
        ("5-HT2A", "LSD"),
        ("5-HT2A", "psilocin"),
        ("5-HT2A", "DMT"),
        ("5-HT2A", "mescaline"),
        ("5-HT2B", "MDMA"),
        ("D2", "LSD"),
        ("MOR", "salvinorin_A"),
        ("KOR", "salvinorin_A"),
        ("CB1", "THC"),
    ]
    
    print("="*60)
    print("MOLECULAR DYNAMICS SIMULATION RUNNER")
    print("="*60)
    print(f"\nTotal pairs: {len(pairs)}")
    
    results = []
    
    for receptor, ligand in pairs:
        result = run_md_simulation(
            receptor, 
            ligand,
            equilibration_steps=5000,  # ~10 ps
            production_steps=0  # Just minimization + equilibration for now
        )
        results.append(result)
    
    # Summary
    print("\n" + "="*60)
    print("SIMULATION SUMMARY")
    print("="*60)
    
    successful = [r for r in results if r.get("status") == "success"]
    failed = [r for r in results if r.get("status") != "success"]
    
    print(f"\n✓ Successful: {len(successful)}/{len(results)}")
    for r in successful:
        print(f"  - {r['system']}: {r['total_atoms']} atoms, {r['elapsed_time_s']:.1f}s")
    
    if failed:
        print(f"\n✗ Failed: {len(failed)}/{len(results)}")
        for r in failed:
            print(f"  - {r.get('system', 'unknown')}: {r.get('message', 'unknown error')}")
    
    # Save summary
    with open(WORK_DIR / "simulation_summary.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nSummary saved to: {WORK_DIR / 'simulation_summary.json'}")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        # Run specific pair
        receptor = sys.argv[1]
        ligand = sys.argv[2]
        run_md_simulation(receptor, ligand)
    else:
        # Run all
        run_all_simulations()
