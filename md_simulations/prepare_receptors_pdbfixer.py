#!/usr/bin/env python3
"""
Prepare ALL receptors using PDBFixer to handle missing atoms.
"""

import sys
from pathlib import Path

def prepare_receptor(input_pdb, output_pdb, target_chain='A'):
    """Use PDBFixer to clean and prepare receptor"""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, ForceField, Modeller, NoCutoff, HBonds
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.unit import kelvin, picoseconds
    import openmm
    
    # Load and fix with PDBFixer
    fixer = PDBFixer(filename=str(input_pdb))
    
    # Keep only target chain
    chains_to_remove = []
    for chain in fixer.topology.chains():
        if chain.id != target_chain:
            chains_to_remove.append(chain.id)
    fixer.removeChains(chainIds=chains_to_remove)
    
    # Find and fix missing residues/atoms
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    # Save fixed structure
    fixed_pdb = output_pdb.parent / f"{output_pdb.stem}_fixed.pdb"
    with open(fixed_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    # Now minimize with OpenMM
    pdb = PDBFile(str(fixed_pdb))
    forcefield = ForceField('amber14-all.xml')
    
    modeller = Modeller(pdb.topology, pdb.positions)
    
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 0.002*picoseconds)
    
    try:
        platform = Platform.getPlatformByName('OpenCL')
    except:
        platform = Platform.getPlatformByName('CPU')
    
    simulation = openmm.app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=500)
    
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    return True

def main():
    base_dir = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/receptors_full")
    
    sys.path.insert(0, str(base_dir.parent))
    from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE
    
    success, failed = 0, []
    
    receptors = sorted([d for d in base_dir.iterdir() if d.is_dir()])
    print(f"Processing {len(receptors)} receptors with PDBFixer...")
    print("="*60)
    
    for i, receptor_dir in enumerate(receptors, 1):
        name = receptor_dir.name
        raw_pdb = receptor_dir / f"{name}_raw.pdb"
        final_pdb = receptor_dir / f"{name}_prepared.pdb"
        
        if not raw_pdb.exists():
            continue
        
        if final_pdb.exists():
            print(f"  [{i}/{len(receptors)}] {name}... (cached)")
            success += 1
            continue
        
        chain = 'A'
        if name in RECEPTOR_PDB_DATABASE:
            chain = RECEPTOR_PDB_DATABASE[name].get('chain', 'A') or 'A'
        
        print(f"  [{i}/{len(receptors)}] {name} (chain {chain})...", end=" ", flush=True)
        
        try:
            prepare_receptor(raw_pdb, final_pdb, chain)
            print("✓")
            success += 1
        except Exception as e:
            print(f"✗ ({str(e)[:50]})")
            failed.append(name)
    
    print("="*60)
    print(f"Success: {success}, Failed: {len(failed)}")
    if failed:
        print(f"Failed: {failed}")

if __name__ == "__main__":
    main()
