#!/usr/bin/env python3
"""
Fix receptor PDB files by:
1. Extracting only the receptor chain (chain A typically)
2. Removing heteroatoms (ligands, waters)
3. Adding missing hydrogens using OpenMM Modeller
4. Energy minimizing
"""

import os
import sys
from pathlib import Path

def clean_pdb_extract_protein(input_pdb, output_pdb, target_chain='A'):
    """Extract only protein atoms from target chain, remove heteroatoms"""
    
    STANDARD_AA = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX',
        'MSE', 'SEC'  # Modified residues
    }
    
    lines = []
    atom_num = 1
    
    with open(input_pdb) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                res_name = line[17:20].strip()
                
                # Only keep protein residues from target chain
                if chain == target_chain and res_name in STANDARD_AA:
                    # Convert MSE to MET
                    if res_name == 'MSE':
                        line = line[:17] + 'MET' + line[20:]
                        # SE -> SD for selenium
                        if line[12:16].strip() == 'SE':
                            line = line[:12] + ' SD ' + line[16:]
                    
                    # Renumber atoms
                    line = f'ATOM  {atom_num:5d}' + line[11:]
                    lines.append(line)
                    atom_num += 1
    
    lines.append('TER\n')
    lines.append('END\n')
    
    with open(output_pdb, 'w') as f:
        f.writelines(lines)
    
    return len(lines) - 2

def prepare_receptor_openmm(input_pdb, output_pdb):
    """Use OpenMM to add hydrogens and minimize"""
    
    from openmm.app import PDBFile, ForceField, Modeller, NoCutoff, HBonds
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picoseconds
    import openmm
    
    # Load PDB
    pdb = PDBFile(str(input_pdb))
    
    # Force field
    forcefield = ForceField('amber14-all.xml')
    
    # Create modeller and add hydrogens
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    
    # Create system (vacuum, no cutoff for speed)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    
    # Integrator
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 0.002*picoseconds)
    
    # Platform
    try:
        platform = Platform.getPlatformByName('OpenCL')
    except:
        platform = Platform.getPlatformByName('CPU')
    
    # Simulation
    simulation = openmm.app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize
    simulation.minimizeEnergy(maxIterations=500)
    
    # Save
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    return True

def main():
    base_dir = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/receptors_full")
    
    # Get chain info from database
    sys.path.insert(0, str(base_dir.parent))
    from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE
    
    success = 0
    failed = []
    
    receptors = list(base_dir.iterdir())
    print(f"Processing {len(receptors)} receptors...")
    
    for receptor_dir in sorted(receptors):
        if not receptor_dir.is_dir():
            continue
            
        name = receptor_dir.name
        raw_pdb = receptor_dir / f"{name}_raw.pdb"
        clean_pdb = receptor_dir / f"{name}_clean.pdb"
        final_pdb = receptor_dir / f"{name}_prepared.pdb"
        
        if not raw_pdb.exists():
            continue
        
        # Skip if already prepared
        if final_pdb.exists():
            print(f"  {name}: already prepared")
            success += 1
            continue
        
        # Get chain from database
        chain = 'A'
        if name in RECEPTOR_PDB_DATABASE:
            chain = RECEPTOR_PDB_DATABASE[name].get('chain', 'A') or 'A'
        
        print(f"  {name} (chain {chain})...", end=" ", flush=True)
        
        try:
            # Step 1: Clean PDB
            n_atoms = clean_pdb_extract_protein(raw_pdb, clean_pdb, chain)
            if n_atoms < 100:
                print(f"only {n_atoms} atoms")
                failed.append(name)
                continue
            
            # Step 2: Add hydrogens and minimize
            prepare_receptor_openmm(clean_pdb, final_pdb)
            print(f"✓ ({n_atoms} heavy atoms)")
            success += 1
            
        except Exception as e:
            print(f"✗ ({str(e)[:60]})")
            failed.append(name)
    
    print(f"\n{'='*60}")
    print(f"Success: {success}, Failed: {len(failed)}")
    if failed:
        print(f"Failed: {failed[:10]}...")

if __name__ == "__main__":
    main()
