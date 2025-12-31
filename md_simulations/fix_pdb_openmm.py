#!/usr/bin/env python3
"""
Custom PDB Fixer for OpenMM Simulations
========================================
Handles common issues with experimental PDB structures:
- Missing terminal atoms
- Incomplete residues
- Non-standard residue names
"""
import openmm.app as app
from pathlib import Path
import numpy as np


def fix_pdb_for_openmm(input_pdb: Path, output_pdb: Path = None) -> Path:
    """
    Fix a PDB file for use with OpenMM's Amber force field.
    
    Issues addressed:
    1. Remove incomplete terminal residues
    2. Fix histidine naming
    3. Remove alternate conformations
    4. Cap chain termini if needed
    """
    if output_pdb is None:
        output_pdb = input_pdb.with_suffix('.fixed.pdb')
    
    STANDARD_AA = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
        'HIS', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 
        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    # Heavy atoms required for each residue type
    REQUIRED_HEAVY_ATOMS = {
        'ALA': {'N', 'CA', 'C', 'O', 'CB'},
        'GLY': {'N', 'CA', 'C', 'O'},
        'VAL': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'},
        'LEU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'},
        'ILE': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'},
        'PRO': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'},
        'PHE': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TYR': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
        'TRP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'SER': {'N', 'CA', 'C', 'O', 'CB', 'OG'},
        'THR': {'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'},
        'CYS': {'N', 'CA', 'C', 'O', 'CB', 'SG'},
        'MET': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'},
        'ASN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'},
        'GLN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'},
        'ASP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'},
        'GLU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
        'LYS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'},
        'ARG': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
        'HIS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HID': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIE': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    }
    
    # Read input
    with open(input_pdb) as f:
        lines = f.readlines()
    
    # Parse residue information
    residues = {}  # (chain, resnum) -> {'name': str, 'atoms': set}
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            alt_loc = line[16]
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = line[22:27].strip()
            
            # Skip alternate conformations
            if alt_loc not in [' ', 'A']:
                continue
            
            # Skip non-standard residues
            if res_name not in STANDARD_AA:
                continue
            
            key = (chain, res_num)
            if key not in residues:
                residues[key] = {'name': res_name, 'atoms': set(), 'lines': []}
            residues[key]['atoms'].add(atom_name)
            residues[key]['lines'].append(line)
    
    # Find complete residues
    complete_residues = set()
    for key, data in residues.items():
        res_name = data['name']
        if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
            required = REQUIRED_HEAVY_ATOMS.get('HIS', set())
        else:
            required = REQUIRED_HEAVY_ATOMS.get(res_name, set())
        
        # Check if backbone is complete
        backbone = {'N', 'CA', 'C', 'O'}
        heavy_atoms = {a for a in data['atoms'] if not a.startswith('H')}
        
        # Allow some missing side chain atoms, but backbone must be complete
        if backbone.issubset(heavy_atoms):
            # Check that most side chain atoms are present
            missing = required - heavy_atoms
            if len(missing) <= 2:  # Allow up to 2 missing atoms
                complete_residues.add(key)
    
    # Get chain order
    chains = {}
    for key in sorted(residues.keys()):
        chain = key[0]
        if chain not in chains:
            chains[chain] = []
        chains[chain].append(key)
    
    # For each chain, find the longest contiguous stretch
    output_residues = []
    for chain, keys in chains.items():
        # Sort by residue number
        keys_sorted = sorted(keys, key=lambda k: int(k[1].replace('A', '').replace('B', '')))
        
        # Find complete residues in this chain
        complete_in_chain = [k for k in keys_sorted if k in complete_residues]
        
        if len(complete_in_chain) < 10:  # Need at least 10 residues
            continue
        
        # Skip first and last few residues (often incomplete)
        if len(complete_in_chain) > 10:
            complete_in_chain = complete_in_chain[2:-2]
        
        output_residues.extend(complete_in_chain)
    
    print(f"  Keeping {len(output_residues)} residues from {len(residues)} total")
    
    # Write output
    output_lines = []
    atom_num = 0
    
    for key in output_residues:
        data = residues[key]
        for line in data['lines']:
            atom_name = line[12:16].strip()
            res_name = data['name']
            
            # Fix histidine naming (default to delta protonated)
            if res_name == 'HIS':
                res_name = 'HID'
            
            atom_num += 1
            
            # Reconstruct line
            new_line = (
                f"ATOM  {atom_num:5d} "
                f"{line[12:16]}"  # Atom name
                f" "  # Alt loc (cleared)
                f"{res_name:3s} "
                f"{line[21]}"  # Chain
                f"{line[22:27]}"  # Residue number
                f"   "
                f"{line[30:54]}"  # Coordinates
                f"{line[54:60]}"  # Occupancy
                f"{line[60:66]}"  # B-factor
                f"          "
                f"{line[76:78] if len(line) > 76 else '  '}"  # Element
                f"\n"
            )
            output_lines.append(new_line)
    
    output_lines.append("TER\nEND\n")
    
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    print(f"  Output: {output_pdb} ({atom_num} atoms)")
    return output_pdb


def main():
    """Fix all PDB structures"""
    structures_dir = Path(__file__).parent / "structures"
    
    print("="*60)
    print("FIXING PDB STRUCTURES FOR OPENMM")
    print("="*60)
    
    # Original PDB files to fix
    pdb_files = [
        "6WHA.pdb",  # 5-HT2A
        "5TVN.pdb",  # 5-HT2B
        "6CM4.pdb",  # D2
        "5C1M.pdb",  # MOR
        "4DJH.pdb",  # KOR
        "5TGZ.pdb",  # CB1
    ]
    
    for pdb_name in pdb_files:
        pdb_path = structures_dir / pdb_name
        if not pdb_path.exists():
            print(f"\n{pdb_name}: Not found")
            continue
        
        print(f"\n{pdb_name}:")
        output_path = structures_dir / pdb_name.replace('.pdb', '_openmm.pdb')
        fix_pdb_for_openmm(pdb_path, output_path)


if __name__ == "__main__":
    main()
