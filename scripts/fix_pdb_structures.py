#!/usr/bin/env python3
"""
PDB Structure Cleaner and Fixer
================================
Fixes common issues with experimental PDB structures:
- Non-standard residue names
- Missing atoms
- Alternate conformations
- Incomplete residues
"""

import os
import sys
from pathlib import Path

def clean_pdb_for_openmm(input_pdb: str, output_pdb: str = None) -> str:
    """
    Clean PDB file for use with OpenMM/Amber force field.
    
    Fixes:
    - Removes heterogens (ligands, waters, ions)
    - Standardizes histidine names (HIS -> HID/HIE/HIP)
    - Removes alternate conformations (keeps A)
    - Removes incomplete residues
    - Renumbers atoms
    """
    if output_pdb is None:
        output_pdb = input_pdb.replace('.pdb', '_fixed.pdb')
    
    print(f"Cleaning PDB: {input_pdb}")
    
    # Standard amino acid names
    STANDARD_AA = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
        'THR', 'TRP', 'TYR', 'VAL',
        # Protonation variants
        'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN'
    }
    
    # Required atoms for each amino acid (backbone + key side chain)
    REQUIRED_BACKBONE = {'N', 'CA', 'C', 'O'}
    
    # Required side chain atoms for residues that need them
    REQUIRED_SIDECHAIN = {
        'HIS': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HID': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIE': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIP': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'PHE': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TYR': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
        'TRP': {'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'ARG': {'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
        'LYS': {'CB', 'CG', 'CD', 'CE', 'NZ'},
        'ASP': {'CB', 'CG', 'OD1', 'OD2'},
        'GLU': {'CB', 'CG', 'CD', 'OE1', 'OE2'},
        'ASN': {'CB', 'CG', 'OD1', 'ND2'},
        'GLN': {'CB', 'CG', 'CD', 'OE1', 'NE2'},
    }
    
    # Read input
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    # First pass: collect residue information
    residues = {}  # (chain, resnum) -> list of atoms
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            atom_name = line[12:16].strip()
            alt_loc = line[16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = line[22:26].strip()
            
            # Skip alternate conformations (keep A or first one)
            if alt_loc and alt_loc != 'A':
                continue
            
            # Skip non-protein
            if res_name not in STANDARD_AA:
                continue
            
            key = (chain, res_num)
            if key not in residues:
                residues[key] = {'name': res_name, 'atoms': set()}
            residues[key]['atoms'].add(atom_name)
    
    # Identify complete residues (have backbone atoms AND required sidechain)
    complete_residues = set()
    for key, data in residues.items():
        # Must have backbone
        if not REQUIRED_BACKBONE.issubset(data['atoms']):
            continue
        
        # Check side chain requirements
        res_name = data['name']
        if res_name in REQUIRED_SIDECHAIN:
            required = REQUIRED_SIDECHAIN[res_name]
            if not required.issubset(data['atoms']):
                print(f"  Skipping incomplete {res_name} at {key}: missing {required - data['atoms']}")
                continue
        
        complete_residues.add(key)
    
    print(f"  Total residues: {len(residues)}")
    print(f"  Complete residues: {len(complete_residues)}")
    
    # Second pass: write cleaned PDB
    output_lines = []
    atom_num = 0
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            alt_loc = line[16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = line[22:26].strip()
            
            key = (chain, res_num)
            
            # Skip alternate conformations
            if alt_loc and alt_loc != 'A':
                continue
            
            # Skip non-protein
            if res_name not in STANDARD_AA:
                continue
            
            # Skip incomplete residues
            if key not in complete_residues:
                continue
            
            atom_num += 1
            
            # Fix histidine naming - default to HID (delta protonated)
            if res_name == 'HIS':
                res_name = 'HID'
            
            # Reconstruct line with fixes
            # Remove alternate location indicator
            new_line = (
                f"ATOM  {atom_num:5d} {line[12:16]}"  # Atom name (keep original spacing)
                f" "  # Remove alt loc
                f"{res_name:3s} "
                f"{chain}"
                f"{res_num:>4s}"
                f"    "  # Insertion code
                f"{line[30:54]}"  # Coordinates
                f"{line[54:60]}"  # Occupancy
                f"{line[60:66]}"  # B-factor
                f"          "
                f"{line[76:78] if len(line) > 76 else '  '}"  # Element
                f"\n"
            )
            output_lines.append(new_line)
        
        elif line.startswith('TER'):
            output_lines.append(f"TER   {atom_num+1:5d}\n")
        
        elif line.startswith('END'):
            output_lines.append('END\n')
    
    # Write output
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    print(f"  Output: {output_pdb}")
    print(f"  Atoms written: {atom_num}")
    
    return output_pdb


def main():
    """Clean all PDB structures in the structures directory"""
    structures_dir = Path(__file__).parent.parent / "md_simulations" / "structures"
    
    if not structures_dir.exists():
        print(f"Structures directory not found: {structures_dir}")
        return
    
    print("="*60)
    print("PDB STRUCTURE CLEANING")
    print("="*60)
    
    for pdb_file in structures_dir.glob("*.pdb"):
        # Skip already cleaned files
        if "_fixed" in pdb_file.name or "_cleaned" in pdb_file.name:
            continue
        
        output_file = pdb_file.with_name(pdb_file.stem + "_fixed.pdb")
        
        try:
            clean_pdb_for_openmm(str(pdb_file), str(output_file))
            print()
        except Exception as e:
            print(f"  Error: {e}\n")


if __name__ == "__main__":
    main()
