#!/usr/bin/env python3
"""
Validate Chai-1 poses against published binding site data from PDB structures.
Compares ligand positions with known binding pocket residues from cryo-EM/X-ray structures.

Key references:
- Kim et al. (2020) Cell 182:1574-1588 - 5-HT2A + 25-CN-NBOH (PDB: 6WHA)
- Wacker et al. (2017) Cell 168:377-389 - 5-HT2B + ergotamine (PDB: 5TUD)
- Wang et al. (2013) Science - D2 receptor structure
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from pathlib import Path
import json

# Published binding site residues from literature
# Format: residue_number (Ballesteros-Weinstein numbering where available)
KNOWN_BINDING_SITES = {
    '5-HT2A': {
        'orthosteric': [
            # From Kim et al. 2020 (6WHA) - key residues in orthosteric pocket
            ('D155', 'Asp', 155, 'Salt bridge with amine'),  # 3.32
            ('S159', 'Ser', 159, 'H-bond'),
            ('T160', 'Thr', 160, 'H-bond'),
            ('S242', 'Ser', 242, 'H-bond donor'),  # 5.43
            ('W336', 'Trp', 336, 'Aromatic stacking'),  # 6.48 - toggle switch
            ('F339', 'Phe', 339, 'Hydrophobic'),  # 6.51
            ('F340', 'Phe', 340, 'Aromatic cage'),  # 6.52
            ('N343', 'Asn', 343, 'H-bond'),  # 6.55
            ('Y370', 'Tyr', 370, 'Aromatic stacking'),  # 7.43
        ],
        'extended_pocket': [
            ('V156', 'Val', 156),
            ('I163', 'Ile', 163),
            ('L229', 'Leu', 229),
            ('G238', 'Gly', 238),
            ('L362', 'Leu', 362),
            ('V366', 'Val', 366),
        ],
        'source': 'PDB 6WHA (Kim et al. 2020 Cell)',
    },
    '5-HT2B': {
        'orthosteric': [
            ('D135', 'Asp', 135, 'Salt bridge'),  # 3.32
            ('S139', 'Ser', 139, 'H-bond'),
            ('S222', 'Ser', 222, 'H-bond'),  # 5.43
            ('W337', 'Trp', 337, 'Aromatic'),  # 6.48
            ('F340', 'Phe', 340, 'Hydrophobic'),  # 6.51
            ('F341', 'Phe', 341, 'Aromatic cage'),  # 6.52
        ],
        'source': 'PDB 5TUD (Wacker et al. 2017 Cell)',
    },
    '5-HT1A': {
        'orthosteric': [
            ('D116', 'Asp', 116, 'Salt bridge'),  # 3.32
            ('S199', 'Ser', 199, 'H-bond'),  # 5.43
            ('W358', 'Trp', 358, 'Aromatic'),  # 6.48
            ('F361', 'Phe', 361, 'Aromatic'),  # 6.51
            ('F362', 'Phe', 362, 'Aromatic'),  # 6.52
        ],
        'source': 'Homology model + mutagenesis data',
    },
    'D2': {
        'orthosteric': [
            ('D114', 'Asp', 114, 'Salt bridge'),  # 3.32
            ('S193', 'Ser', 193, 'H-bond'),  # 5.43
            ('S194', 'Ser', 194, 'H-bond'),  # 5.44
            ('W386', 'Trp', 386, 'Aromatic'),  # 6.48
            ('F389', 'Phe', 389, 'Hydrophobic'),  # 6.51
            ('F390', 'Phe', 390, 'Aromatic'),  # 6.52
            ('H393', 'His', 393, 'H-bond'),  # 6.55
        ],
        'source': 'PDB 6CM4 (Wang et al. 2018 Nature)',
    },
    'MOR': {
        'orthosteric': [
            ('D147', 'Asp', 147, 'Salt bridge'),  # 3.32
            ('Y148', 'Tyr', 148, 'Aromatic'),
            ('W293', 'Trp', 293, 'Aromatic'),  # 6.48
            ('H297', 'His', 297, 'H-bond'),  # 6.52
            ('Y326', 'Tyr', 326, 'Aromatic'),  # 7.43
        ],
        'source': 'PDB 5C1M (Huang et al. 2015 Nature)',
    },
    'CB1': {
        'orthosteric': [
            ('F170', 'Phe', 170, 'Aromatic'),
            ('F174', 'Phe', 174, 'Aromatic'),
            ('F177', 'Phe', 177, 'Aromatic'),
            ('F189', 'Phe', 189, 'Hydrophobic'),
            ('W279', 'Trp', 279, 'Aromatic'),
            ('W356', 'Trp', 356, 'Aromatic'),  # 6.48
        ],
        'source': 'PDB 5TGZ (Hua et al. 2017 Nature)',
    },
    'SERT': {
        'orthosteric': [
            ('D98', 'Asp', 98, 'Salt bridge'),
            ('Y95', 'Tyr', 95, 'Aromatic'),
            ('I172', 'Ile', 172, 'Hydrophobic'),
            ('F335', 'Phe', 335, 'Aromatic'),
            ('S336', 'Ser', 336, 'H-bond'),
            ('F341', 'Phe', 341, 'Aromatic'),
        ],
        'source': 'PDB 5I6X (Coleman et al. 2016 Nature)',
    },
}

def parse_cif_atoms(cif_file):
    """Parse CIF and return protein + ligand atoms separately."""
    protein_atoms = []
    ligand_atoms = []
    
    with open(cif_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) >= 18:
                    try:
                        x = float(parts[10])
                        y = float(parts[11])
                        z = float(parts[12])
                        residue = parts[5]
                        resnum = parts[6] if parts[6].isdigit() else parts[7]
                        chain = parts[17]
                        element = parts[2]
                        
                        atom = {
                            'x': x, 'y': y, 'z': z,
                            'residue': residue,
                            'resnum': int(resnum) if resnum.isdigit() else 0,
                            'chain': chain,
                            'element': element
                        }
                        
                        if 'LIG' in residue or chain == 'B':
                            ligand_atoms.append(atom)
                        else:
                            protein_atoms.append(atom)
                    except:
                        continue
    
    return protein_atoms, ligand_atoms

def get_ligand_center(ligand_atoms):
    """Get center of mass of ligand."""
    if not ligand_atoms:
        return None
    coords = np.array([[a['x'], a['y'], a['z']] for a in ligand_atoms])
    return coords.mean(axis=0)

def get_binding_site_residues(protein_atoms, ligand_center, cutoff=8.0):
    """Find protein residues within cutoff distance of ligand center."""
    if ligand_center is None:
        return {}
    
    residues = {}
    for atom in protein_atoms:
        coord = np.array([atom['x'], atom['y'], atom['z']])
        dist = np.linalg.norm(coord - ligand_center)
        if dist < cutoff:
            key = (atom['residue'], atom['resnum'])
            if key not in residues or dist < residues[key]['min_dist']:
                residues[key] = {
                    'residue': atom['residue'],
                    'resnum': atom['resnum'],
                    'min_dist': dist
                }
    
    return residues

def validate_pose(pair_name, cif_file, known_site):
    """Validate a predicted pose against known binding site."""
    protein_atoms, ligand_atoms = parse_cif_atoms(cif_file)
    
    if not ligand_atoms:
        return {'valid': False, 'error': 'No ligand found'}
    
    ligand_center = get_ligand_center(ligand_atoms)
    nearby_residues = get_binding_site_residues(protein_atoms, ligand_center, cutoff=10.0)
    
    # Check overlap with known binding site residues
    known_resnums = [r[2] for r in known_site.get('orthosteric', [])]
    found_known = []
    missing_known = []
    
    nearby_resnums = [r['resnum'] for r in nearby_residues.values()]
    
    for resinfo in known_site.get('orthosteric', []):
        resnum = resinfo[2]
        resname = resinfo[1]
        # Allow some flexibility in numbering (+/- 5)
        found = any(abs(resnum - nr) <= 5 for nr in nearby_resnums)
        if found:
            found_known.append(f"{resname}{resnum}")
        else:
            missing_known.append(f"{resname}{resnum}")
    
    # Calculate overlap score
    total_known = len(known_site.get('orthosteric', []))
    overlap_score = len(found_known) / total_known if total_known > 0 else 0
    
    return {
        'valid': overlap_score >= 0.5,
        'overlap_score': overlap_score,
        'found_residues': found_known,
        'missing_residues': missing_known,
        'ligand_atoms': len(ligand_atoms),
        'nearby_residues': len(nearby_residues),
        'ligand_center': ligand_center.tolist() if ligand_center is not None else None,
    }

def main():
    output_base = Path('all_pairs_results')
    
    print("="*80)
    print("POSE VALIDATION AGAINST PUBLISHED BINDING SITE DATA")
    print("="*80)
    print("Using best_pose.cif (selected from all 5 Chai-1 models)")
    print()
    
    # Key pairs to validate
    validation_pairs = [
        ('5-HT2A_DMT', '5-HT2A'),
        ('5-HT2A_LSD', '5-HT2A'),
        ('5-HT2A_Psilocin', '5-HT2A'),
        ('5-HT2A_Mescaline', '5-HT2A'),
        ('5-HT2A_2C-B', '5-HT2A'),
        ('5-HT2A_MDMA', '5-HT2A'),
        ('5-HT2A_5-MeO-DMT', '5-HT2A'),
        ('5-HT2B_DMT', '5-HT2B'),
        ('5-HT2B_LSD', '5-HT2B'),
        ('5-HT1A_DMT', '5-HT1A'),
        ('5-HT1A_LSD', '5-HT1A'),
        ('D2_LSD', 'D2'),
        ('D2_DMT', 'D2'),
        ('D2_Dopamine', 'D2'),
        ('MOR_Morphine', 'MOR'),
        ('MOR_Fentanyl', 'MOR'),
        ('CB1_THC-O-acetate', 'CB1'),
        ('SERT_MDMA', 'SERT'),
        ('SERT_Serotonin', 'SERT'),
    ]
    
    results = {}
    valid_count = 0
    total_count = 0
    
    for pair_name, receptor in validation_pairs:
        # Use best_pose.cif if available, otherwise fall back to model_idx_0
        best_pose = output_base / pair_name / 'chai1' / 'best_pose.cif'
        cif_file = best_pose if best_pose.exists() else output_base / pair_name / 'chai1' / 'pred.model_idx_0.cif'
        
        if not cif_file.exists():
            print(f"✗ {pair_name}: CIF not found")
            continue
        
        if receptor not in KNOWN_BINDING_SITES:
            print(f"? {pair_name}: No reference data for {receptor}")
            continue
        
        known_site = KNOWN_BINDING_SITES[receptor]
        result = validate_pose(pair_name, cif_file, known_site)
        results[pair_name] = result
        total_count += 1
        
        if result['valid']:
            valid_count += 1
            status = "✓ VALID"
        else:
            status = "✗ INVALID"
        
        print(f"\n{status} - {pair_name}")
        print(f"  Overlap score: {result['overlap_score']:.1%}")
        print(f"  Found key residues: {', '.join(result['found_residues'][:5])}...")
        if result['missing_residues']:
            print(f"  Missing: {', '.join(result['missing_residues'][:3])}...")
        print(f"  Ligand atoms: {result['ligand_atoms']}, Nearby protein residues: {result['nearby_residues']}")
    
    # Summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    print(f"\nTotal validated: {total_count}")
    print(f"Valid poses: {valid_count} ({valid_count/total_count*100:.0f}%)")
    print(f"Invalid poses: {total_count - valid_count}")
    
    print("\n" + "-"*80)
    print("REFERENCE DATA SOURCES:")
    print("-"*80)
    for receptor, data in KNOWN_BINDING_SITES.items():
        print(f"  {receptor}: {data['source']}")
    
    # Save results
    with open(output_base / 'pose_validation_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {output_base / 'pose_validation_results.json'}")

if __name__ == '__main__':
    main()
