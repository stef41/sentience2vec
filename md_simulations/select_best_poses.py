#!/usr/bin/env python3
"""
Best Pose Selection Script
--------------------------
For each receptor-ligand pair, evaluate all 5 Chai-1 model predictions
and select the best pose based on binding pocket scoring.

Strategy:
1. For receptors with known binding sites (from PDB structures), score by
   proximity to key residues
2. For other receptors, use a generic ligand-protein contact score
3. Save best_pose.cif symlink/copy in each pair's chai1 directory
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from pathlib import Path
import json
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from gemmi import cif
import os

# Known binding site residues from literature/PDB structures
KNOWN_BINDING_SITES = {
    '5-HT2A': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE', 'TYR', 'ASN'],
        'critical_types': ['ASP', 'TRP', 'PHE'],  # Must have salt bridge + aromatic
        'source': 'PDB 6WHA'
    },
    '5-HT2B': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'PDB 5TUD'
    },
    '5-HT1A': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'Homology model'
    },
    '5-HT1B': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE'],
        'critical_types': ['ASP', 'TRP'],
        'source': 'PDB 4IAR'
    },
    '5-HT2C': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'PDB 6BQH'
    },
    'D1': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE', 'HIS'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'PDB 7CKZ'
    },
    'D2': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE', 'HIS'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'PDB 6CM4'
    },
    'D3': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE', 'HIS'],
        'critical_types': ['ASP', 'TRP', 'PHE'],
        'source': 'PDB 3PBL'
    },
    'D4': {
        'key_residues': ['ASP', 'SER', 'TRP', 'PHE'],
        'critical_types': ['ASP', 'TRP'],
        'source': 'PDB 5WIU'
    },
    'MOR': {
        'key_residues': ['ASP', 'TYR', 'TRP', 'HIS', 'MET'],
        'critical_types': ['ASP', 'TRP', 'TYR'],
        'source': 'PDB 4DKL'
    },
    'DOR': {
        'key_residues': ['ASP', 'TYR', 'TRP', 'HIS'],
        'critical_types': ['ASP', 'TRP', 'TYR'],
        'source': 'PDB 4N6H'
    },
    'KOR': {
        'key_residues': ['ASP', 'TYR', 'TRP', 'HIS'],
        'critical_types': ['ASP', 'TRP', 'TYR'],
        'source': 'PDB 4DJH'
    },
    'CB1': {
        'key_residues': ['PHE', 'TRP', 'SER', 'CYS', 'LEU'],
        'critical_types': ['TRP', 'PHE'],
        'source': 'PDB 5TGZ'
    },
    'CB2': {
        'key_residues': ['PHE', 'TRP', 'SER', 'CYS'],
        'critical_types': ['TRP', 'PHE'],
        'source': 'PDB 5ZTY'
    },
    'SERT': {
        'key_residues': ['ASP', 'TYR', 'ILE', 'PHE', 'SER'],
        'critical_types': ['ASP', 'PHE', 'TYR'],
        'source': 'PDB 5I6X'
    },
    'DAT': {
        'key_residues': ['ASP', 'PHE', 'SER', 'VAL'],
        'critical_types': ['ASP', 'PHE'],
        'source': 'PDB 4XP4'
    },
    'NET': {
        'key_residues': ['ASP', 'PHE', 'SER'],
        'critical_types': ['ASP', 'PHE'],
        'source': 'Homology model'
    },
    'NMDA': {
        'key_residues': ['ASN', 'GLU', 'TYR', 'PHE', 'TRP'],
        'critical_types': ['GLU', 'TYR', 'PHE'],
        'source': 'PDB 4PE5'
    },
    'GABAA': {
        'key_residues': ['PHE', 'TYR', 'SER', 'THR'],
        'critical_types': ['PHE', 'TYR'],
        'source': 'PDB 6HUP'
    },
    'GABAB': {
        'key_residues': ['TYR', 'SER', 'TRP', 'GLU'],
        'critical_types': ['TYR', 'TRP'],
        'source': 'PDB 7C7S'
    },
}


def parse_cif_with_gemmi(cif_path):
    """Parse CIF file to extract protein and ligand atoms using gemmi."""
    try:
        doc = cif.read(str(cif_path))
        block = doc[0]
        
        atom_site = block.find('_atom_site.',
            ['group_PDB', 'type_symbol', 'label_atom_id', 'label_comp_id',
             'label_seq_id', 'label_asym_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'])
        
        protein_atoms = []
        ligand_coords = []
        
        for row in atom_site:
            try:
                x, y, z = float(row[6]), float(row[7]), float(row[8])
                chain = row[5]
                res_name = row[3]
                res_num = row[4]
                
                if chain == 'B' or 'LIG' in res_name:
                    ligand_coords.append([x, y, z])
                else:
                    protein_atoms.append({
                        'residue': res_name,
                        'resnum': res_num,
                        'coords': np.array([x, y, z])
                    })
            except:
                continue
        
        return protein_atoms, np.array(ligand_coords) if ligand_coords else np.array([])
    except Exception as e:
        return [], np.array([])


def score_pose(cif_path, receptor_name):
    """
    Score a pose based on binding pocket residue proximity.
    
    Returns a score where higher is better:
    - Base score: Number of unique residue types within 6Å of ligand
    - Bonus for key residues if receptor is known
    - Bonus for critical residues (doubled)
    """
    protein_atoms, ligand_coords = parse_cif_with_gemmi(cif_path)
    
    if len(ligand_coords) == 0 or len(protein_atoms) == 0:
        return -1, {}
    
    ligand_center = ligand_coords.mean(axis=0)
    distance_cutoff = 6.0  # Angstroms
    
    # Find nearby residue types
    nearby_residue_types = set()
    nearby_count = 0
    
    for atom in protein_atoms:
        dist = np.linalg.norm(atom['coords'] - ligand_center)
        if dist < distance_cutoff:
            nearby_residue_types.add(atom['residue'])
            nearby_count += 1
    
    # Calculate score
    base_score = len(nearby_residue_types) * 10 + nearby_count * 0.1
    
    # Get receptor base name (e.g., "5-HT2A" from "5-HT2A_LSD")
    receptor_base = receptor_name.split('_')[0] if '_' in receptor_name else receptor_name
    
    # Bonus scoring for known receptors
    details = {
        'nearby_types': list(nearby_residue_types),
        'nearby_count': nearby_count,
        'base_score': base_score
    }
    
    if receptor_base in KNOWN_BINDING_SITES:
        site_info = KNOWN_BINDING_SITES[receptor_base]
        
        key_found = sum(1 for r in site_info['key_residues'] if r in nearby_residue_types)
        critical_found = sum(1 for r in site_info['critical_types'] if r in nearby_residue_types)
        
        # Bonus: +20 for each key residue, +50 for each critical residue
        key_bonus = key_found * 20
        critical_bonus = critical_found * 50
        
        total_score = base_score + key_bonus + critical_bonus
        
        details['key_residues_found'] = [r for r in site_info['key_residues'] if r in nearby_residue_types]
        details['critical_found'] = [r for r in site_info['critical_types'] if r in nearby_residue_types]
        details['key_bonus'] = key_bonus
        details['critical_bonus'] = critical_bonus
        details['receptor_matched'] = receptor_base
        
        return total_score, details
    
    return base_score, details


def process_pair(pair_dir):
    """Process a single pair directory and select the best pose."""
    pair_name = pair_dir.name
    chai1_dir = pair_dir / 'chai1'
    
    if not chai1_dir.exists():
        return None
    
    # Find all model files
    model_files = sorted(chai1_dir.glob('pred.model_idx_*.cif'))
    
    if not model_files:
        return None
    
    # Score each model
    best_score = -float('inf')
    best_model = None
    best_details = {}
    all_scores = {}
    
    for model_file in model_files:
        model_idx = model_file.stem.split('_')[-1]  # e.g., "0" from "pred.model_idx_0"
        score, details = score_pose(model_file, pair_name)
        all_scores[model_idx] = score
        
        if score > best_score:
            best_score = score
            best_model = model_file
            best_details = details
    
    if best_model is None:
        return None
    
    # Create symlink or copy to best_pose.cif
    best_pose_path = chai1_dir / 'best_pose.cif'
    try:
        if best_pose_path.exists() or best_pose_path.is_symlink():
            best_pose_path.unlink()
        
        # Use relative symlink
        best_pose_path.symlink_to(best_model.name)
    except Exception as e:
        # Fall back to copy
        try:
            shutil.copy2(best_model, best_pose_path)
        except:
            pass
    
    return {
        'pair': pair_name,
        'best_model': best_model.name,
        'best_score': best_score,
        'all_scores': all_scores,
        'details': best_details,
        'model_0_was_best': best_model.name == 'pred.model_idx_0.cif'
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Select best poses from Chai-1 predictions')
    parser.add_argument('--workers', type=int, default=32, help='Number of parallel workers')
    parser.add_argument('--sample', type=int, default=None, help='Process only N pairs (for testing)')
    args = parser.parse_args()
    
    base_dir = Path('all_pairs_results')
    
    print("="*80)
    print("BEST POSE SELECTION FROM CHAI-1 PREDICTIONS")
    print("="*80)
    print(f"Strategy: Score all 5 models per pair, select best binding pocket placement")
    print()
    
    # Get all pair directories
    pair_dirs = [d for d in sorted(base_dir.iterdir()) 
                 if d.is_dir() and not d.name.startswith('.') 
                 and d.name not in ('gifs', 'figures')]
    
    if args.sample:
        pair_dirs = pair_dirs[:args.sample]
    
    print(f"Processing {len(pair_dirs):,} receptor-ligand pairs...")
    print(f"Using {args.workers} parallel workers")
    print()
    
    results = []
    improved_count = 0
    total_processed = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_pair, pd): pd for pd in pair_dirs}
        
        with tqdm(total=len(futures), desc="Selecting best poses") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result:
                    results.append(result)
                    total_processed += 1
                    if not result['model_0_was_best']:
                        improved_count += 1
                pbar.update(1)
    
    # Summary statistics
    print()
    print("="*80)
    print("SELECTION SUMMARY")
    print("="*80)
    print(f"Total pairs processed: {total_processed:,}")
    print(f"Pairs where model_idx_0 was best: {total_processed - improved_count:,} ({(total_processed-improved_count)/total_processed*100:.1f}%)")
    print(f"Pairs improved by selecting alternative model: {improved_count:,} ({improved_count/total_processed*100:.1f}%)")
    print()
    
    # Model selection distribution
    model_counts = {}
    for r in results:
        model = r['best_model']
        model_counts[model] = model_counts.get(model, 0) + 1
    
    print("Best model distribution:")
    for model in sorted(model_counts.keys()):
        count = model_counts[model]
        print(f"  {model}: {count:,} ({count/total_processed*100:.1f}%)")
    
    # Show some examples of improved pairs
    print()
    print("-"*80)
    print("Examples of pairs improved by model selection:")
    print("-"*80)
    
    improved = [r for r in results if not r['model_0_was_best']]
    for r in improved[:10]:
        print(f"  {r['pair']}: {r['best_model']} (score: {r['best_score']:.0f})")
        if 'critical_found' in r['details']:
            print(f"    Critical residues: {r['details']['critical_found']}")
    
    if len(improved) > 10:
        print(f"  ... and {len(improved)-10:,} more")
    
    # Save results
    output_file = base_dir / 'best_pose_selection.json'
    with open(output_file, 'w') as f:
        json.dump({
            'summary': {
                'total_processed': total_processed,
                'improved_count': improved_count,
                'model_distribution': model_counts
            },
            'pairs': results
        }, f, indent=2)
    
    print()
    print(f"Results saved to: {output_file}")
    print()
    print("✓ Best pose selection complete!")
    print("  Each pair now has a 'best_pose.cif' symlink in its chai1 directory")


if __name__ == '__main__':
    main()
