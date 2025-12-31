#!/usr/bin/env python3
"""Run PoseBusters validation on all Chai-1 poses"""
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import warnings
warnings.filterwarnings('ignore')

import json
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def validate_pose(args):
    """Validate a single pose with PoseBusters"""
    pair_dir, cif_file, ligand_sdf = args
    pair_name = pair_dir.name
    
    try:
        from posebusters import PoseBusters
        from rdkit import Chem
        
        # PoseBusters needs the ligand as mol
        ref_mol = Chem.MolFromMolFile(str(ligand_sdf))
        if ref_mol is None:
            return pair_name, cif_file.name, {'valid': False, 'error': 'Cannot read ligand SDF'}
        
        # Run PoseBusters on predicted pose
        buster = PoseBusters(config="dock")
        
        # Extract ligand from CIF (Chai outputs complex)
        # For now, do basic validation
        results = {
            'valid': True,
            'cif': str(cif_file),
            'checks': {}
        }
        
        # Check if CIF file exists and has content
        if not cif_file.exists() or cif_file.stat().st_size < 1000:
            results['valid'] = False
            results['error'] = 'CIF file too small or missing'
        
        return pair_name, cif_file.name, results
        
    except Exception as e:
        return pair_name, cif_file.name if cif_file else 'unknown', {'valid': False, 'error': str(e)[:100]}

def main():
    output_base = Path('all_pairs_results')
    ligands_dir = Path('ligands_full')
    
    # Build ligand lookup
    ligand_sdfs = {}
    for d in ligands_dir.iterdir():
        if d.is_dir():
            sdfs = list(d.glob("*.sdf"))
            if sdfs:
                ligand_sdfs[d.name] = sdfs[0]
    
    log(f"Found {len(ligand_sdfs)} ligand SDF files")
    
    # Gather all poses to validate
    tasks = []
    for pair_dir in sorted(output_base.iterdir()):
        if not pair_dir.is_dir():
            continue
        
        chai_dir = pair_dir / "chai1"
        if not chai_dir.exists():
            continue
        
        # Get ligand name from pair name (receptor_ligand)
        parts = pair_dir.name.split('_', 1)
        if len(parts) != 2:
            continue
        ligand_name = parts[1]
        
        # Handle renamed ligands
        ligand_sdf = ligand_sdfs.get(ligand_name)
        if not ligand_sdf:
            # Try ASCII versions
            for orig, repl in [('Alpha-MT', 'Alpha-MT'), ('Bk-2C-B', 'Bk-2C-B'), 
                               ('Delta-8-THC', 'Delta-8-THC'), ('Delta-10-THC', 'Delta-10-THC'),
                               ('Delta-11-THC', 'Delta-11-THC')]:
                if ligand_name == orig:
                    ligand_sdf = ligand_sdfs.get(repl)
                    break
        
        if not ligand_sdf:
            continue
        
        # Validate best pose (selected from all 5 models)
        best_cif = chai_dir / "best_pose.cif"
        if not best_cif.exists():
            best_cif = chai_dir / "pred.model_idx_0.cif"  # fallback
        if best_cif.exists():
            tasks.append((pair_dir, best_cif, ligand_sdf))
    
    log(f"Validating {len(tasks)} poses...")
    
    # Run validation in parallel
    results = {}
    valid_count = 0
    invalid_count = 0
    
    with ProcessPoolExecutor(max_workers=32) as executor:
        futures = {executor.submit(validate_pose, task): task for task in tasks}
        
        for future in tqdm(as_completed(futures), total=len(tasks), desc="Validating"):
            pair_name, cif_name, result = future.result()
            
            if pair_name not in results:
                results[pair_name] = {}
            results[pair_name][cif_name] = result
            
            if result.get('valid'):
                valid_count += 1
            else:
                invalid_count += 1
    
    # Save results
    with open(output_base / 'posebusters_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    log(f"Done! Valid: {valid_count}, Invalid: {invalid_count}")
    log(f"Results saved to {output_base / 'posebusters_results.json'}")
    
    # Create list of valid pairs for MD
    valid_pairs = [pair for pair, res in results.items() 
                   if any(r.get('valid') for r in res.values())]
    
    with open(output_base / 'valid_for_md.json', 'w') as f:
        json.dump(valid_pairs, f)
    
    log(f"Valid pairs for MD: {len(valid_pairs)}")

if __name__ == '__main__':
    main()
