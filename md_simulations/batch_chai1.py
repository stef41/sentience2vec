#!/usr/bin/env python3
"""
Robust batch Chai-1 for all receptor-ligand pairs.
Run with: nohup python batch_chai1.py > batch.log 2>&1 &
"""
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['TOKENIZERS_PARALLELISM'] = 'false'

import json
import time
import sys
from pathlib import Path
from datetime import datetime

def log(msg):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)

def extract_sequence(pdb_path):
    aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
          'GLY':'G','HIS':'H','HID':'H','HIE':'H','ILE':'I','LEU':'L','LYS':'K',
          'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    seq, seen = [], set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith('ATOM'):
                    rn, ri = line[17:20].strip(), line[22:26].strip()
                    if ri not in seen and rn in aa:
                        seq.append(aa[rn])
                        seen.add(ri)
    except:
        return None
    return ''.join(seq) if seq else None

def extract_smiles(sdf_path):
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_path))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

def main():
    log("="*70)
    log("BATCH CHAI-1: ALL RECEPTOR-LIGAND PAIRS")
    log("="*70)
    
    # Patch timeout before importing chai_lab
    try:
        import chai_lab.utils.timeout as timeout_module
        def no_timeout(seconds):
            def decorator(func):
                return func
            return decorator
        timeout_module.timeout = no_timeout
        log("✓ Patched chai_lab timeout")
    except Exception as e:
        log(f"⚠ Could not patch timeout: {e}")
    
    # Load pairs
    pairs_file = Path('all_pairs_results/pairs_list.json')
    with open(pairs_file) as f:
        all_pairs = json.load(f)
    log(f"Total pairs: {len(all_pairs):,}")
    
    output_base = Path('all_pairs_results')
    output_base.mkdir(exist_ok=True)
    
    # Check completed
    done = set()
    for d in output_base.iterdir():
        if d.is_dir() and (d / "chai1").exists():
            cifs = list((d / "chai1").glob("*.cif"))
            if cifs:
                done.add(d.name)
    
    remaining = [(r, l, rp, lp) for r, l, rp, lp in all_pairs if f"{r}_{l}" not in done]
    log(f"Already done: {len(done)}, Remaining: {len(remaining)}")
    
    if not remaining:
        log("All pairs completed!")
        return
    
    # Setup
    import torch
    n_gpus = torch.cuda.device_count()
    log(f"GPUs: {n_gpus} x {torch.cuda.get_device_name(0)}")
    
    from chai_lab.chai1 import run_inference
    
    results = []
    errors = []
    start_time = time.time()
    
    for idx, (r_name, l_name, r_path, l_path) in enumerate(remaining):
        pair_name = f"{r_name}_{l_name}"
        output_dir = output_base / pair_name
        gpu_id = idx % n_gpus
        
        # Progress
        pct = (idx + 1) / len(remaining) * 100
        elapsed = time.time() - start_time
        rate = (idx + 1) / elapsed if elapsed > 0 else 0
        eta = (len(remaining) - idx - 1) / rate / 3600 if rate > 0 else 0
        
        log(f"[{idx+1}/{len(remaining)} {pct:.1f}%] {pair_name} (GPU{gpu_id}) ETA: {eta:.1f}h")
        
        try:
            pair_start = time.time()
            output_dir.mkdir(parents=True, exist_ok=True)
            
            seq = extract_sequence(r_path)
            smiles = extract_smiles(l_path)
            
            if not seq:
                errors.append({'pair': pair_name, 'error': 'no_sequence'})
                log(f"  ⚠ Skipped: no sequence")
                continue
                
            if not smiles:
                errors.append({'pair': pair_name, 'error': 'no_smiles'})
                log(f"  ⚠ Skipped: no SMILES")
                continue
            
            if len(smiles) < 5:
                errors.append({'pair': pair_name, 'error': f'smiles_too_short:{len(smiles)}'})
                log(f"  ⚠ Skipped: SMILES too short ({len(smiles)})")
                continue
                
            if len(smiles) > 300:
                errors.append({'pair': pair_name, 'error': f'smiles_too_long:{len(smiles)}'})
                log(f"  ⚠ Skipped: SMILES too long ({len(smiles)})")
                continue
            
            # Create FASTA
            fasta_path = output_dir / "input.fasta"
            with open(fasta_path, 'w') as f:
                f.write(f">protein|name={r_name}\n{seq}\n>ligand|name={l_name}\n{smiles}\n")
            
            # Run Chai-1
            candidates = run_inference(
                fasta_file=fasta_path,
                output_dir=output_dir / "chai1",
                num_trunk_recycles=3,
                num_diffn_timesteps=200,
                seed=42,
                device=f"cuda:{gpu_id}",
                use_esm_embeddings=True,
            )
            
            pair_elapsed = time.time() - pair_start
            cifs = len(list((output_dir / "chai1").glob("*.cif")))
            results.append({'pair': pair_name, 'structures': cifs, 'time': pair_elapsed})
            log(f"  ✓ {cifs} structures in {pair_elapsed:.0f}s")
            
        except Exception as e:
            err = str(e)[:200]
            errors.append({'pair': pair_name, 'error': err})
            log(f"  ✗ Error: {err[:80]}")
        
        # Save progress every 10 pairs
        if (idx + 1) % 10 == 0:
            progress = {
                'timestamp': datetime.now().isoformat(),
                'completed': len(done) + len(results),
                'success': len(results),
                'errors': len(errors),
                'remaining': len(remaining) - idx - 1,
                'rate_per_hour': rate * 3600,
                'eta_hours': eta
            }
            with open(output_base / 'progress.json', 'w') as f:
                json.dump(progress, f, indent=2)
    
    # Final save
    total_time = time.time() - start_time
    final = {
        'timestamp': datetime.now().isoformat(),
        'total_time_hours': total_time / 3600,
        'success': results,
        'errors': errors,
        'summary': {
            'total_pairs': len(all_pairs),
            'successful': len(results),
            'errors': len(errors),
            'skipped_already_done': len(done)
        }
    }
    with open(output_base / 'final_results.json', 'w') as f:
        json.dump(final, f, indent=2)
    
    log("")
    log("="*70)
    log("BATCH COMPLETE")
    log("="*70)
    log(f"Total time: {total_time/3600:.1f} hours")
    log(f"Successful: {len(results)}")
    log(f"Errors: {len(errors)}")
    log(f"Results saved to: {output_base / 'final_results.json'}")

if __name__ == '__main__':
    main()
