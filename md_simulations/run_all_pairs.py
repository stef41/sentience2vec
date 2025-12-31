#!/usr/bin/env python3
"""
Batch Chai-1 inference for all receptor-ligand pairs
Uses all 8 H100 GPUs in parallel
"""
import warnings
warnings.filterwarnings('ignore')
import os
import sys
import json
import time
import torch
import torch.multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
from datetime import datetime

# Set multiprocessing
mp.set_start_method('spawn', force=True)

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def extract_sequence_from_pdb(pdb_path):
    """Extract protein sequence from PDB"""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C'
    }
    seq = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                chain = line[21]
                key = (chain, res_num)
                if key not in seen and res_name in aa_map:
                    seq.append(aa_map[res_name])
                    seen.add(key)
    return ''.join(seq)

def extract_smiles_from_sdf(sdf_path):
    """Extract SMILES from SDF file"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_path))
        if mol:
            return Chem.MolToSmiles(mol)
    except:
        pass
    
    # Fallback: read from file
    with open(sdf_path) as f:
        content = f.read()
    if 'SMILES' in content:
        for line in content.split('\n'):
            if line.strip() and not line.startswith('>') and len(line) < 500:
                if any(c in line for c in ['C', 'N', 'O', 'c', 'n']):
                    return line.strip()
    return None

def process_pair(args):
    """Process single receptor-ligand pair"""
    idx, total, r_name, l_name, r_path, l_path, gpu_id, output_base = args
    
    pair_name = f"{r_name}_{l_name}"
    output_dir = Path(output_base) / pair_name
    
    # Skip if done
    chai_dir = output_dir / "chai1"
    if chai_dir.exists() and list(chai_dir.glob("*.cif")):
        return {'pair': pair_name, 'status': 'cached', 'time': 0}
    
    try:
        start = time.time()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract sequence and SMILES
        seq = extract_sequence_from_pdb(r_path)
        smiles = extract_smiles_from_sdf(l_path)
        
        if not seq or not smiles:
            return {'pair': pair_name, 'status': 'failed', 'error': 'extraction'}
        
        # Prepare FASTA
        fasta_path = output_dir / "input.fasta"
        with open(fasta_path, 'w') as f:
            f.write(f">protein|name={r_name}\n{seq}\n>ligand|name={l_name}\n{smiles}\n")
        
        # Run Chai-1
        from chai_lab.chai1 import run_inference
        
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=chai_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device=f"cuda:{gpu_id}",
            use_esm_embeddings=True,
        )
        
        elapsed = time.time() - start
        cif_count = len(list(chai_dir.glob("*.cif")))
        
        return {'pair': pair_name, 'status': 'success', 'structures': cif_count, 'time': elapsed, 'gpu': gpu_id}
        
    except Exception as e:
        return {'pair': pair_name, 'status': 'failed', 'error': str(e)[:200]}

def worker(gpu_id, pairs_queue, results_queue, output_base):
    """Worker process for a single GPU"""
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    
    while True:
        item = pairs_queue.get()
        if item is None:
            break
        
        idx, total, r_name, l_name, r_path, l_path = item
        result = process_pair((idx, total, r_name, l_name, r_path, l_path, 0, output_base))
        result['gpu'] = gpu_id
        results_queue.put(result)

def main():
    log("="*70)
    log("BATCH CHAI-1: ALL RECEPTOR-LIGAND PAIRS")
    log("="*70)
    
    # Load pairs
    pairs_file = Path('all_pairs_results/pairs_list.json')
    with open(pairs_file) as f:
        all_pairs = json.load(f)
    
    log(f"Loaded {len(all_pairs):,} pairs")
    
    # Check GPUs
    n_gpus = torch.cuda.device_count()
    log(f"GPUs available: {n_gpus}")
    
    output_base = Path('all_pairs_results')
    
    # Check already done
    done = set()
    for d in output_base.iterdir():
        if d.is_dir() and (d / "chai1").exists():
            cifs = list((d / "chai1").glob("*.cif"))
            if cifs:
                done.add(d.name)
    
    remaining = [(r, l, rp, lp) for r, l, rp, lp in all_pairs if f"{r}_{l}" not in done]
    log(f"Already done: {len(done)}, Remaining: {len(remaining)}")
    
    if not remaining:
        log("All pairs already processed!")
        return
    
    # Process sequentially but rotate GPUs
    results = []
    
    log(f"\nProcessing {len(remaining)} pairs across {n_gpus} GPUs...")
    
    from chai_lab.chai1 import run_inference
    
    for idx, (r_name, l_name, r_path, l_path) in enumerate(tqdm(remaining, desc="Pairs", ncols=100)):
        gpu_id = idx % n_gpus
        
        result = process_pair((
            idx, len(remaining), r_name, l_name, r_path, l_path, 
            gpu_id, str(output_base)
        ))
        results.append(result)
        
        # Save progress every 10 pairs
        if (idx + 1) % 10 == 0:
            with open(output_base / 'progress.json', 'w') as f:
                json.dump({
                    'completed': len(results),
                    'total': len(remaining),
                    'results': results[-100:]  # Last 100
                }, f, indent=2)
        
        if result['status'] == 'success':
            tqdm.write(f"  ✓ {result['pair']}: {result.get('structures', 0)} structures ({result.get('time', 0):.1f}s) [GPU{gpu_id}]")
        elif result['status'] == 'cached':
            pass  # Silent for cached
        else:
            tqdm.write(f"  ✗ {result['pair']}: {result.get('error', 'unknown')[:50]}")
    
    # Final save
    with open(output_base / 'all_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    # Summary
    success = sum(1 for r in results if r['status'] == 'success')
    cached = sum(1 for r in results if r['status'] == 'cached')
    failed = sum(1 for r in results if r['status'] == 'failed')
    
    log(f"\n{'='*70}")
    log("BATCH COMPLETE")
    log(f"{'='*70}")
    log(f"Success: {success} | Cached: {cached} | Failed: {failed}")
    log(f"Total time: {sum(r.get('time', 0) for r in results)/3600:.1f} hours")
    log(f"Results: {output_base}")

if __name__ == '__main__':
    main()
