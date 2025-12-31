#!/usr/bin/env python3
"""Retry failed Chai-1 predictions"""
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import warnings
warnings.filterwarnings('ignore')

import json
import time
import torch
import torch.multiprocessing as mp
from pathlib import Path

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def sanitize_name(name):
    """Replace problematic Unicode chars"""
    replacements = {'Α': 'Alpha-', 'Β': 'B', 'Δ': 'Delta-', 'α': 'alpha-', 'β': 'beta-', 'δ': 'delta-'}
    for old, new in replacements.items():
        name = name.replace(old, new)
    return name

def extract_sequence(pdb_path):
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'MSE': 'M', 'HID': 'H', 'HIE': 'H', 'HIP': 'H'
    }
    seq = []
    seen = set()
    with open(pdb_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ATOM'):
                res = line[17:20].strip()
                chain = line[21]
                resnum = line[22:26].strip()
                key = (chain, resnum)
                if key not in seen and res in aa_map:
                    seen.add(key)
                    seq.append(aa_map[res])
    return ''.join(seq)

def extract_smiles(sdf_path):
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol:
        return Chem.MolToSmiles(mol)
    return None

def worker(gpu_id, task_queue, result_queue):
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    import torch
    torch.cuda.set_device(0)
    
    from chai_lab.chai1 import run_inference
    
    output_base = Path('all_pairs_results')
    
    while True:
        task = task_queue.get()
        if task is None:
            break
        
        idx, total, r_name, l_name, r_path, l_path = task
        pair_name = sanitize_name(f"{r_name}_{l_name}")
        pair_dir = output_base / pair_name
        chai_dir = pair_dir / "chai1"
        
        try:
            # Skip if already done
            if chai_dir.exists() and list(chai_dir.glob("*.cif")):
                result_queue.put(('skip', pair_name))
                continue
            
            chai_dir.mkdir(parents=True, exist_ok=True)
            
            seq = extract_sequence(r_path)
            smiles = extract_smiles(l_path)
            
            if not seq or not smiles:
                result_queue.put(('error', pair_name, "No seq/smiles"))
                continue
            
            fasta = f">protein|name=receptor\n{seq}\n>ligand|name=ligand\n{smiles}"
            fasta_file = pair_dir / "input.fasta"
            with open(fasta_file, 'w', encoding='utf-8') as f:
                f.write(fasta)
            
            start = time.time()
            run_inference(
                fasta_file=fasta_file,
                output_dir=chai_dir,
                num_trunk_recycles=3,
                num_diffn_timesteps=200,
                seed=42,
                device=torch.device('cuda:0'),
                use_esm_embeddings=True
            )
            elapsed = time.time() - start
            
            cifs = list(chai_dir.glob("*.cif"))
            result_queue.put(('success', pair_name, len(cifs), elapsed))
            
        except Exception as e:
            result_queue.put(('error', pair_name, str(e)[:100]))

def main():
    mp.set_start_method('spawn', force=True)
    
    output_base = Path('all_pairs_results')
    
    with open(output_base / 'retry_pairs.json') as f:
        retry_pairs = json.load(f)
    
    log(f"Retrying {len(retry_pairs)} failed pairs")
    
    num_gpus = 8
    task_queue = mp.Queue()
    result_queue = mp.Queue()
    
    # Start workers
    workers = []
    for gpu_id in range(num_gpus):
        p = mp.Process(target=worker, args=(gpu_id, task_queue, result_queue))
        p.start()
        workers.append(p)
        log(f"Started worker on GPU {gpu_id}")
    
    # Queue tasks
    for i, (r, l, rp, lp) in enumerate(retry_pairs):
        task_queue.put((i, len(retry_pairs), r, l, rp, lp))
    
    # Send stop signals
    for _ in range(num_gpus):
        task_queue.put(None)
    
    # Collect results
    success = 0
    errors = 0
    skipped = 0
    
    for _ in range(len(retry_pairs)):
        result = result_queue.get()
        if result[0] == 'success':
            success += 1
            log(f"[GPU] ✓ {result[1]}: {result[2]} structures ({result[3]:.0f}s)")
        elif result[0] == 'skip':
            skipped += 1
        else:
            errors += 1
            log(f"[GPU] ✗ {result[1]}: {result[2]}")
    
    for p in workers:
        p.join()
    
    log(f"Done! Success: {success}, Errors: {errors}, Skipped: {skipped}")

if __name__ == '__main__':
    main()
