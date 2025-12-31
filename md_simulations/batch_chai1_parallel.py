#!/usr/bin/env python3
"""
Parallel Chai-1 batch processing using all 8 GPUs simultaneously.
Each GPU runs as an independent worker processing different pairs.
"""
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['TOKENIZERS_PARALLELISM'] = 'false'

import json
import time
import sys
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from queue import Empty

def log(msg, worker_id=None):
    ts = datetime.now().strftime('%H:%M:%S')
    prefix = f"[{ts}]" if worker_id is None else f"[{ts}][GPU{worker_id}]"
    print(f"{prefix} {msg}", flush=True)

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

def worker_process(gpu_id, task_queue, result_queue, output_base):
    """Worker that processes pairs on a specific GPU."""
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    os.environ['PYTHONIOENCODING'] = 'utf-8'
    
    # Patch timeout before importing chai_lab
    try:
        import chai_lab.utils.timeout as timeout_module
        def no_timeout(seconds):
            def decorator(func):
                return func
            return decorator
        timeout_module.timeout = no_timeout
    except:
        pass
    
    import torch
    from chai_lab.chai1 import run_inference
    
    log(f"Worker started on {torch.cuda.get_device_name(0)}", gpu_id)
    
    while True:
        try:
            task = task_queue.get(timeout=5)
        except Empty:
            # Check if we should exit
            if task_queue.empty():
                log("No more tasks, exiting", gpu_id)
                break
            continue
        
        if task is None:  # Poison pill
            log("Received shutdown signal", gpu_id)
            break
        
        idx, total, r_name, l_name, r_path, l_path = task
        pair_name = f"{r_name}_{l_name}"
        output_dir = Path(output_base) / pair_name
        
        try:
            start = time.time()
            output_dir.mkdir(parents=True, exist_ok=True)
            
            seq = extract_sequence(r_path)
            smiles = extract_smiles(l_path)
            
            if not seq:
                result_queue.put({'pair': pair_name, 'status': 'error', 'error': 'no_sequence'})
                continue
            if not smiles:
                result_queue.put({'pair': pair_name, 'status': 'error', 'error': 'no_smiles'})
                continue
            if len(smiles) < 5 or len(smiles) > 300:
                result_queue.put({'pair': pair_name, 'status': 'error', 'error': f'smiles_len:{len(smiles)}'})
                continue
            
            # Sanitize pair name for filesystem
            pair_name_safe = pair_name.encode('ascii', 'replace').decode('ascii').replace('?', '_')
            fasta_path = output_dir / "input.fasta"
            with open(fasta_path, 'w', encoding='utf-8') as f:
                f.write(f">protein|name={r_name}\n{seq}\n>ligand|name={l_name}\n{smiles}\n")
            
            # Run on device 0 since CUDA_VISIBLE_DEVICES is set
            candidates = run_inference(
                fasta_file=fasta_path,
                output_dir=output_dir / "chai1",
                num_trunk_recycles=3,
                num_diffn_timesteps=200,
                seed=42,
                device="cuda:0",
                use_esm_embeddings=True,
            )
            
            elapsed = time.time() - start
            cifs = len(list((output_dir / "chai1").glob("*.cif")))
            log(f"✓ {pair_name}: {cifs} structures ({elapsed:.0f}s)", gpu_id)
            result_queue.put({'pair': pair_name, 'status': 'success', 'structures': cifs, 'time': elapsed})
            
        except Exception as e:
            log(f"✗ {pair_name}: {str(e)[:60]}", gpu_id)
            result_queue.put({'pair': pair_name, 'status': 'error', 'error': str(e)[:200]})

def main():
    mp.set_start_method('spawn', force=True)
    
    log("="*70)
    log("PARALLEL CHAI-1: 8 GPUs")
    log("="*70)
    
    # Load pairs
    pairs_file = Path('all_pairs_results/pairs_list.json')
    with open(pairs_file) as f:
        all_pairs = json.load(f)
    log(f"Total pairs: {len(all_pairs):,}")
    
    output_base = Path('all_pairs_results')
    
    # Check completed
    done = set()
    for d in output_base.iterdir():
        if d.is_dir() and (d / "chai1").exists():
            cifs = list((d / "chai1").glob("*.cif"))
            if cifs:
                done.add(d.name)
    
    remaining = [(i, len(all_pairs), r, l, rp, lp) 
                 for i, (r, l, rp, lp) in enumerate(all_pairs) 
                 if f"{r}_{l}" not in done]
    
    log(f"Already done: {len(done)}, Remaining: {len(remaining)}")
    
    if not remaining:
        log("All pairs completed!")
        return
    
    # Setup queues
    n_gpus = 8
    task_queue = mp.Queue()
    result_queue = mp.Queue()
    
    # Add all tasks
    for task in remaining:
        task_queue.put(task)
    
    # Add poison pills
    for _ in range(n_gpus):
        task_queue.put(None)
    
    log(f"Starting {n_gpus} workers...")
    
    # Start workers
    workers = []
    for gpu_id in range(n_gpus):
        p = mp.Process(target=worker_process, args=(gpu_id, task_queue, result_queue, str(output_base)))
        p.start()
        workers.append(p)
        time.sleep(2)  # Stagger starts
    
    log(f"All {n_gpus} workers started!")
    
    # Monitor progress
    results = []
    start_time = time.time()
    last_save = time.time()
    
    while any(w.is_alive() for w in workers):
        try:
            result = result_queue.get(timeout=1)
            results.append(result)
            
            # Progress update
            completed = len(done) + len(results)
            elapsed = time.time() - start_time
            rate = len(results) / elapsed if elapsed > 0 else 0
            eta = (len(remaining) - len(results)) / rate / 3600 if rate > 0 else 0
            
            success = sum(1 for r in results if r['status'] == 'success')
            errors = sum(1 for r in results if r['status'] == 'error')
            
            if len(results) % 10 == 0:
                log(f"Progress: {completed}/{len(all_pairs)} ({success} ok, {errors} err) | Rate: {rate*3600:.0f}/hr | ETA: {eta:.1f}h")
            
            # Save progress periodically
            if time.time() - last_save > 60:
                progress = {
                    'timestamp': datetime.now().isoformat(),
                    'completed': completed,
                    'success': success,
                    'errors': errors,
                    'remaining': len(remaining) - len(results),
                    'rate_per_hour': rate * 3600,
                    'eta_hours': eta
                }
                with open(output_base / 'progress.json', 'w') as f:
                    json.dump(progress, f, indent=2)
                last_save = time.time()
                
        except Empty:
            continue
    
    # Wait for workers
    for w in workers:
        w.join()
    
    # Drain remaining results
    while not result_queue.empty():
        results.append(result_queue.get_nowait())
    
    # Final save
    total_time = time.time() - start_time
    success_results = [r for r in results if r['status'] == 'success']
    error_results = [r for r in results if r['status'] == 'error']
    
    final = {
        'timestamp': datetime.now().isoformat(),
        'total_time_hours': total_time / 3600,
        'success': success_results,
        'errors': error_results,
        'summary': {
            'total_pairs': len(all_pairs),
            'successful': len(success_results),
            'errors': len(error_results),
            'already_done': len(done)
        }
    }
    with open(output_base / 'final_results.json', 'w') as f:
        json.dump(final, f, indent=2)
    
    log("")
    log("="*70)
    log("COMPLETE")
    log("="*70)
    log(f"Time: {total_time/3600:.1f} hours")
    log(f"Success: {len(success_results)}")
    log(f"Errors: {len(error_results)}")
    log(f"Rate: {len(results)/(total_time/3600):.0f} pairs/hour")

if __name__ == '__main__':
    main()
