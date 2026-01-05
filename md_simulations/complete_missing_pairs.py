#!/usr/bin/env python3
"""
Complete the 228 missing pairs for buprenorphine, morphine, oxycodone
"""
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['TOKENIZERS_PARALLELISM'] = 'false'
os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2,3,4,5,6,7'

import json
import time
import sys
from pathlib import Path
from datetime import datetime
import multiprocessing as mp

def log(msg):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)

def extract_sequence(pdb_path):
    aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
          'GLY':'G','HIS':'H','HID':'H','HIE':'H','HIP':'H','ILE':'I','LEU':'L',
          'LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
          'TYR':'Y','VAL':'V','CYX':'C'}
    seq, seen = [], set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith('ATOM'):
                    rn = line[17:20].strip()
                    ri = line[22:26].strip()
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

def run_chai1_worker(args):
    gpu_id, receptor_name, ligand_name, receptor_path, ligand_path, output_base = args
    
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    
    pair_name = f"{receptor_name}_{ligand_name}"
    output_dir = Path(output_base) / pair_name
    chai_dir = output_dir / "chai1"
    
    if chai_dir.exists() and list(chai_dir.glob("*.cif")):
        return {'pair': pair_name, 'status': 'cached'}
    
    try:
        try:
            import chai_lab.utils.timeout as timeout_module
            def no_timeout(seconds):
                def decorator(func):
                    return func
                return decorator
            timeout_module.timeout = no_timeout
        except:
            pass
        
        from chai_lab.chai1 import run_inference
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        seq = extract_sequence(receptor_path)
        smiles = extract_smiles(ligand_path)
        
        if not seq or not smiles:
            return {'pair': pair_name, 'status': 'error', 'error': 'extraction'}
        
        fasta_path = output_dir / "input.fasta"
        with open(fasta_path, 'w') as f:
            f.write(f">protein|name={receptor_name}\n{seq}\n")
            f.write(f">ligand|name={ligand_name}\n{smiles}\n")
        
        start = time.time()
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=chai_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device="cuda:0",
            use_esm_embeddings=True,
        )
        elapsed = time.time() - start
        
        cif_count = len(list(chai_dir.glob("*.cif")))
        return {'pair': pair_name, 'status': 'success', 'structures': cif_count, 'time': elapsed}
        
    except Exception as e:
        return {'pair': pair_name, 'status': 'error', 'error': str(e)[:200]}

def main():
    log("="*70)
    log("COMPLETING MISSING PAIRS: buprenorphine, morphine, oxycodone")
    log("="*70)
    
    with open('all_pairs_results/prepared_receptors.json') as f:
        receptors = json.load(f)
    
    missing_ligands = {
        'buprenorphine': 'ligands_full/buprenorphine/buprenorphine.sdf',
        'morphine': 'ligands_full/morphine/morphine.sdf',
        'oxycodone': 'ligands_full/oxycodone/oxycodone.sdf'
    }
    
    pairs = []
    for r_name, r_path in receptors.items():
        for l_name, l_path in missing_ligands.items():
            pairs.append((r_name, l_name, r_path, l_path))
    
    log(f"Processing {len(pairs)} pairs on 8 GPUs")
    
    output_base = Path('all_pairs_results')
    
    results = []
    ctx = mp.get_context('spawn')
    
    with ctx.Pool(processes=8) as pool:
        tasks = []
        for idx, (r_name, l_name, r_path, l_path) in enumerate(pairs):
            gpu_id = idx % 8
            tasks.append((gpu_id, r_name, l_name, r_path, l_path, str(output_base)))
        
        for i, result in enumerate(pool.imap_unordered(run_chai1_worker, tasks)):
            results.append(result)
            if (i + 1) % 10 == 0:
                success = sum(1 for r in results if r['status'] == 'success')
                cached = sum(1 for r in results if r['status'] == 'cached')
                errors = sum(1 for r in results if r['status'] == 'error')
                log(f"Progress: {i+1}/{len(pairs)} (success={success}, cached={cached}, errors={errors})")
    
    success = sum(1 for r in results if r['status'] == 'success')
    cached = sum(1 for r in results if r['status'] == 'cached')
    errors = sum(1 for r in results if r['status'] == 'error')
    
    log("="*70)
    log(f"COMPLETE: {success} success, {cached} cached, {errors} errors")
    log("="*70)
    
    with open('all_pairs_results/missing_pairs_results.json', 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == '__main__':
    mp.set_start_method('spawn', force=True)
    main()
