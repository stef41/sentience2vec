#!/usr/bin/env python3
"""
Complete regeneration pipeline for all receptor-ligand pairs.

This script:
1. Downloads receptor structures from PDB
2. Prepares receptors (fix, clean, add hydrogens)
3. Creates ligand SDF files from SMILES
4. Runs Chai-1 predictions for all pairs
5. Generates GIF visualizations

Usage: python regenerate_all.py [--skip-download] [--skip-chai] [--gifs-only]
"""
import warnings
warnings.filterwarnings('ignore')

import os
os.environ['TOKENIZERS_PARALLELISM'] = 'false'

import sys
import json
import time
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp

def log(msg, level='INFO'):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{ts}] [{level}] {msg}", flush=True)

# ============================================================================
# STEP 1: Download and prepare receptors
# ============================================================================

def download_pdb(pdb_id, output_dir):
    """Download PDB structure from RCSB"""
    import urllib.request
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = output_dir / f"{pdb_id}.pdb"
    if output_path.exists():
        return output_path
    try:
        urllib.request.urlretrieve(url, output_path)
        return output_path
    except Exception as e:
        log(f"Failed to download {pdb_id}: {e}", 'ERROR')
        return None

def prepare_receptor(pdb_path, receptor_name, output_dir):
    """Prepare receptor structure for docking"""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        
        output_path = output_dir / f"{receptor_name}_prepared.pdb"
        if output_path.exists():
            return output_path
            
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        
        # Remove heterogens and water
        fixer.removeHeterogens(keepWater=False)
        
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        
        return output_path
    except Exception as e:
        log(f"Failed to prepare {receptor_name}: {e}", 'ERROR')
        return None

def prepare_all_receptors(receptor_db, output_base):
    """Download and prepare all receptors"""
    from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE
    
    pdb_dir = output_base / "pdb_downloads"
    receptors_dir = output_base / "receptors_full"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    receptors_dir.mkdir(parents=True, exist_ok=True)
    
    prepared = {}
    
    for receptor_name, info in RECEPTOR_PDB_DATABASE.items():
        pdb_id = info.get('pdb_id')
        if not pdb_id:
            log(f"Skipping {receptor_name} - no PDB structure", 'WARN')
            continue
        
        receptor_out_dir = receptors_dir / receptor_name
        receptor_out_dir.mkdir(parents=True, exist_ok=True)
        
        # Download
        pdb_path = download_pdb(pdb_id, pdb_dir)
        if not pdb_path:
            continue
        
        # Prepare
        prepared_path = prepare_receptor(pdb_path, receptor_name, receptor_out_dir)
        if prepared_path:
            prepared[receptor_name] = str(prepared_path)
            log(f"✓ Prepared {receptor_name} from {pdb_id}")
    
    return prepared

# ============================================================================
# STEP 2: Create ligand SDF files
# ============================================================================

def create_ligand_sdf(name, smiles, output_dir):
    """Create SDF file from SMILES"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        output_path = output_dir / f"{name}.sdf"
        if output_path.exists():
            return output_path
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        
        return output_path
    except Exception as e:
        log(f"Failed to create SDF for {name}: {e}", 'ERROR')
        return None

def prepare_all_ligands(output_base):
    """Create SDF files for all ligands"""
    from ligand_smiles_database import LIGAND_DATABASE
    
    ligands_dir = output_base / "ligands_full"
    ligands_dir.mkdir(parents=True, exist_ok=True)
    
    prepared = {}
    
    for name, info in LIGAND_DATABASE.items():
        smiles = info.get('smiles')
        if not smiles:
            continue
        
        ligand_dir = ligands_dir / name
        ligand_dir.mkdir(parents=True, exist_ok=True)
        
        sdf_path = create_ligand_sdf(name, smiles, ligand_dir)
        if sdf_path:
            prepared[name] = str(sdf_path)
            log(f"✓ Created ligand {name}")
    
    return prepared

# ============================================================================
# STEP 3: Run Chai-1 predictions
# ============================================================================

def extract_sequence(pdb_path):
    """Extract protein sequence from PDB"""
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
    """Extract SMILES from SDF"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_path))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

def run_chai1_worker(args):
    """Worker function for Chai-1 inference"""
    gpu_id, receptor_name, ligand_name, receptor_path, ligand_path, output_base = args
    
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    
    pair_name = f"{receptor_name}_{ligand_name}"
    output_dir = Path(output_base) / pair_name
    chai_dir = output_dir / "chai1"
    
    # Skip if exists
    if chai_dir.exists() and list(chai_dir.glob("*.cif")):
        return {'pair': pair_name, 'status': 'cached'}
    
    try:
        # Patch timeout
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
        
        # Create FASTA
        fasta_path = output_dir / "input.fasta"
        with open(fasta_path, 'w') as f:
            f.write(f">protein|name={receptor_name}\n{seq}\n")
            f.write(f">ligand|name={ligand_name}\n{smiles}\n")
        
        # Run inference
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

def run_all_chai1(receptors, ligands, output_base, n_gpus=8):
    """Run Chai-1 for all receptor-ligand pairs"""
    import torch
    
    n_gpus = min(n_gpus, torch.cuda.device_count())
    log(f"Running Chai-1 predictions on {n_gpus} GPUs")
    
    # Build pairs list
    pairs = []
    for r_name, r_path in receptors.items():
        for l_name, l_path in ligands.items():
            pairs.append((r_name, l_name, r_path, l_path))
    
    log(f"Total pairs to process: {len(pairs)}")
    
    # Save pairs list
    pairs_list = [[r, l, rp, lp] for r, l, rp, lp in pairs]
    with open(output_base / "pairs_list.json", 'w') as f:
        json.dump(pairs_list, f)
    
    # Process with multiprocessing
    results = []
    
    # Use spawn to avoid CUDA issues
    ctx = mp.get_context('spawn')
    
    with ctx.Pool(processes=n_gpus) as pool:
        tasks = []
        for idx, (r_name, l_name, r_path, l_path) in enumerate(pairs):
            gpu_id = idx % n_gpus
            tasks.append((gpu_id, r_name, l_name, r_path, l_path, str(output_base)))
        
        for i, result in enumerate(pool.imap_unordered(run_chai1_worker, tasks)):
            results.append(result)
            if (i + 1) % 100 == 0:
                success = sum(1 for r in results if r['status'] == 'success')
                cached = sum(1 for r in results if r['status'] == 'cached')
                errors = sum(1 for r in results if r['status'] == 'error')
                log(f"Progress: {i+1}/{len(pairs)} (success={success}, cached={cached}, errors={errors})")
    
    return results

# ============================================================================
# STEP 4: Generate GIFs
# ============================================================================

def generate_gif(pair_dir, output_path):
    """Generate GIF for a receptor-ligand pair"""
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import imageio
    
    # Find CIF file
    chai_dir = pair_dir / "chai1"
    cif_files = list(chai_dir.glob("*.cif"))
    if not cif_files:
        return False
    
    cif_path = cif_files[0]
    
    # Parse CIF
    coords = {'protein': [], 'ligand': []}
    elements = {'protein': [], 'ligand': []}
    
    try:
        with open(cif_path) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    parts = line.split()
                    if len(parts) >= 15:
                        try:
                            x = float(parts[10])
                            y = float(parts[11])
                            z = float(parts[12])
                            element = parts[2][0]
                            chain = parts[6] if len(parts) > 6 else 'A'
                            res_name = parts[5] if len(parts) > 5 else 'UNK'
                            
                            if chain == 'B' or res_name in ['LIG', 'LIG2', 'UNL']:
                                coords['ligand'].append([x, y, z])
                                elements['ligand'].append(element)
                            else:
                                coords['protein'].append([x, y, z])
                                elements['protein'].append(element)
                        except:
                            continue
    except:
        return False
    
    if not coords['ligand']:
        return False
    
    protein_coords = np.array(coords['protein']) if coords['protein'] else np.array([])
    ligand_coords = np.array(coords['ligand'])
    
    # Color maps
    ligand_colors = {'C': '#00FF00', 'N': '#00FFFF', 'O': '#FF69B4', 'S': '#FFD700', 'F': '#FF4500', 'H': '#FFFFFF'}
    protein_colors = {'C': '#404040', 'N': '#303050', 'O': '#503030', 'S': '#505030', 'H': '#303030'}
    
    # Generate frames
    frames = []
    n_frames = 36
    
    fig = plt.figure(figsize=(8, 8), facecolor='black')
    
    for i in range(n_frames):
        fig.clear()
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        
        angle = i * 10
        
        # Plot protein (translucent)
        if len(protein_coords) > 0:
            p_colors = [protein_colors.get(e, '#404040') for e in elements['protein']]
            ax.scatter(protein_coords[:, 0], protein_coords[:, 1], protein_coords[:, 2],
                      c=p_colors, s=20, alpha=0.15, depthshade=True)
        
        # Plot ligand (bright)
        l_colors = [ligand_colors.get(e, '#00FF00') for e in elements['ligand']]
        ax.scatter(ligand_coords[:, 0], ligand_coords[:, 1], ligand_coords[:, 2],
                  c=l_colors, s=150, alpha=1.0, edgecolors='white', linewidths=0.5)
        
        ax.view_init(elev=20, azim=angle)
        ax.set_axis_off()
        
        # Get center
        if len(protein_coords) > 0:
            all_coords = np.vstack([protein_coords, ligand_coords])
        else:
            all_coords = ligand_coords
        center = all_coords.mean(axis=0)
        max_range = np.max(np.abs(all_coords - center)) * 1.2
        
        ax.set_xlim(center[0] - max_range, center[0] + max_range)
        ax.set_ylim(center[1] - max_range, center[1] + max_range)
        ax.set_zlim(center[2] - max_range, center[2] + max_range)
        
        fig.canvas.draw()
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        frames.append(image)
    
    plt.close(fig)
    
    # Save GIF
    imageio.mimsave(str(output_path), frames, fps=10, loop=0)
    return True

def generate_all_gifs(output_base, n_workers=8):
    """Generate GIFs for all pairs"""
    gifs_dir = output_base / "gifs"
    gifs_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all pair directories with CIF files
    pairs = []
    for pair_dir in output_base.iterdir():
        if pair_dir.is_dir() and pair_dir.name != 'gifs':
            chai_dir = pair_dir / "chai1"
            if chai_dir.exists() and list(chai_dir.glob("*.cif")):
                gif_path = gifs_dir / f"{pair_dir.name}.gif"
                if not gif_path.exists():
                    pairs.append((pair_dir, gif_path))
    
    log(f"Generating {len(pairs)} GIFs")
    
    created = 0
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(generate_gif, p[0], p[1]): p for p in pairs}
        for i, future in enumerate(as_completed(futures)):
            pair_dir, gif_path = futures[future]
            try:
                if future.result():
                    created += 1
            except Exception as e:
                log(f"Error generating GIF for {pair_dir.name}: {e}", 'ERROR')
            
            if (i + 1) % 100 == 0:
                log(f"GIF progress: {i+1}/{len(pairs)} ({created} created)")
    
    return created

# ============================================================================
# MAIN
# ============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Regenerate all receptor-ligand data')
    parser.add_argument('--skip-download', action='store_true', help='Skip receptor download/prep')
    parser.add_argument('--skip-ligands', action='store_true', help='Skip ligand preparation')
    parser.add_argument('--skip-chai', action='store_true', help='Skip Chai-1 predictions')
    parser.add_argument('--gifs-only', action='store_true', help='Only generate GIFs')
    parser.add_argument('--n-gpus', type=int, default=8, help='Number of GPUs')
    args = parser.parse_args()
    
    output_base = Path('all_pairs_results')
    output_base.mkdir(parents=True, exist_ok=True)
    
    log("=" * 70)
    log("REGENERATING ALL RECEPTOR-LIGAND DATA")
    log("=" * 70)
    
    if args.gifs_only:
        log("GIFs only mode - skipping structure preparation")
        created = generate_all_gifs(output_base)
        log(f"Created {created} GIFs")
        return
    
    # Step 1: Prepare receptors
    if not args.skip_download:
        log("STEP 1: Preparing receptors...")
        from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE
        receptors = prepare_all_receptors(RECEPTOR_PDB_DATABASE, output_base.parent)
        log(f"Prepared {len(receptors)} receptors")
        
        with open(output_base / "prepared_receptors.json", 'w') as f:
            json.dump(receptors, f, indent=2)
    else:
        with open(output_base / "prepared_receptors.json") as f:
            receptors = json.load(f)
    
    # Step 2: Prepare ligands
    if not args.skip_ligands:
        log("STEP 2: Preparing ligands...")
        ligands = prepare_all_ligands(output_base.parent)
        log(f"Prepared {len(ligands)} ligands")
        
        with open(output_base / "prepared_ligands.json", 'w') as f:
            json.dump(ligands, f, indent=2)
    else:
        with open(output_base / "prepared_ligands.json") as f:
            ligands = json.load(f)
    
    # Step 3: Run Chai-1
    if not args.skip_chai:
        log("STEP 3: Running Chai-1 predictions...")
        results = run_all_chai1(receptors, ligands, output_base, n_gpus=args.n_gpus)
        
        success = sum(1 for r in results if r['status'] == 'success')
        cached = sum(1 for r in results if r['status'] == 'cached')
        errors = sum(1 for r in results if r['status'] == 'error')
        log(f"Chai-1 complete: {success} success, {cached} cached, {errors} errors")
        
        with open(output_base / "chai1_results.json", 'w') as f:
            json.dump(results, f, indent=2)
    
    # Step 4: Generate GIFs
    log("STEP 4: Generating GIFs...")
    created = generate_all_gifs(output_base)
    log(f"Created {created} GIFs")
    
    log("=" * 70)
    log("REGENERATION COMPLETE")
    log("=" * 70)

if __name__ == '__main__':
    mp.set_start_method('spawn', force=True)
    main()
