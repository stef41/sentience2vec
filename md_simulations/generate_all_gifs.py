#!/usr/bin/env python3
"""
Generate rotating 3D molecular structure GIFs for ALL receptor-ligand pairs.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import imageio
from PIL import Image
import io
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Color scheme for atoms
ATOM_COLORS = {
    'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'S': '#FFFF30',
    'H': '#FFFFFF', 'P': '#FF8000', 'CL': '#1FF01F', 'BR': '#A62929',
    'F': '#90E050', 'I': '#940094',
}

def parse_pdb(pdb_file):
    """Parse PDB file and extract atom coordinates."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip()
                    if not element:
                        element = line[12:14].strip()[0]
                    element = element.upper()
                    atoms.append((x, y, z, element))
                except:
                    continue
    return atoms

def parse_cif(cif_file):
    """Parse CIF file and extract atom coordinates."""
    atoms = []
    x_col, y_col, z_col, elem_col = None, None, None, None
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('_atom_site.'):
            col_idx = 0
            while i < len(lines) and lines[i].strip().startswith('_atom_site.'):
                col_name = lines[i].strip()
                if 'Cartn_x' in col_name:
                    x_col = col_idx
                elif 'Cartn_y' in col_name:
                    y_col = col_idx
                elif 'Cartn_z' in col_name:
                    z_col = col_idx
                elif 'type_symbol' in col_name:
                    elem_col = col_idx
                col_idx += 1
                i += 1
            continue
        
        if x_col is not None and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop') and line:
            parts = line.split()
            if len(parts) > max(x_col, y_col, z_col):
                try:
                    x = float(parts[x_col])
                    y = float(parts[y_col])
                    z = float(parts[z_col])
                    elem = parts[elem_col].upper() if elem_col and len(parts) > elem_col else 'C'
                    atoms.append((x, y, z, elem))
                except:
                    pass
        i += 1
    
    return atoms

def create_rotation_gif(atoms, output_path, n_frames=24, dpi=80):
    """Create a rotating 3D view GIF."""
    if not atoms:
        return False
    
    coords = np.array([(a[0], a[1], a[2]) for a in atoms])
    elements = [a[3] for a in atoms]
    
    # Center coordinates
    center = coords.mean(axis=0)
    coords = coords - center
    
    # Subsample if too many atoms
    if len(atoms) > 2000:
        step = len(atoms) // 2000
        coords = coords[::step]
        elements = elements[::step]
    
    colors = [ATOM_COLORS.get(e, '#808080') for e in elements]
    
    frames = []
    for i in range(n_frames):
        fig = plt.figure(figsize=(4, 4), facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        
        angle = i * (360 / n_frames)
        ax.view_init(elev=20, azim=angle)
        
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                   c=colors, s=8, alpha=0.8, edgecolors='none')
        
        ax.set_xlim([coords[:, 0].min(), coords[:, 0].max()])
        ax.set_ylim([coords[:, 1].min(), coords[:, 1].max()])
        ax.set_zlim([coords[:, 2].min(), coords[:, 2].max()])
        ax.set_axis_off()
        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=dpi, facecolor='black',
                    bbox_inches='tight', pad_inches=0.1)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)
    
    frames[0].save(output_path, save_all=True, append_images=frames[1:],
                   duration=100, loop=0, optimize=True)
    return True

def process_pair(pair_info):
    """Process a single pair - worker function for parallel processing."""
    pair_dir, gif_path = pair_info
    pair_name = pair_dir.name
    
    if gif_path.exists():
        return (pair_name, 'skipped')
    
    # Try minimized PDB first, then CIF
    pdb_file = pair_dir / 'md' / 'minimized.pdb'
    cif_file = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
    
    atoms = None
    if pdb_file.exists():
        atoms = parse_pdb(pdb_file)
    elif cif_file.exists():
        atoms = parse_cif(cif_file)
    
    if atoms and len(atoms) > 10:
        try:
            if create_rotation_gif(atoms, gif_path):
                return (pair_name, 'created')
        except Exception as e:
            return (pair_name, f'error: {e}')
    
    return (pair_name, 'no_structure')

def main():
    output_base = Path('all_pairs_results')
    gif_dir = output_base / 'gifs'
    gif_dir.mkdir(exist_ok=True)
    
    # Get all pair directories
    all_pairs = [d for d in output_base.iterdir() 
                 if d.is_dir() and d.name != 'gifs']
    
    print(f"Found {len(all_pairs)} receptor-ligand pairs")
    print("="*60)
    
    # Prepare work items
    work_items = []
    for pair_dir in sorted(all_pairs):
        gif_path = gif_dir / f"{pair_dir.name}.gif"
        work_items.append((pair_dir, gif_path))
    
    # Process in parallel
    created = 0
    skipped = 0
    no_structure = 0
    errors = 0
    
    n_workers = min(8, multiprocessing.cpu_count())
    print(f"Processing with {n_workers} workers...")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_pair, item): item for item in work_items}
        
        for i, future in enumerate(as_completed(futures)):
            pair_name, status = future.result()
            
            if status == 'created':
                created += 1
            elif status == 'skipped':
                skipped += 1
            elif status == 'no_structure':
                no_structure += 1
            else:
                errors += 1
            
            if (i + 1) % 100 == 0:
                print(f"Progress: {i+1}/{len(work_items)} | Created: {created} | Skipped: {skipped}")
    
    print(f"\n{'='*60}")
    print(f"COMPLETE:")
    print(f"  Created: {created}")
    print(f"  Skipped (existing): {skipped}")
    print(f"  No structure: {no_structure}")
    print(f"  Errors: {errors}")
    print(f"  Total GIFs: {len(list(gif_dir.glob('*.gif')))}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
