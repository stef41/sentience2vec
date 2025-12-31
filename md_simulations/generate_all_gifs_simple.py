#!/usr/bin/env python3
"""
Generate rotating 3D molecular structure GIFs for ALL receptor-ligand pairs.
Simple sequential version with clear progress.
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
import sys

ATOM_COLORS = {
    'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'S': '#FFFF30',
    'H': '#FFFFFF', 'P': '#FF8000', 'CL': '#1FF01F', 'BR': '#A62929',
    'F': '#90E050', 'I': '#940094',
}

def parse_pdb(pdb_file):
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
                    atoms.append((x, y, z, element.upper()))
                except:
                    continue
    return atoms

def parse_cif(cif_file):
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
                if 'Cartn_x' in col_name: x_col = col_idx
                elif 'Cartn_y' in col_name: y_col = col_idx
                elif 'Cartn_z' in col_name: z_col = col_idx
                elif 'type_symbol' in col_name: elem_col = col_idx
                col_idx += 1
                i += 1
            continue
        if x_col is not None and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop') and line:
            parts = line.split()
            if len(parts) > max(x_col, y_col, z_col):
                try:
                    x, y, z = float(parts[x_col]), float(parts[y_col]), float(parts[z_col])
                    elem = parts[elem_col].upper() if elem_col and len(parts) > elem_col else 'C'
                    atoms.append((x, y, z, elem))
                except: pass
        i += 1
    return atoms

def create_gif(atoms, output_path, n_frames=24, dpi=80):
    if not atoms:
        return False
    coords = np.array([(a[0], a[1], a[2]) for a in atoms])
    elements = [a[3] for a in atoms]
    coords = coords - coords.mean(axis=0)
    
    if len(atoms) > 2000:
        step = len(atoms) // 2000
        coords, elements = coords[::step], elements[::step]
    
    colors = [ATOM_COLORS.get(e, '#808080') for e in elements]
    frames = []
    
    for i in range(n_frames):
        fig = plt.figure(figsize=(4, 4), facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        ax.view_init(elev=20, azim=i * (360 / n_frames))
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=8, alpha=0.8, edgecolors='none')
        ax.set_axis_off()
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=dpi, facecolor='black', bbox_inches='tight', pad_inches=0.1)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)
    
    frames[0].save(output_path, save_all=True, append_images=frames[1:], duration=100, loop=0, optimize=True)
    return True

def main():
    output_base = Path('all_pairs_results')
    gif_dir = output_base / 'gifs'
    gif_dir.mkdir(exist_ok=True)
    
    all_pairs = sorted([d for d in output_base.iterdir() if d.is_dir() and d.name != 'gifs'])
    total = len(all_pairs)
    
    print(f"Processing {total} pairs...", flush=True)
    
    created = skipped = no_struct = 0
    
    for i, pair_dir in enumerate(all_pairs):
        gif_path = gif_dir / f"{pair_dir.name}.gif"
        
        if gif_path.exists():
            skipped += 1
            continue
        
        pdb = pair_dir / 'md' / 'minimized.pdb'
        cif = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
        
        atoms = None
        if pdb.exists(): atoms = parse_pdb(pdb)
        elif cif.exists(): atoms = parse_cif(cif)
        
        if atoms and len(atoms) > 10:
            try:
                if create_gif(atoms, gif_path):
                    created += 1
            except Exception as e:
                pass
        else:
            no_struct += 1
        
        if (i + 1) % 50 == 0:
            print(f"[{i+1}/{total}] Created: {created}, Skipped: {skipped}, No struct: {no_struct}", flush=True)
    
    print(f"\nDONE: Created {created}, Skipped {skipped}, No struct {no_struct}")
    print(f"Total GIFs: {len(list(gif_dir.glob('*.gif')))}")

if __name__ == '__main__':
    main()
