#!/usr/bin/env python3
"""Test GIF generation with ligand highlighting on a few key pairs."""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
from PIL import Image
import io

PROTEIN_COLORS = {
    'C': '#8888AA', 'N': '#4060C0', 'O': '#C04040', 'S': '#C0C040',
    'H': '#AAAAAA', 'P': '#C08040',
}

LIGAND_COLORS = {
    'C': '#00FF00', 'N': '#00CCFF', 'O': '#FF3333', 'S': '#FFFF00',
    'H': '#AAFFAA', 'P': '#FF8800',
}

def parse_cif_with_ligand(cif_file):
    protein_atoms, ligand_atoms = [], []
    x_col, y_col, z_col, elem_col, comp_col = None, None, None, None, None
    
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
                elif 'label_comp_id' in col_name: comp_col = col_idx
                col_idx += 1
                i += 1
            continue
        
        if x_col is not None and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop') and line:
            parts = line.split()
            if len(parts) > max(filter(None, [x_col, y_col, z_col])):
                try:
                    x, y, z = float(parts[x_col]), float(parts[y_col]), float(parts[z_col])
                    elem = parts[elem_col].upper() if elem_col and len(parts) > elem_col else 'C'
                    comp = parts[comp_col] if comp_col and len(parts) > comp_col else ''
                    
                    if 'LIG' in comp.upper():
                        ligand_atoms.append((x, y, z, elem))
                    else:
                        protein_atoms.append((x, y, z, elem))
                except:
                    pass
        i += 1
    
    return protein_atoms, ligand_atoms

def create_gif(protein_atoms, ligand_atoms, output_path, title=""):
    all_atoms = protein_atoms + ligand_atoms
    if not all_atoms:
        return False
    
    all_coords = np.array([(a[0], a[1], a[2]) for a in all_atoms])
    center = all_coords.mean(axis=0)
    
    # Protein: subsample, muted colors
    if protein_atoms:
        p_coords = np.array([(a[0], a[1], a[2]) for a in protein_atoms]) - center
        p_elements = [a[3] for a in protein_atoms]
        if len(protein_atoms) > 1000:
            step = len(protein_atoms) // 1000
            p_coords = p_coords[::step]
            p_elements = p_elements[::step]
        p_colors = [PROTEIN_COLORS.get(e, '#404040') for e in p_elements]
    else:
        p_coords, p_colors = np.array([]).reshape(0,3), []
    
    # Ligand: keep all, bright colors
    if ligand_atoms:
        l_coords = np.array([(a[0], a[1], a[2]) for a in ligand_atoms]) - center
        l_elements = [a[3] for a in ligand_atoms]
        l_colors = [LIGAND_COLORS.get(e, '#00FF00') for e in l_elements]
    else:
        l_coords, l_colors = np.array([]).reshape(0,3), []
    
    frames = []
    n_frames = 36
    
    for i in range(n_frames):
        fig = plt.figure(figsize=(6, 6), facecolor='#0a0a0a')
        ax = fig.add_subplot(111, projection='3d', facecolor='#0a0a0a')
        ax.view_init(elev=15, azim=i * (360 / n_frames))
        
        # Protein (visible, larger)
        if len(p_coords) > 0:
            ax.scatter(p_coords[:, 0], p_coords[:, 1], p_coords[:, 2],
                       c=p_colors, s=15, alpha=0.7, edgecolors='none')
        
        # Ligand (large, bright, with edge)
        if len(l_coords) > 0:
            ax.scatter(l_coords[:, 0], l_coords[:, 1], l_coords[:, 2],
                       c=l_colors, s=120, alpha=1.0, edgecolors='white', linewidths=0.8)
        
        ax.set_axis_off()
        if title:
            ax.set_title(title, color='white', fontsize=10, pad=-5)
        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, facecolor='#0a0a0a', bbox_inches='tight', pad_inches=0.1)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)
    
    frames[0].save(output_path, save_all=True, append_images=frames[1:], duration=80, loop=0, optimize=True)
    return True

# Test on key pairs
base = Path('all_pairs_results')
gif_dir = base / 'gifs'

test_pairs = ['5-HT2A_DMT', '5-HT2A_LSD', '5-HT2A_Psilocin', 'D2_LSD', 'CB1_THC-O-acetate', 'MOR_Morphine']

for pair in test_pairs:
    cif = base / pair / 'chai1' / 'pred.model_idx_0.cif'
    if cif.exists():
        protein, ligand = parse_cif_with_ligand(cif)
        print(f"{pair}: {len(protein)} protein, {len(ligand)} ligand atoms")
        
        out = gif_dir / f"{pair}_highlighted.gif"
        receptor, lig_name = pair.split('_', 1)
        title = f"{receptor} + {lig_name}"
        
        if create_gif(protein, ligand, out, title):
            print(f"  -> Created {out}")
    else:
        print(f"{pair}: no CIF found")
