#!/usr/bin/env python3
"""
Generate rotating 3D molecular structure GIFs with highlighted ligand.
Ligand atoms are shown larger and in bright colors.
"""
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
import sys

# Colors for protein atoms (visible but not overwhelming)
PROTEIN_COLORS = {
    'C': '#8888AA', 'N': '#4060C0', 'O': '#C04040', 'S': '#C0C040',
    'H': '#AAAAAA', 'P': '#C08040', 'CL': '#40C040', 'BR': '#905050',
    'F': '#80C060', 'I': '#804080',
}

# Colors for ligand atoms (bright)
LIGAND_COLORS = {
    'C': '#00FF00',  # Bright green
    'N': '#00CCFF',  # Cyan
    'O': '#FF3333',  # Bright red
    'S': '#FFFF00',  # Yellow
    'H': '#FFFFFF',  # White
    'P': '#FF8800',  # Orange
    'CL': '#33FF33', 'BR': '#FF6666', 'F': '#AAFFAA', 'I': '#FF00FF',
}

def parse_pdb_with_ligand(pdb_file):
    """Parse PDB and identify ligand atoms (last chain or HETATM)."""
    protein_atoms = []
    ligand_atoms = []
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # First pass: find all chains
    chains = set()
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21:22].strip() or 'A'
            chains.add(chain)
    
    # Assume last chain alphabetically is the ligand, or chain B
    ligand_chain = max(chains) if len(chains) > 1 else None
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip()
                if not element:
                    element = line[12:14].strip()[0]
                element = element.upper()
                chain = line[21:22].strip() or 'A'
                resname = line[17:20].strip()
                
                is_hetatm = line.startswith('HETATM')
                is_ligand_chain = (chain == ligand_chain) if ligand_chain else False
                is_non_standard = resname not in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HOH','WAT']
                
                if is_hetatm or is_ligand_chain or is_non_standard:
                    ligand_atoms.append((x, y, z, element))
                else:
                    protein_atoms.append((x, y, z, element))
            except:
                continue
    
    return protein_atoms, ligand_atoms

def parse_cif_with_ligand(cif_file):
    """Parse CIF and identify ligand atoms (LIG residues or chain B)."""
    protein_atoms = []
    ligand_atoms = []
    
    x_col, y_col, z_col, elem_col, comp_col, chain_col = None, None, None, None, None, None
    
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
                elif 'label_asym_id' in col_name: chain_col = col_idx
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
                    chain = parts[chain_col] if chain_col and len(parts) > chain_col else 'A'
                    
                    # Ligand identification: LIG in residue name or chain B
                    is_ligand = 'LIG' in comp.upper() or chain == 'B'
                    
                    if is_ligand:
                        ligand_atoms.append((x, y, z, elem))
                    else:
                        protein_atoms.append((x, y, z, elem))
                except:
                    pass
        i += 1
    
    return protein_atoms, ligand_atoms

def create_gif(protein_atoms, ligand_atoms, output_path, n_frames=24, dpi=80):
    """Create rotating GIF with highlighted ligand."""
    if not protein_atoms and not ligand_atoms:
        return False
    
    all_atoms = protein_atoms + ligand_atoms
    all_coords = np.array([(a[0], a[1], a[2]) for a in all_atoms])
    center = all_coords.mean(axis=0)
    
    # Process protein atoms
    if protein_atoms:
        p_coords = np.array([(a[0], a[1], a[2]) for a in protein_atoms]) - center
        p_elements = [a[3] for a in protein_atoms]
        # Subsample protein if too large
        if len(protein_atoms) > 1500:
            step = len(protein_atoms) // 1500
            p_coords = p_coords[::step]
            p_elements = p_elements[::step]
        p_colors = [PROTEIN_COLORS.get(e, '#404040') for e in p_elements]
    else:
        p_coords = np.array([]).reshape(0, 3)
        p_colors = []
    
    # Process ligand atoms (keep all, no subsampling)
    if ligand_atoms:
        l_coords = np.array([(a[0], a[1], a[2]) for a in ligand_atoms]) - center
        l_elements = [a[3] for a in ligand_atoms]
        l_colors = [LIGAND_COLORS.get(e, '#00FF00') for e in l_elements]
    else:
        l_coords = np.array([]).reshape(0, 3)
        l_colors = []
    
    frames = []
    
    for i in range(n_frames):
        fig = plt.figure(figsize=(5, 5), facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        ax.view_init(elev=20, azim=i * (360 / n_frames))
        
        # Draw protein (larger, more visible)
        if len(p_coords) > 0:
            ax.scatter(p_coords[:, 0], p_coords[:, 1], p_coords[:, 2],
                       c=p_colors, s=12, alpha=0.6, edgecolors='none')
        
        # Draw ligand (large, bright, opaque)
        if len(l_coords) > 0:
            ax.scatter(l_coords[:, 0], l_coords[:, 1], l_coords[:, 2],
                       c=l_colors, s=80, alpha=1.0, edgecolors='white', linewidths=0.5)
        
        # Set axis limits based on all atoms
        if len(all_coords) > 0:
            all_centered = all_coords - center
            margin = 5
            ax.set_xlim([all_centered[:, 0].min() - margin, all_centered[:, 0].max() + margin])
            ax.set_ylim([all_centered[:, 1].min() - margin, all_centered[:, 1].max() + margin])
            ax.set_zlim([all_centered[:, 2].min() - margin, all_centered[:, 2].max() + margin])
        
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
    
    print(f"Processing {total} pairs (with ligand highlighting)...", flush=True)
    
    created = skipped = no_struct = 0
    
    for i, pair_dir in enumerate(all_pairs):
        gif_path = gif_dir / f"{pair_dir.name}.gif"
        
        # Always regenerate to get new ligand highlighting
        # if gif_path.exists():
        #     skipped += 1
        #     continue
        
        pdb = pair_dir / 'md' / 'minimized.pdb'
        cif = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
        
        protein_atoms, ligand_atoms = [], []
        
        if pdb.exists():
            protein_atoms, ligand_atoms = parse_pdb_with_ligand(pdb)
        elif cif.exists():
            protein_atoms, ligand_atoms = parse_cif_with_ligand(cif)
        
        if protein_atoms or ligand_atoms:
            try:
                if create_gif(protein_atoms, ligand_atoms, gif_path):
                    created += 1
                    if ligand_atoms:
                        print(f"  ✓ {pair_dir.name}: {len(protein_atoms)} protein + {len(ligand_atoms)} ligand atoms", flush=True)
            except Exception as e:
                print(f"  ✗ {pair_dir.name}: {e}", flush=True)
        else:
            no_struct += 1
        
        if (i + 1) % 100 == 0:
            print(f"[{i+1}/{total}] Created: {created}, No struct: {no_struct}", flush=True)
    
    print(f"\nDONE: Created {created}, No struct {no_struct}")
    print(f"Total GIFs: {len(list(gif_dir.glob('*.gif')))}")

if __name__ == '__main__':
    main()
