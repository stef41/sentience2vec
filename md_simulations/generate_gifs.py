#!/usr/bin/env python3
"""
Generate rotating 3D molecular structure GIFs for key receptor-ligand pairs.
Uses matplotlib 3D for visualization.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import imageio
from PIL import Image
import io
import os

# Color scheme for atoms
ATOM_COLORS = {
    'C': '#909090',  # Gray
    'N': '#3050F8',  # Blue
    'O': '#FF0D0D',  # Red
    'S': '#FFFF30',  # Yellow
    'H': '#FFFFFF',  # White
    'P': '#FF8000',  # Orange
    'CL': '#1FF01F', # Green
    'BR': '#A62929', # Brown
    'F': '#90E050',  # Light green
    'I': '#940094',  # Purple
}

ATOM_RADII = {
    'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2,
    'P': 1.8, 'CL': 1.75, 'BR': 1.85, 'F': 1.47, 'I': 1.98
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
    in_atom_site = False
    
    # Column indices
    x_col, y_col, z_col, elem_col = None, None, None, None
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('_atom_site.'):
            # Find column positions
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
            in_atom_site = True
            continue
        
        if in_atom_site and line and not line.startswith('_') and not line.startswith('#'):
            if line.startswith('loop_') or line.startswith('data_'):
                in_atom_site = False
            else:
                parts = line.split()
                if len(parts) > max(x_col or 0, y_col or 0, z_col or 0, elem_col or 0):
                    try:
                        x = float(parts[x_col])
                        y = float(parts[y_col])
                        z = float(parts[z_col])
                        elem = parts[elem_col].upper() if elem_col else 'C'
                        atoms.append((x, y, z, elem))
                    except:
                        pass
        i += 1
    
    return atoms

def center_atoms(atoms):
    """Center atoms at origin."""
    if not atoms:
        return atoms
    coords = np.array([(a[0], a[1], a[2]) for a in atoms])
    center = coords.mean(axis=0)
    return [(a[0]-center[0], a[1]-center[1], a[2]-center[2], a[3]) for a in atoms]

def create_rotation_gif(atoms, output_path, title="", n_frames=36, figsize=(8, 8), dpi=100):
    """Create a rotating 3D GIF of the molecular structure."""
    if not atoms:
        print(f"No atoms to visualize for {output_path}")
        return False
    
    atoms = center_atoms(atoms)
    
    # Separate by element for coloring
    coords = np.array([(a[0], a[1], a[2]) for a in atoms])
    elements = [a[3] for a in atoms]
    
    # Subsample if too many atoms (for speed)
    if len(atoms) > 5000:
        indices = np.random.choice(len(atoms), 5000, replace=False)
        coords = coords[indices]
        elements = [elements[i] for i in indices]
    
    frames = []
    
    for frame in range(n_frames):
        fig = plt.figure(figsize=figsize, facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        
        # Plot atoms by element
        for elem in set(elements):
            mask = [e == elem for e in elements]
            elem_coords = coords[mask]
            color = ATOM_COLORS.get(elem, '#FFFFFF')
            size = ATOM_RADII.get(elem, 1.5) * 15
            ax.scatter(elem_coords[:, 0], elem_coords[:, 1], elem_coords[:, 2],
                      c=color, s=size, alpha=0.8, edgecolors='none')
        
        # Rotate view
        angle = frame * (360 / n_frames)
        ax.view_init(elev=20, azim=angle)
        
        # Style
        ax.set_xlim(coords[:, 0].min() - 5, coords[:, 0].max() + 5)
        ax.set_ylim(coords[:, 1].min() - 5, coords[:, 1].max() + 5)
        ax.set_zlim(coords[:, 2].min() - 5, coords[:, 2].max() + 5)
        
        ax.set_axis_off()
        ax.set_title(title, color='white', fontsize=14, fontweight='bold', pad=20)
        
        # Save frame to buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=dpi, facecolor='black', 
                    bbox_inches='tight', pad_inches=0.1)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)
    
    # Save GIF
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=100,  # ms per frame
        loop=0
    )
    
    print(f"✓ Created: {output_path}")
    return True

def main():
    output_base = Path('all_pairs_results')
    gif_dir = output_base / 'gifs'
    gif_dir.mkdir(exist_ok=True)
    
    # Key pairs to visualize
    key_pairs = [
        '5-HT2A_DMT',
        '5-HT2A_LSD', 
        '5-HT2A_Psilocin',
        '5-HT2A_Mescaline',
        '5-HT2A_MDMA',
        '5-HT2A_2C-B',
        '5-HT1A_DMT',
        'D2_LSD',
        'SERT_MDMA',
        'CB1_THC-O-acetate',
        'MOR_Morphine',
        'NMDA_Ketamine',
    ]
    
    print("="*60)
    print("GENERATING 3D MOLECULAR STRUCTURE GIFS")
    print("="*60)
    
    created = 0
    
    for pair_name in key_pairs:
        pair_dir = output_base / pair_name
        
        # Try minimized PDB first, then CIF
        pdb_file = pair_dir / 'md' / 'minimized.pdb'
        cif_file = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
        
        atoms = None
        source = None
        
        if pdb_file.exists():
            atoms = parse_pdb(pdb_file)
            source = 'minimized'
        elif cif_file.exists():
            atoms = parse_cif(cif_file)
            source = 'chai1'
        
        if atoms:
            receptor, ligand = pair_name.split('_', 1)
            title = f"{receptor} + {ligand}\n({source}, {len(atoms)} atoms)"
            gif_path = gif_dir / f"{pair_name}.gif"
            
            print(f"\nProcessing {pair_name} ({len(atoms)} atoms)...")
            if create_rotation_gif(atoms, gif_path, title=title):
                created += 1
        else:
            print(f"✗ No structure found for {pair_name}")
    
    print(f"\n{'='*60}")
    print(f"Created {created} GIFs in {gif_dir}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
