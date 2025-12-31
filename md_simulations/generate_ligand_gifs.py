#!/usr/bin/env python3
"""
Generate rotating 3D molecular structure GIFs with ligand highlighted.
Protein shown as transparent surface, ligand as bright colored spheres.
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

# Color scheme for ligand atoms (bright colors)
LIGAND_COLORS = {
    'C': '#00FF00',  # Bright green for carbon
    'N': '#0080FF',  # Bright blue
    'O': '#FF0000',  # Bright red
    'S': '#FFFF00',  # Yellow
    'H': '#FFFFFF',  # White
    'P': '#FF8000',  # Orange
    'CL': '#00FF80', # Cyan-green
    'BR': '#A62929', # Brown
    'F': '#90E050',  # Light green
}

# Colors for protein
PROTEIN_COLOR = '#6699DD'  # Light blue

def parse_cif_with_ligand(cif_file):
    """Parse CIF file and separate protein from ligand atoms."""
    protein_atoms = []
    ligand_atoms = []
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            parts = line.split()
            if len(parts) >= 12:
                try:
                    # CIF format: record type, id, atom_type, atom_name, ?, residue, resnum, ...
                    # coords are typically around positions 10, 11, 12
                    x = float(parts[10])
                    y = float(parts[11])
                    z = float(parts[12])
                    
                    # Get element (usually position 2 or from atom name)
                    element = parts[2].upper()
                    if len(element) > 2:
                        element = element[0]
                    
                    # Check if ligand (LIG in residue name or chain B)
                    residue = parts[5] if len(parts) > 5 else ''
                    chain = parts[17] if len(parts) > 17 else ''
                    
                    is_ligand = 'LIG' in residue or chain == 'B'
                    
                    if is_ligand:
                        ligand_atoms.append((x, y, z, element))
                    else:
                        protein_atoms.append((x, y, z, element))
                except (ValueError, IndexError):
                    continue
    
    return protein_atoms, ligand_atoms

def center_atoms(protein_atoms, ligand_atoms):
    """Center all atoms around the ligand."""
    if ligand_atoms:
        # Center on ligand
        ligand_coords = np.array([(a[0], a[1], a[2]) for a in ligand_atoms])
        center = ligand_coords.mean(axis=0)
    elif protein_atoms:
        # Fall back to protein center
        protein_coords = np.array([(a[0], a[1], a[2]) for a in protein_atoms])
        center = protein_coords.mean(axis=0)
    else:
        return [], []
    
    protein_centered = [(a[0]-center[0], a[1]-center[1], a[2]-center[2], a[3]) for a in protein_atoms]
    ligand_centered = [(a[0]-center[0], a[1]-center[1], a[2]-center[2], a[3]) for a in ligand_atoms]
    
    return protein_centered, ligand_centered

def create_rotation_gif(protein_atoms, ligand_atoms, output_path, title="", 
                        n_frames=36, figsize=(10, 10), dpi=100, zoom=30):
    """Create a rotating 3D GIF with ligand highlighted."""
    if not ligand_atoms:
        print(f"  No ligand atoms found!")
        return False
    
    protein_atoms, ligand_atoms = center_atoms(protein_atoms, ligand_atoms)
    
    # Get coordinates
    ligand_coords = np.array([(a[0], a[1], a[2]) for a in ligand_atoms])
    ligand_elements = [a[3] for a in ligand_atoms]
    
    # Subsample protein for speed (keep atoms near ligand)
    if protein_atoms:
        protein_coords = np.array([(a[0], a[1], a[2]) for a in protein_atoms])
        
        # Calculate distance from ligand center
        ligand_center = ligand_coords.mean(axis=0)
        distances = np.linalg.norm(protein_coords - ligand_center, axis=1)
        
        # Keep only atoms within zoom distance of ligand
        mask = distances < zoom
        protein_coords = protein_coords[mask]
        
        # Further subsample if still too many
        if len(protein_coords) > 3000:
            indices = np.random.choice(len(protein_coords), 3000, replace=False)
            protein_coords = protein_coords[indices]
    else:
        protein_coords = np.array([]).reshape(0, 3)
    
    frames = []
    
    for frame in range(n_frames):
        fig = plt.figure(figsize=figsize, facecolor='black')
        ax = fig.add_subplot(111, projection='3d', facecolor='black')
        
        # Plot protein atoms (visible, medium size)
        if len(protein_coords) > 0:
            ax.scatter(protein_coords[:, 0], protein_coords[:, 1], protein_coords[:, 2],
                      c=PROTEIN_COLOR, s=30, alpha=0.6, edgecolors='none')
        
        # Plot ligand atoms (bright, larger)
        for elem in set(ligand_elements):
            mask = [e == elem for e in ligand_elements]
            elem_coords = ligand_coords[mask]
            color = LIGAND_COLORS.get(elem, '#FF00FF')
            # Larger size for ligand atoms
            size = 200 if elem != 'H' else 80
            ax.scatter(elem_coords[:, 0], elem_coords[:, 1], elem_coords[:, 2],
                      c=color, s=size, alpha=1.0, edgecolors='white', linewidths=0.5)
        
        # Rotate view
        angle = frame * (360 / n_frames)
        ax.view_init(elev=15, azim=angle)
        
        # Set view limits centered on ligand
        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(-zoom, zoom)
        ax.set_zlim(-zoom, zoom)
        
        ax.set_axis_off()
        
        # Title with ligand info
        ax.set_title(title, color='white', fontsize=16, fontweight='bold', pad=20)
        
        # Add legend
        ax.text2D(0.02, 0.98, f'Ligand: {len(ligand_atoms)} atoms', 
                  transform=ax.transAxes, color='#00FF00', fontsize=10,
                  verticalalignment='top')
        ax.text2D(0.02, 0.94, f'Binding site: {len(protein_coords)} atoms', 
                  transform=ax.transAxes, color=PROTEIN_COLOR, fontsize=10,
                  verticalalignment='top')
        
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
        duration=100,
        loop=0
    )
    
    print(f"  ✓ Created: {output_path.name}")
    return True

def main():
    output_base = Path('all_pairs_results')
    gif_dir = output_base / 'gifs'
    gif_dir.mkdir(exist_ok=True)
    
    # Key pairs to visualize
    key_pairs = [
        ('5-HT2A_DMT', '5-HT2A + DMT'),
        ('5-HT2A_LSD', '5-HT2A + LSD'),
        ('5-HT2A_Psilocin', '5-HT2A + Psilocin'),
        ('5-HT2A_Mescaline', '5-HT2A + Mescaline'),
        ('5-HT2A_MDMA', '5-HT2A + MDMA'),
        ('5-HT2A_2C-B', '5-HT2A + 2C-B'),
        ('5-HT2A_5-MeO-DMT', '5-HT2A + 5-MeO-DMT'),
        ('5-HT1A_DMT', '5-HT1A + DMT'),
        ('D2_LSD', 'D2 Dopamine + LSD'),
        ('D2_DMT', 'D2 Dopamine + DMT'),
        ('SERT_MDMA', 'SERT + MDMA'),
        ('CB1_THC-O-acetate', 'CB1 + THC-O'),
        ('MOR_Morphine', 'μ-Opioid + Morphine'),
        ('GABA-A_Diazepam', 'GABA-A + Diazepam'),
    ]
    
    print("="*70)
    print("GENERATING 3D MOLECULAR GIFS WITH LIGAND HIGHLIGHTING")
    print("="*70)
    
    created = 0
    
    for pair_name, title in key_pairs:
        pair_dir = output_base / pair_name
        cif_file = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
        
        if not cif_file.exists():
            print(f"\n✗ {pair_name}: CIF not found")
            continue
        
        print(f"\nProcessing {pair_name}...")
        
        protein_atoms, ligand_atoms = parse_cif_with_ligand(cif_file)
        print(f"  Protein: {len(protein_atoms)} atoms, Ligand: {len(ligand_atoms)} atoms")
        
        if ligand_atoms:
            gif_path = gif_dir / f"{pair_name}_ligand.gif"
            if create_rotation_gif(protein_atoms, ligand_atoms, gif_path, title=title):
                created += 1
        else:
            print(f"  ✗ No ligand atoms found in CIF")
    
    print(f"\n{'='*70}")
    print(f"Created {created} GIFs in {gif_dir}")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
