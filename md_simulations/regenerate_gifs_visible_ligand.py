#!/usr/bin/env python3
"""
Regenerate GIFs with HIGHLY VISIBLE ligand molecules.
- Ligand: Large spheres, bright neon colors, glowing effect
- Protein: Translucent, muted colors, smaller dots
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
from PIL import Image, ImageFilter, ImageEnhance
import io
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Protein: very muted, translucent
PROTEIN_COLORS = {
    'C': '#444466', 'N': '#445588', 'O': '#664444', 'S': '#666644',
    'H': '#555555', 'P': '#665544', 'CL': '#446644', 'BR': '#554444',
    'F': '#556644', 'I': '#554466',
}

# Ligand: SUPER BRIGHT neon colors
LIGAND_COLORS = {
    'C': '#00FF00',   # Neon green
    'N': '#00FFFF',   # Neon cyan
    'O': '#FF0066',   # Neon pink/red
    'S': '#FFFF00',   # Neon yellow
    'H': '#FFFFFF',   # White
    'P': '#FF8800',   # Neon orange
    'CL': '#66FF66',  # Light green
    'BR': '#FF6666',  # Light red
    'F': '#88FF88',   # Pale green
    'I': '#FF66FF',   # Neon magenta
}

# Standard amino acids (for identifying what's protein vs ligand)
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HOH', 'WAT', 'NA', 'CL', 'MG', 'ZN', 'CA', 'K', 'ACE', 'NME', 'NH2'
}

def parse_pdb_with_ligand(pdb_file):
    """Parse PDB file, separating protein from ligand atoms."""
    protein_atoms = []
    ligand_atoms = []
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip().upper()
            if not element:
                element = line[12:14].strip()[0].upper()
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip() or 'A'
            
            # Skip hydrogens for cleaner visualization
            if element == 'H':
                continue
            
            # Identify ligand: HETATM + non-standard residue, or chain B
            is_hetatm = line.startswith('HETATM')
            is_non_standard = resname not in STANDARD_RESIDUES
            is_lig_chain = chain in ['B', 'L', 'X']
            is_lig_residue = 'LIG' in resname or resname == 'UNK' or resname == 'UNL'
            
            if is_lig_residue or is_lig_chain or (is_hetatm and is_non_standard):
                ligand_atoms.append((x, y, z, element))
            else:
                protein_atoms.append((x, y, z, element))
        except:
            continue
    
    return protein_atoms, ligand_atoms

def parse_cif_with_ligand(cif_file):
    """Parse CIF file, separating protein from ligand atoms."""
    protein_atoms = []
    ligand_atoms = []
    
    col_map = {}
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('_atom_site.'):
            col_idx = 0
            while i < len(lines) and lines[i].strip().startswith('_atom_site.'):
                col_name = lines[i].strip().replace('_atom_site.', '')
                col_map[col_name] = col_idx
                col_idx += 1
                i += 1
            continue
        
        if col_map and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop') and line:
            parts = line.split()
            try:
                x = float(parts[col_map.get('Cartn_x', 0)])
                y = float(parts[col_map.get('Cartn_y', 1)])
                z = float(parts[col_map.get('Cartn_z', 2)])
                elem = parts[col_map.get('type_symbol', 3)].upper() if 'type_symbol' in col_map else 'C'
                comp = parts[col_map.get('label_comp_id', 4)].upper() if 'label_comp_id' in col_map else ''
                chain = parts[col_map.get('label_asym_id', 5)] if 'label_asym_id' in col_map else 'A'
                record_type = parts[0] if parts else ''  # ATOM or HETATM
                
                if elem == 'H':
                    i += 1
                    continue
                
                # Ligand detection: HETATM, LIG in comp name, chain B, or non-standard residue
                is_hetatm = record_type == 'HETATM'
                is_ligand = ('LIG' in comp or comp == 'UNK' or comp == 'UNL' or 
                            chain == 'B' or is_hetatm or comp not in STANDARD_RESIDUES)
                
                if is_ligand:
                    ligand_atoms.append((x, y, z, elem))
                else:
                    protein_atoms.append((x, y, z, elem))
            except:
                pass
        i += 1
    
    return protein_atoms, ligand_atoms

def create_visible_ligand_gif(protein_atoms, ligand_atoms, output_path, n_frames=30, dpi=100):
    """Create GIF with HIGHLY VISIBLE ligand."""
    if not ligand_atoms:
        # If no ligand detected, treat last 100 atoms as ligand (fallback)
        if len(protein_atoms) > 100:
            ligand_atoms = protein_atoms[-50:]
            protein_atoms = protein_atoms[:-50]
        else:
            return False
    
    all_atoms = protein_atoms + ligand_atoms
    if not all_atoms:
        return False
    
    all_coords = np.array([(a[0], a[1], a[2]) for a in all_atoms])
    center = all_coords.mean(axis=0)
    
    # Process protein
    if protein_atoms:
        p_coords = np.array([(a[0], a[1], a[2]) for a in protein_atoms]) - center
        p_elements = [a[3] for a in protein_atoms]
        # Subsample heavily for performance
        if len(protein_atoms) > 800:
            step = len(protein_atoms) // 800
            p_coords = p_coords[::step]
            p_elements = p_elements[::step]
        p_colors = [PROTEIN_COLORS.get(e, '#333344') for e in p_elements]
    else:
        p_coords = np.array([]).reshape(0, 3)
        p_colors = []
    
    # Process ligand (keep all atoms)
    l_coords = np.array([(a[0], a[1], a[2]) for a in ligand_atoms]) - center
    l_elements = [a[3] for a in ligand_atoms]
    l_colors = [LIGAND_COLORS.get(e, '#00FF00') for e in l_elements]
    
    frames = []
    
    for i in range(n_frames):
        fig = plt.figure(figsize=(5, 5), facecolor='#0a0a12')
        ax = fig.add_subplot(111, projection='3d', facecolor='#0a0a12')
        ax.view_init(elev=15, azim=i * (360 / n_frames))
        
        # Draw protein: small, translucent, muted
        if len(p_coords) > 0:
            ax.scatter(p_coords[:, 0], p_coords[:, 1], p_coords[:, 2],
                       c=p_colors, s=6, alpha=0.25, edgecolors='none', depthshade=True)
        
        # Draw ligand: LARGE, BRIGHT, with white edge for "glow"
        if len(l_coords) > 0:
            # Glow effect - larger pale version behind
            ax.scatter(l_coords[:, 0], l_coords[:, 1], l_coords[:, 2],
                       c=l_colors, s=200, alpha=0.3, edgecolors='none')
            # Main ligand spheres
            ax.scatter(l_coords[:, 0], l_coords[:, 1], l_coords[:, 2],
                       c=l_colors, s=120, alpha=1.0, edgecolors='white', linewidths=1)
        
        # Set limits
        margin = 3
        ax.set_xlim([all_coords[:, 0].min() - center[0] - margin, all_coords[:, 0].max() - center[0] + margin])
        ax.set_ylim([all_coords[:, 1].min() - center[1] - margin, all_coords[:, 1].max() - center[1] + margin])
        ax.set_zlim([all_coords[:, 2].min() - center[2] - margin, all_coords[:, 2].max() - center[2] + margin])
        ax.set_axis_off()
        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=dpi, facecolor='#0a0a12', bbox_inches='tight', pad_inches=0.05)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)
    
    # Save GIF
    frames[0].save(output_path, save_all=True, append_images=frames[1:], 
                   duration=80, loop=0, optimize=True)
    return True

def process_pair(args):
    """Worker function for parallel processing."""
    pair_dir, gif_path, force = args
    pair_name = pair_dir.name
    
    if gif_path.exists() and not force:
        return (pair_name, 'skipped')
    
    # Prioritize CIF (has proper ligand labels) over PDB
    cif = pair_dir / 'chai1' / 'pred.model_idx_0.cif'
    pdb = pair_dir / 'md' / 'minimized.pdb'
    
    protein_atoms, ligand_atoms = [], []
    
    if cif.exists():
        protein_atoms, ligand_atoms = parse_cif_with_ligand(cif)
    elif pdb.exists():
        protein_atoms, ligand_atoms = parse_pdb_with_ligand(pdb)
    
    if protein_atoms or ligand_atoms:
        try:
            if create_visible_ligand_gif(protein_atoms, ligand_atoms, gif_path):
                return (pair_name, 'created', len(ligand_atoms))
        except Exception as e:
            return (pair_name, f'error: {e}', 0)
    
    return (pair_name, 'no_structure', 0)

def main():
    output_base = Path('all_pairs_results')
    gif_dir = output_base / 'gifs'
    gif_dir.mkdir(exist_ok=True)
    
    force_regenerate = '--force' in sys.argv
    
    all_pairs = sorted([d for d in output_base.iterdir() if d.is_dir() and d.name != 'gifs'])
    total = len(all_pairs)
    
    print(f"🎬 Regenerating {total} GIFs with VISIBLE LIGANDS")
    print(f"   Force regenerate: {force_regenerate}")
    print("="*60)
    
    # Prepare work items
    work_items = [(p, gif_dir / f"{p.name}.gif", force_regenerate) for p in all_pairs]
    
    created = skipped = errors = 0
    
    # Use multiprocessing
    n_workers = max(1, multiprocessing.cpu_count() - 2)
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_pair, item): item for item in work_items}
        
        for i, future in enumerate(as_completed(futures)):
            result = future.result()
            name, status = result[0], result[1]
            
            if status == 'skipped':
                skipped += 1
            elif 'created' in status:
                created += 1
                lig_count = result[2] if len(result) > 2 else 0
                print(f"✓ {name} ({lig_count} ligand atoms)")
            else:
                errors += 1
                print(f"✗ {name}: {status}")
            
            if (i + 1) % 100 == 0:
                print(f"Progress: {i+1}/{total} ({created} created, {skipped} skipped)")
    
    print("="*60)
    print(f"Done! Created: {created}, Skipped: {skipped}, Errors: {errors}")

if __name__ == '__main__':
    main()
