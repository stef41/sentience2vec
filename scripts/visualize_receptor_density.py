#!/usr/bin/env python3
"""
Visualize Neuroreceptor Density Distribution

Creates visualizations of receptor density data including:
1. Heatmaps of receptor-region densities
2. Bar charts for individual receptors
3. Brain region comparison across receptors
"""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import defaultdict

# Set style
plt.style.use('dark_background')
plt.rcParams['figure.facecolor'] = '#1a1a2e'
plt.rcParams['axes.facecolor'] = '#16213e'
plt.rcParams['axes.edgecolor'] = '#e94560'
plt.rcParams['axes.labelcolor'] = '#eee'
plt.rcParams['text.color'] = '#eee'
plt.rcParams['xtick.color'] = '#eee'
plt.rcParams['ytick.color'] = '#eee'
plt.rcParams['grid.color'] = '#333'

# Receptor family colors
RECEPTOR_COLORS = {
    '5-HT1A': '#8b5cf6', '5-HT2A': '#a78bfa', '5-HT2B': '#c4b5fd', '5-HT2C': '#ddd6fe',
    'D1': '#f97316', 'D2': '#fb923c', 'D3': '#fdba74',
    'MOR': '#ec4899', 'KOR': '#f472b6', 'DOR': '#f9a8d4',
    'CB1': '#10b981',
    'NMDA_GluN2A': '#06b6d4', 'NMDA_GluN2B': '#22d3ee', 'AMPA_GluA1': '#67e8f9',
    'GABA_A_alpha1': '#eab308', 'GABA_A_alpha2': '#facc15',
    'alpha2A': '#ef4444',
    'Sigma1': '#84cc16',
    'nAChR_alpha4beta2': '#6366f1',
    'M1_muscarinic': '#14b8a6'
}

# Key brain regions for consciousness
CONSCIOUSNESS_REGIONS = [
    'prefrontal_cortex_layer_v', 'prefrontal_cortex_layer_ii',
    'anterior_cingulate_ba24', 'posterior_cingulate_ba23',
    'claustrum', 'insula_anterior',
    'thalamus', 'thalamus_pulvinar', 'thalamus_medial',
    'hippocampus_ca1', 'hippocampus_ca3',
    'visual_cortex_v1_layer_iv', 'visual_cortex_v2',
    'nucleus_accumbens', 'nucleus_accumbens_shell',
    'amygdala_basolateral', 'amygdala_central',
    'dorsal_raphe_nucleus', 'locus_coeruleus', 'vta'
]


def load_density_data(data_dir: Path) -> dict:
    """Load receptor density data."""
    with open(data_dir / 'receptor_densities' / 'receptor_density_distribution.json') as f:
        return json.load(f)


def plot_receptor_comparison(data: dict, output_dir: Path):
    """Create comparison of key receptors in consciousness-related regions."""
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Select key psychedelic receptors
    key_receptors = ['5-HT2A', '5-HT1A', 'D2', 'CB1', 'KOR', 'NMDA_GluN2A']
    
    # Build matrix
    regions_found = set()
    for receptor in key_receptors:
        if receptor in data['receptors']:
            for entry in data['receptors'][receptor]['density_distribution']:
                regions_found.add(entry['region_id'])
    
    # Filter to consciousness-related regions
    plot_regions = [r for r in CONSCIOUSNESS_REGIONS if r in regions_found][:15]
    
    # Create density matrix
    matrix = np.zeros((len(key_receptors), len(plot_regions)))
    
    for i, receptor in enumerate(key_receptors):
        if receptor in data['receptors']:
            region_densities = {e['region_id']: e['density_fmol_mg'] 
                               for e in data['receptors'][receptor]['density_distribution']}
            for j, region in enumerate(plot_regions):
                matrix[i, j] = region_densities.get(region, 0)
    
    # Normalize by row (receptor) for visualization
    matrix_norm = matrix / (matrix.max(axis=1, keepdims=True) + 1e-10)
    
    # Plot heatmap
    im = ax.imshow(matrix_norm, aspect='auto', cmap='magma')
    
    ax.set_xticks(range(len(plot_regions)))
    ax.set_xticklabels([r.replace('_', '\n') for r in plot_regions], 
                       rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(key_receptors)))
    ax.set_yticklabels(key_receptors, fontsize=11)
    
    # Add actual values as text
    for i in range(len(key_receptors)):
        for j in range(len(plot_regions)):
            if matrix[i, j] > 0:
                color = 'white' if matrix_norm[i, j] > 0.5 else 'gray'
                ax.text(j, i, f'{int(matrix[i, j])}', ha='center', va='center', 
                       fontsize=8, color=color)
    
    ax.set_title('Receptor Density in Consciousness-Related Brain Regions\n(fmol/mg tissue)',
                 fontsize=14, pad=20)
    ax.set_xlabel('Brain Region', fontsize=12, labelpad=10)
    ax.set_ylabel('Receptor', fontsize=12, labelpad=10)
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Relative Density', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'receptor_density_heatmap.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'receptor_density_heatmap.svg', bbox_inches='tight')
    plt.close()
    print(f"✓ Saved heatmap: {output_dir / 'receptor_density_heatmap.png'}")


def plot_individual_receptor(data: dict, receptor: str, output_dir: Path):
    """Create bar chart for a single receptor."""
    if receptor not in data['receptors']:
        print(f"  Receptor {receptor} not found")
        return
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    regions = data['receptors'][receptor]['density_distribution']
    names = [r['region_name'].replace('_', ' ').title() for r in regions[:15]]
    densities = [r['density_fmol_mg'] for r in regions[:15]]
    
    colors = [RECEPTOR_COLORS.get(receptor, '#888')] * len(names)
    bars = ax.barh(range(len(names)), densities, color=colors, alpha=0.8, edgecolor='white')
    
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    ax.invert_yaxis()
    
    ax.set_xlabel('Density (fmol/mg tissue)', fontsize=12)
    ax.set_title(f'{receptor} Receptor Distribution in Brain', fontsize=14, pad=15)
    
    # Add values on bars
    for bar, density in zip(bars, densities):
        ax.text(bar.get_width() + max(densities)*0.02, bar.get_y() + bar.get_height()/2,
                f'{density}', va='center', fontsize=9, color='white')
    
    ax.set_xlim(0, max(densities) * 1.15)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    filename = f'receptor_{receptor.replace("-", "").lower()}_distribution.png'
    plt.savefig(output_dir / filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {output_dir / filename}")


def plot_region_receptor_profile(data: dict, region: str, output_dir: Path):
    """Show all receptors present in a brain region."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    receptor_densities = []
    for receptor, info in data['receptors'].items():
        for entry in info['density_distribution']:
            if entry['region_id'] == region:
                receptor_densities.append({
                    'receptor': receptor,
                    'density': entry['density_fmol_mg']
                })
                break
    
    if not receptor_densities:
        print(f"  No data for region: {region}")
        return
    
    # Sort by density
    receptor_densities.sort(key=lambda x: x['density'], reverse=True)
    
    names = [r['receptor'] for r in receptor_densities]
    densities = [r['density'] for r in receptor_densities]
    colors = [RECEPTOR_COLORS.get(r['receptor'], '#888') for r in receptor_densities]
    
    bars = ax.barh(range(len(names)), densities, color=colors, alpha=0.8, edgecolor='white')
    
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    ax.invert_yaxis()
    
    ax.set_xlabel('Density (fmol/mg tissue)', fontsize=12)
    ax.set_title(f'Receptor Profile: {region.replace("_", " ").title()}', fontsize=14, pad=15)
    
    for bar, density in zip(bars, densities):
        ax.text(bar.get_width() + max(densities)*0.02, bar.get_y() + bar.get_height()/2,
                f'{density}', va='center', fontsize=9, color='white')
    
    ax.set_xlim(0, max(densities) * 1.15)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    filename = f'region_{region}_profile.png'
    plt.savefig(output_dir / filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {output_dir / filename}")


def plot_psychedelic_targets(data: dict, output_dir: Path):
    """Focused visualization of psychedelic receptor targets."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    targets = [
        ('5-HT2A', 'Primary Psychedelic Target'),
        ('5-HT1A', 'Anxiolytic/Introspection'),
        ('KOR', 'Salvinorin A Target'),
        ('CB1', 'Cannabinoid Target')
    ]
    
    for ax, (receptor, title) in zip(axes.flat, targets):
        if receptor not in data['receptors']:
            continue
            
        regions = data['receptors'][receptor]['density_distribution'][:10]
        names = [r['region_name'].replace('_', ' ').title()[:20] for r in regions]
        densities = [r['density_fmol_mg'] for r in regions]
        
        color = RECEPTOR_COLORS.get(receptor, '#888')
        bars = ax.barh(range(len(names)), densities, color=color, alpha=0.8, edgecolor='white')
        
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=9)
        ax.invert_yaxis()
        
        ax.set_title(f'{receptor}: {title}', fontsize=12, pad=10, color=color)
        ax.set_xlabel('Density (fmol/mg)', fontsize=10)
        ax.grid(axis='x', alpha=0.3)
        
        for bar, density in zip(bars, densities):
            ax.text(bar.get_width() + max(densities)*0.02, bar.get_y() + bar.get_height()/2,
                    f'{density}', va='center', fontsize=8, color='white')
        ax.set_xlim(0, max(densities) * 1.15)
    
    plt.suptitle('Psychedelic Receptor Target Distribution', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'psychedelic_targets.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'psychedelic_targets.svg', bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {output_dir / 'psychedelic_targets.png'}")


def main():
    """Generate all visualizations."""
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data'
    output_dir = data_dir / 'receptor_densities' / 'visualizations'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("GENERATING RECEPTOR DENSITY VISUALIZATIONS")
    print("="*60)
    
    # Load data
    print("\nLoading data...")
    data = load_density_data(data_dir)
    
    # Generate plots
    print("\nGenerating visualizations...")
    
    # Main comparison heatmap
    plot_receptor_comparison(data, output_dir)
    
    # Individual receptor plots
    key_receptors = ['5-HT2A', '5-HT1A', 'D2', 'CB1', 'MOR', 'KOR']
    for receptor in key_receptors:
        plot_individual_receptor(data, receptor, output_dir)
    
    # Region profiles
    key_regions = ['claustrum', 'prefrontal_cortex_layer_v', 'hippocampus_ca1', 'nucleus_accumbens']
    for region in key_regions:
        plot_region_receptor_profile(data, region, output_dir)
    
    # Psychedelic targets summary
    plot_psychedelic_targets(data, output_dir)
    
    print("\n" + "="*60)
    print("VISUALIZATION COMPLETE")
    print("="*60)
    print(f"\nOutput directory: {output_dir}")


if __name__ == '__main__':
    main()
