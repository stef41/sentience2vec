"""
Visualization utilities for consciousness research
Creates plots of receptor density maps and substance-region activation
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Optional
from pathlib import Path


def plot_receptor_density_heatmap(
    receptor_data: Dict,
    receptors: List[str] = None,
    output_path: Optional[str] = None,
    figsize: tuple = (14, 10)
):
    """
    Create heatmap of receptor densities across brain regions
    """
    from src.pharmacology.receptor_mapping import RECEPTOR_BRAIN_REGIONS
    
    if receptors is None:
        receptors = ['5-HT1A', '5-HT2A', '5-HT2C', 'D1', 'D2', 'NMDA_GluN2B', 
                     'GABA_A_alpha1', 'MOR', 'KOR', 'CB1']
    
    # Collect all regions across selected receptors
    all_regions = set()
    for receptor in receptors:
        info = RECEPTOR_BRAIN_REGIONS.get(receptor, {})
        density_map = info.get('density_map', {})
        all_regions.update(density_map.keys())
    
    all_regions = sorted(all_regions)
    
    # Build matrix
    matrix = np.zeros((len(receptors), len(all_regions)))
    
    for i, receptor in enumerate(receptors):
        info = RECEPTOR_BRAIN_REGIONS.get(receptor, {})
        density_map = info.get('density_map', {})
        for j, region in enumerate(all_regions):
            if region in density_map:
                matrix[i, j] = density_map[region].get('density_fmol_mg', 0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Log scale for better visualization
    matrix_log = np.log10(matrix + 1)
    
    im = ax.imshow(matrix_log, cmap='YlOrRd', aspect='auto')
    
    # Labels
    ax.set_xticks(range(len(all_regions)))
    ax.set_xticklabels([r.replace('_', ' ')[:25] for r in all_regions], 
                       rotation=90, fontsize=6)
    ax.set_yticks(range(len(receptors)))
    ax.set_yticklabels(receptors, fontsize=10)
    
    ax.set_xlabel('Brain Region', fontsize=12)
    ax.set_ylabel('Receptor', fontsize=12)
    ax.set_title('Receptor Density Map (log10 fmol/mg)', fontsize=14)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log10(density + 1)', fontsize=10)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_path}")
    
    return fig, ax


def plot_substance_activation_profile(
    substance: str,
    output_path: Optional[str] = None,
    top_n: int = 20,
    figsize: tuple = (12, 8)
):
    """
    Plot brain region activation profile for a substance
    """
    from src.pharmacology.receptor_mapping import PharmacologyAnalyzer
    
    analyzer = PharmacologyAnalyzer()
    regions = analyzer.get_brain_regions_for_substance(substance)
    
    if not regions:
        print(f"No data for {substance}")
        return None, None
    
    # Get top N regions
    top_regions = list(regions.items())[:top_n]
    
    region_names = [r[0].replace('_', ' ') for r in top_regions]
    activations = [r[1]['total_activation'] for r in top_regions]
    
    # Color by primary receptor
    colors = []
    for r in top_regions:
        receptors = [x['receptor'] for x in r[1]['receptors']]
        if '5-HT2A' in receptors:
            colors.append('#FF6B6B')  # Red - psychedelic
        elif '5-HT1A' in receptors:
            colors.append('#4ECDC4')  # Teal
        elif 'D2' in receptors or 'D1' in receptors:
            colors.append('#45B7D1')  # Blue - dopamine
        elif 'NMDA' in str(receptors):
            colors.append('#96CEB4')  # Green - dissociative
        elif 'KOR' in receptors:
            colors.append('#9B59B6')  # Purple - kappa
        elif 'CB1' in receptors:
            colors.append('#F39C12')  # Orange - cannabinoid
        else:
            colors.append('#95A5A6')  # Gray
    
    fig, ax = plt.subplots(figsize=figsize)
    
    bars = ax.barh(range(len(region_names)), activations, color=colors)
    ax.set_yticks(range(len(region_names)))
    ax.set_yticklabels(region_names, fontsize=9)
    ax.invert_yaxis()
    
    ax.set_xlabel('Activation Score (density × affinity)', fontsize=12)
    ax.set_title(f'{substance} Brain Region Activation Profile', fontsize=14)
    
    # Add receptor labels
    for i, (bar, region_data) in enumerate(zip(bars, top_regions)):
        receptors = ', '.join([x['receptor'] for x in region_data[1]['receptors']])
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                receptors, va='center', fontsize=7, alpha=0.7)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_path}")
    
    return fig, ax


def plot_effect_network(
    substance: str,
    output_path: Optional[str] = None,
    figsize: tuple = (14, 10)
):
    """
    Create network visualization of substance -> receptor -> effect relationships
    """
    import networkx as nx
    from src.pharmacology.receptor_mapping import PharmacologyAnalyzer, RECEPTOR_BRAIN_REGIONS
    
    analyzer = PharmacologyAnalyzer()
    profile = analyzer.get_receptor_profile(substance)
    
    if not profile:
        print(f"No data for {substance}")
        return None, None
    
    G = nx.DiGraph()
    
    # Add substance node
    G.add_node(substance, node_type='substance', size=3000)
    
    # Add receptor and effect nodes
    for binding in profile.bindings:
        receptor = binding.receptor
        ki = binding.ki_nm
        
        # Size based on affinity
        receptor_size = max(500, 2000 / (np.log10(ki + 1) + 1))
        G.add_node(receptor, node_type='receptor', size=receptor_size, ki=ki)
        G.add_edge(substance, receptor, weight=1/np.log10(ki+1))
        
        # Add effects
        receptor_info = RECEPTOR_BRAIN_REGIONS.get(receptor, {})
        effects = receptor_info.get('associated_effects', [])
        
        for effect in effects[:4]:  # Limit effects per receptor
            G.add_node(effect, node_type='effect', size=400)
            G.add_edge(receptor, effect, weight=0.5)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # Node colors by type
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        ntype = G.nodes[node].get('node_type', 'unknown')
        if ntype == 'substance':
            node_colors.append('#E74C3C')
            node_sizes.append(3000)
        elif ntype == 'receptor':
            node_colors.append('#3498DB')
            node_sizes.append(G.nodes[node].get('size', 1000))
        else:
            node_colors.append('#2ECC71')
            node_sizes.append(400)
    
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                          alpha=0.9, ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.4, arrows=True, 
                          edge_color='gray', ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
    
    ax.set_title(f'{substance} Effect Network', fontsize=14)
    ax.axis('off')
    
    # Legend
    legend_elements = [
        plt.scatter([], [], c='#E74C3C', s=150, label='Substance'),
        plt.scatter([], [], c='#3498DB', s=100, label='Receptor'),
        plt.scatter([], [], c='#2ECC71', s=50, label='Effect')
    ]
    ax.legend(handles=legend_elements, loc='upper left')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_path}")
    
    return fig, ax


def generate_all_visualizations(output_dir: str = "./figures"):
    """Generate all standard visualizations"""
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating visualizations...")
    
    # 1. Receptor density heatmap
    print("  - Receptor density heatmap")
    plot_receptor_density_heatmap(
        {},
        output_path=str(output_dir / "receptor_density_heatmap.png")
    )
    
    # 2. Substance profiles
    substances = ['LSD', 'DMT', 'Ketamine', 'Salvinorin_A', 'Cannabis_THC']
    for sub in substances:
        print(f"  - {sub} activation profile")
        plot_substance_activation_profile(
            sub,
            output_path=str(output_dir / f"{sub.lower()}_activation.png")
        )
    
    # 3. Effect networks
    for sub in ['LSD', 'Ketamine']:
        print(f"  - {sub} effect network")
        plot_effect_network(
            sub,
            output_path=str(output_dir / f"{sub.lower()}_network.png")
        )
    
    print(f"\nAll visualizations saved to {output_dir}")


if __name__ == "__main__":
    # Quick test
    print("Testing visualization module...")
    
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        
        fig, ax = plot_substance_activation_profile('LSD', top_n=15)
        if fig:
            fig.savefig('/tmp/lsd_test.png', dpi=100)
            print("Test visualization saved to /tmp/lsd_test.png")
        
    except ImportError as e:
        print(f"Matplotlib not available: {e}")
        print("Install with: pip install matplotlib")
