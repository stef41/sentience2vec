#!/usr/bin/env python3
"""
Create dimensionality reduction plot of drugs based on receptor binding profiles.
"""
import re
import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# Try to import UMAP, fall back to t-SNE if not available
try:
    from umap import UMAP
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    print("UMAP not available, using t-SNE")

# Color scheme for drug classes
CLASS_COLORS = {
    'psychedelic': '#FF6B6B',      # Red
    'dissociative': '#4ECDC4',     # Cyan
    'entactogen': '#FFE66D',       # Yellow
    'empathogen': '#FFE66D',       # Yellow (same as entactogen)
    'cannabinoid': '#95E86C',      # Green
    'opioid': '#A855F7',           # Purple
    'stimulant': '#FF8C42',        # Orange
    'depressant': '#6B7280',       # Gray
    'deliriant': '#EC4899',        # Pink
    'nootropic': '#3B82F6',        # Blue
    'other': '#9CA3AF',            # Light gray
}

def parse_data_js(filepath):
    """Parse the data.js file to extract drug binding profiles."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find drugs section
    drugs_match = re.search(r'drugs:\s*\{', content)
    if not drugs_match:
        raise ValueError("Could not find drugs section")
    
    # Extract drugs manually
    drugs = {}
    
    # Find each drug entry
    drug_pattern = re.compile(r'"([^"]+)":\s*\{\s*"receptors":\s*\{([^}]+)\}[^}]*"class":\s*"([^"]+)"', re.DOTALL)
    
    for match in drug_pattern.finditer(content):
        drug_name = match.group(1)
        receptors_str = match.group(2)
        drug_class = match.group(3)
        
        # Parse receptors
        receptors = {}
        receptor_pattern = re.compile(r'"([^"]+)":\s*(-?[\d.]+)')
        for rec_match in receptor_pattern.finditer(receptors_str):
            receptors[rec_match.group(1)] = float(rec_match.group(2))
        
        if receptors:
            drugs[drug_name] = {
                'receptors': receptors,
                'class': drug_class
            }
    
    return drugs

def create_feature_matrix(drugs):
    """Create feature matrix from drug binding profiles."""
    # Get all unique receptors
    all_receptors = set()
    for drug in drugs.values():
        all_receptors.update(drug['receptors'].keys())
    all_receptors = sorted(all_receptors)
    
    # Create matrix
    drug_names = list(drugs.keys())
    X = np.zeros((len(drug_names), len(all_receptors)))
    
    for i, drug_name in enumerate(drug_names):
        for j, receptor in enumerate(all_receptors):
            X[i, j] = drugs[drug_name]['receptors'].get(receptor, 0)
    
    return X, drug_names, all_receptors

def main():
    print("Loading drug binding data...")
    drugs = parse_data_js('/data/users/zacharie/whc/consciousness_study/website/data.js')
    print(f"Found {len(drugs)} drugs with binding profiles")
    
    # Filter out "other" class to show meaningful clusters
    drugs_filtered = {name: data for name, data in drugs.items() if data['class'] != 'other'}
    print(f"After filtering 'other': {len(drugs_filtered)} drugs with known classes")
    
    # Only keep drugs with at least 5 receptor measurements (rich profiles)
    drugs_rich = {name: data for name, data in drugs_filtered.items() 
                  if len(data['receptors']) >= 5}
    print(f"After filtering sparse profiles: {len(drugs_rich)} drugs with ≥5 receptors")
    
    # Show class distribution
    from collections import Counter
    class_counts = Counter(d['class'] for d in drugs_rich.values())
    print("Class distribution:")
    for cls, count in class_counts.most_common():
        print(f"  {cls}: {count}")
    
    drugs = drugs_rich
    
    # Create feature matrix
    X, drug_names, receptors = create_feature_matrix(drugs)
    print(f"Feature matrix: {X.shape[0]} drugs x {X.shape[1]} receptors")
    
    # Normalize binding energies (more negative = stronger binding)
    # Convert to positive (higher = stronger)
    X = -X
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Apply dimensionality reduction
    print("Running dimensionality reduction...")
    
    # Use t-SNE with adjusted perplexity for smaller dataset
    perplexity = min(30, len(drug_names) // 4)
    reducer = TSNE(n_components=2, perplexity=perplexity, random_state=42, n_iter=2000)
    X_2d = reducer.fit_transform(X_scaled)
    method = "t-SNE"
    
    # Get drug classes
    classes = [drugs[name]['class'] for name in drug_names]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(16, 12), facecolor='#0a0a1a')
    ax.set_facecolor('#0a0a1a')
    
    # Plot each class
    unique_classes = list(set(classes))
    for drug_class in unique_classes:
        mask = [c == drug_class for c in classes]
        color = CLASS_COLORS.get(drug_class, '#9CA3AF')
        
        x = X_2d[mask, 0]
        y = X_2d[mask, 1]
        
        ax.scatter(x, y, c=color, label=drug_class.capitalize(), 
                   s=80, alpha=0.7, edgecolors='white', linewidths=0.5)
    
    # Add labels for key drugs
    key_drugs = ['LSD', 'DMT', 'Psilocybin', 'Mescaline', 'MDMA', 'Ketamine', 
                 'THC', 'Morphine', 'Amphetamine', 'Cocaine', '2C-B', 'Salvinorin_A',
                 'Ibogaine', 'Ayahuasca', '5-MeO-DMT', 'DXM', 'PCP', 'Nitrous']
    
    for i, name in enumerate(drug_names):
        if any(key in name for key in key_drugs):
            ax.annotate(name, (X_2d[i, 0], X_2d[i, 1]), 
                        fontsize=8, color='white', alpha=0.9,
                        xytext=(5, 5), textcoords='offset points')
    
    ax.set_title(f'Drug Receptor Binding Profile Similarity ({method})', 
                 fontsize=16, color='white', pad=20)
    ax.set_xlabel(f'{method} Dimension 1', fontsize=12, color='white')
    ax.set_ylabel(f'{method} Dimension 2', fontsize=12, color='white')
    
    # Style axes
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('#333')
    
    # Legend
    legend = ax.legend(loc='upper right', fontsize=10, framealpha=0.3)
    legend.get_frame().set_facecolor('#1a1a2e')
    for text in legend.get_texts():
        text.set_color('white')
    
    plt.tight_layout()
    
    # Save
    output_path = '/data/users/zacharie/whc/consciousness_study/website/drug_binding_map.png'
    plt.savefig(output_path, dpi=150, facecolor='#0a0a1a', bbox_inches='tight')
    print(f"Saved to {output_path}")
    
    # Also save interactive HTML version
    create_interactive_plot(X_2d, drug_names, classes, drugs, method)
    
    plt.close()

def create_interactive_plot(X_2d, drug_names, classes, drugs, method):
    """Create an interactive HTML version of the plot."""
    
    html_content = '''<!DOCTYPE html>
<html>
<head>
    <title>Drug Binding Profile Map</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body { 
            margin: 0; 
            background: #0a0a1a; 
            font-family: Arial, sans-serif;
            color: white;
        }
        #plot { width: 100vw; height: 100vh; }
        h1 {
            position: absolute;
            top: 10px;
            left: 20px;
            margin: 0;
            font-size: 24px;
            z-index: 100;
        }
    </style>
</head>
<body>
    <h1>Drug Receptor Binding Profile Similarity (''' + method + ''')</h1>
    <div id="plot"></div>
    <script>
        const classColors = {
            'psychedelic': '#FF6B6B',
            'dissociative': '#4ECDC4',
            'entactogen': '#FFE66D',
            'cannabinoid': '#95E86C',
            'opioid': '#A855F7',
            'stimulant': '#FF8C42',
            'depressant': '#6B7280',
            'deliriant': '#EC4899',
            'nootropic': '#3B82F6',
            'other': '#9CA3AF'
        };
        
        const data = ''' + json.dumps([
            {
                'name': drug_names[i],
                'class': classes[i],
                'x': float(X_2d[i, 0]),
                'y': float(X_2d[i, 1]),
                'receptors': drugs[drug_names[i]]['receptors']
            }
            for i in range(len(drug_names))
        ]) + ''';
        
        // Group by class
        const traces = {};
        data.forEach(d => {
            if (!traces[d.class]) {
                traces[d.class] = {
                    x: [], y: [], text: [], hovertext: [],
                    mode: 'markers',
                    type: 'scatter',
                    name: d.class.charAt(0).toUpperCase() + d.class.slice(1),
                    marker: {
                        color: classColors[d.class] || '#9CA3AF',
                        size: 12,
                        line: { color: 'white', width: 1 }
                    }
                };
            }
            traces[d.class].x.push(d.x);
            traces[d.class].y.push(d.y);
            traces[d.class].text.push(d.name);
            
            // Top 5 binding receptors
            const sorted = Object.entries(d.receptors)
                .sort((a, b) => a[1] - b[1])
                .slice(0, 5);
            const hoverText = d.name + '<br>' + 
                sorted.map(([r, e]) => `${r}: ${e.toFixed(0)} kJ/mol`).join('<br>');
            traces[d.class].hovertext.push(hoverText);
        });
        
        const layout = {
            paper_bgcolor: '#0a0a1a',
            plot_bgcolor: '#0a0a1a',
            font: { color: 'white' },
            xaxis: { 
                title: "''' + method + ''' Dimension 1",
                gridcolor: '#333',
                zerolinecolor: '#444'
            },
            yaxis: { 
                title: "''' + method + ''' Dimension 2",
                gridcolor: '#333',
                zerolinecolor: '#444'
            },
            hovermode: 'closest',
            legend: {
                bgcolor: 'rgba(26, 26, 46, 0.8)',
                bordercolor: '#333'
            },
            margin: { t: 60 }
        };
        
        Plotly.newPlot('plot', Object.values(traces), layout, {responsive: true});
    </script>
</body>
</html>'''
    
    output_path = '/data/users/zacharie/whc/consciousness_study/website/drug_binding_map.html'
    with open(output_path, 'w') as f:
        f.write(html_content)
    print(f"Saved interactive plot to {output_path}")

if __name__ == '__main__':
    main()
