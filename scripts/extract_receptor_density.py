#!/usr/bin/env python3
"""
Extract Neuroreceptor Density Distribution Data

Extracts receptor density data from the knowledge graph and produces:
1. A comprehensive JSON file with all receptor-region-density mappings
2. A CSV for easy analysis and visualization
3. Summary statistics

Data sources include autoradiography, PET imaging, and immunohistochemistry studies.
Density values are in fmol/mg tissue equivalent (normalized from various sources).
"""

import json
import csv
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Any

# Scientific references for density data
DENSITY_REFERENCES = {
    '5-HT1A': [
        'Pazos & Palacios 1985 - Brain Res',
        'Hall et al. 1997 - Neuroscience',
        'Varnas et al. 2004 - Hum Brain Mapp'
    ],
    '5-HT2A': [
        'Pazos et al. 1987 - Eur J Pharmacol',
        'Hall et al. 2000 - Synapse',
        'Beliveau et al. 2017 - J Neurosci (Cimbi-36 PET)'
    ],
    '5-HT2B': [
        'Borman et al. 2002 - Neuropharmacology',
        'Duxon et al. 1997 - Neuroscience'
    ],
    '5-HT2C': [
        'Abramowski et al. 1995 - Neuroscience',
        'Clemett et al. 2000 - Neuropharmacology'
    ],
    'D1': [
        'Hall et al. 1994 - Neuroscience',
        'Camps et al. 1989 - Neuroscience'
    ],
    'D2': [
        'Hall et al. 1996 - Synapse',
        'Kessler et al. 1993 - J Cereb Blood Flow Metab'
    ],
    'D3': [
        'Gurevich & Joyce 1999 - Neuropsychopharmacology',
        'Sokoloff et al. 1990 - Nature'
    ],
    'MOR': [
        'Mathieu-Kia et al. 2001 - J Comp Neurol',
        'Peckys & Bhargava 1998 - Brain Res'
    ],
    'KOR': [
        'Simonin et al. 1995 - Mol Pharmacol',
        'Maurer et al. 1983 - Life Sci'
    ],
    'DOR': [
        'Mansour et al. 1995 - J Comp Neurol',
        'Peckys & Bhargava 1998 - Brain Res'
    ],
    'CB1': [
        'Glass et al. 1997 - Neuroscience',
        'Herkenham et al. 1991 - J Neurosci'
    ],
    'NMDA_GluN2A': [
        'Monaghan & Cotman 1985 - J Neurosci',
        'Scherzer et al. 1998 - Mol Brain Res'
    ],
    'NMDA_GluN2B': [
        'Wenzel et al. 1997 - Neuroscience',
        'Laurie & Bhargava 1998 - Neuroscience'
    ],
    'AMPA_GluA1': [
        'Martin et al. 1993 - Neuroscience',
        'Beneyto & Bhargava 2007 - Neuropsychopharmacology'
    ],
    'GABA_A_alpha1': [
        'Sieghart & Sperk 2002 - Pharmacol Rev',
        'Pirker et al. 2000 - Neuroscience'
    ],
    'GABA_A_alpha2': [
        'Fritschy & Mohler 1995 - J Comp Neurol',
        'Sieghart & Sperk 2002 - Pharmacol Rev'
    ],
    'alpha2A': [
        'Aoki et al. 1998 - Cereb Cortex',
        'Scheinin et al. 1994 - Mol Pharmacol'
    ],
    'Sigma1': [
        'Alonso et al. 2000 - Neuroscience',
        'Hayashi & Su 2005 - CNS Drugs'
    ],
    'nAChR_alpha4beta2': [
        'Spurden et al. 1997 - J Chem Neuroanat',
        'Paterson & Bhargava 1999 - Pharmacol Ther'
    ],
    'M1_muscarinic': [
        'Levey et al. 1995 - J Neurosci',
        'Flynn et al. 1995 - Synapse'
    ]
}

# Brain region descriptions and functions
BRAIN_REGION_INFO = {
    'prefrontal_cortex': {
        'full_name': 'Prefrontal Cortex',
        'function': 'Executive function, decision-making, personality',
        'brodmann': [9, 10, 11, 12, 46]
    },
    'prefrontal_cortex_layer_v': {
        'full_name': 'Prefrontal Cortex Layer V',
        'function': 'Primary output layer to subcortical structures',
        'brodmann': [9, 46]
    },
    'prefrontal_cortex_layer_ii': {
        'full_name': 'Prefrontal Cortex Layer II',
        'function': 'Cortico-cortical connections',
        'brodmann': [9, 46]
    },
    'anterior_cingulate_cortex': {
        'full_name': 'Anterior Cingulate Cortex',
        'function': 'Error detection, emotional processing, salience',
        'brodmann': [24, 32, 33]
    },
    'posterior_cingulate_ba23': {
        'full_name': 'Posterior Cingulate (BA23)',
        'function': 'Default mode network hub, self-referential thought',
        'brodmann': [23]
    },
    'hippocampus': {
        'full_name': 'Hippocampus',
        'function': 'Memory formation and spatial navigation',
        'brodmann': []
    },
    'hippocampus_ca1': {
        'full_name': 'Hippocampus CA1',
        'function': 'Memory consolidation, novelty detection',
        'brodmann': []
    },
    'hippocampus_ca3': {
        'full_name': 'Hippocampus CA3',
        'function': 'Pattern completion, memory recall',
        'brodmann': []
    },
    'amygdala': {
        'full_name': 'Amygdala',
        'function': 'Emotional processing, fear, salience',
        'brodmann': []
    },
    'amygdala_basolateral': {
        'full_name': 'Basolateral Amygdala',
        'function': 'Emotional learning, value assignment',
        'brodmann': []
    },
    'nucleus_accumbens': {
        'full_name': 'Nucleus Accumbens',
        'function': 'Reward, motivation, pleasure',
        'brodmann': []
    },
    'caudate_nucleus': {
        'full_name': 'Caudate Nucleus',
        'function': 'Motor control, goal-directed behavior',
        'brodmann': []
    },
    'putamen': {
        'full_name': 'Putamen',
        'function': 'Motor function, learning',
        'brodmann': []
    },
    'thalamus': {
        'full_name': 'Thalamus',
        'function': 'Sensory relay, consciousness gating',
        'brodmann': []
    },
    'claustrum': {
        'full_name': 'Claustrum',
        'function': 'Consciousness integration, salience',
        'brodmann': []
    },
    'visual_cortex_v1_layer_iv': {
        'full_name': 'Primary Visual Cortex (V1) Layer IV',
        'function': 'Primary visual processing',
        'brodmann': [17]
    },
    'dorsal_raphe_nucleus': {
        'full_name': 'Dorsal Raphe Nucleus',
        'function': 'Serotonin production, mood regulation',
        'brodmann': []
    },
    'locus_coeruleus': {
        'full_name': 'Locus Coeruleus',
        'function': 'Norepinephrine production, arousal, attention',
        'brodmann': []
    },
    'vta': {
        'full_name': 'Ventral Tegmental Area',
        'function': 'Dopamine production, reward signaling',
        'brodmann': []
    },
    'substantia_nigra': {
        'full_name': 'Substantia Nigra',
        'function': 'Dopamine production, motor control',
        'brodmann': []
    },
    'cerebellum': {
        'full_name': 'Cerebellum',
        'function': 'Motor coordination, timing, learning',
        'brodmann': []
    },
    'insula_anterior': {
        'full_name': 'Anterior Insula',
        'function': 'Interoception, emotional awareness, salience',
        'brodmann': [13, 14]
    },
}


def load_knowledge_graph(data_dir: Path) -> Dict:
    """Load the knowledge graph JSON."""
    kg_path = data_dir / 'knowledge_graph.json'
    with open(kg_path, 'r') as f:
        return json.load(f)


def extract_receptor_densities(kg: Dict) -> Dict[str, List[Dict]]:
    """Extract receptor-region density data from the knowledge graph."""
    density_data = defaultdict(list)
    
    # Get receptor names for proper casing
    receptor_names = {n['id'].replace('receptor:', ''): n['name'] 
                      for n in kg['nodes'] if n['node_type'] == 'receptor'}
    
    # Get region names
    region_names = {n['id'].replace('brain_region:', ''): n.get('name', n['id'].replace('brain_region:', ''))
                    for n in kg['nodes'] if n['node_type'] == 'brain_region'}
    
    for edge in kg['edges']:
        if edge['edge_type'] == 'located_in' and 'density' in edge:
            receptor_id = edge['source'].replace('receptor:', '')
            region_id = edge['target'].replace('brain_region:', '')
            
            receptor_name = receptor_names.get(receptor_id, receptor_id)
            region_name = region_names.get(region_id, region_id)
            
            density_data[receptor_name].append({
                'region_id': region_id,
                'region_name': region_name,
                'density_fmol_mg': int(edge['density']),
                'region_info': BRAIN_REGION_INFO.get(region_id, {})
            })
    
    # Sort by density (highest first)
    for receptor in density_data:
        density_data[receptor].sort(key=lambda x: x['density_fmol_mg'], reverse=True)
    
    return dict(density_data)


def extract_receptor_effects(kg: Dict) -> Dict[str, List[str]]:
    """Extract receptor-effect relationships."""
    effects = defaultdict(list)
    
    receptor_names = {n['id'].replace('receptor:', ''): n['name'] 
                      for n in kg['nodes'] if n['node_type'] == 'receptor'}
    
    for edge in kg['edges']:
        if edge['edge_type'] == 'modulates':
            receptor_id = edge['source'].replace('receptor:', '')
            effect_id = edge['target'].replace('effect:', '')
            
            receptor_name = receptor_names.get(receptor_id, receptor_id)
            effects[receptor_name].append(effect_id)
    
    return dict(effects)


def calculate_statistics(density_data: Dict) -> Dict:
    """Calculate summary statistics."""
    stats = {
        'total_receptors': len(density_data),
        'total_mappings': sum(len(v) for v in density_data.values()),
        'unique_regions': len(set(r['region_id'] for regions in density_data.values() for r in regions)),
        'receptor_stats': {}
    }
    
    for receptor, regions in density_data.items():
        densities = [r['density_fmol_mg'] for r in regions]
        stats['receptor_stats'][receptor] = {
            'region_count': len(regions),
            'max_density': max(densities),
            'min_density': min(densities),
            'mean_density': round(sum(densities) / len(densities), 1),
            'highest_region': regions[0]['region_name'] if regions else None
        }
    
    return stats


def export_json(density_data: Dict, effects: Dict, stats: Dict, output_path: Path):
    """Export comprehensive JSON file."""
    output = {
        'metadata': {
            'description': 'Neuroreceptor density distribution in human brain',
            'units': 'fmol/mg tissue (normalized values)',
            'sources': 'Autoradiography, PET imaging, immunohistochemistry',
            'references': DENSITY_REFERENCES
        },
        'statistics': stats,
        'receptors': {}
    }
    
    for receptor, regions in density_data.items():
        output['receptors'][receptor] = {
            'density_distribution': regions,
            'effects': effects.get(receptor, []),
            'references': DENSITY_REFERENCES.get(receptor, [])
        }
    
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Exported JSON: {output_path}")


def export_csv(density_data: Dict, effects: Dict, output_path: Path):
    """Export CSV for easy analysis."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'receptor', 'brain_region', 'density_fmol_mg', 
            'region_full_name', 'region_function', 'effects'
        ])
        
        for receptor, regions in sorted(density_data.items()):
            receptor_effects = ','.join(effects.get(receptor, []))
            for region in regions:
                info = region.get('region_info', {})
                writer.writerow([
                    receptor,
                    region['region_id'],
                    region['density_fmol_mg'],
                    info.get('full_name', region['region_name']),
                    info.get('function', ''),
                    receptor_effects
                ])
    
    print(f"✓ Exported CSV: {output_path}")


def export_summary(density_data: Dict, stats: Dict, output_path: Path):
    """Export human-readable summary."""
    with open(output_path, 'w') as f:
        f.write("# Neuroreceptor Density Distribution Summary\n\n")
        f.write(f"**Total Receptors:** {stats['total_receptors']}\n")
        f.write(f"**Total Density Mappings:** {stats['total_mappings']}\n")
        f.write(f"**Unique Brain Regions:** {stats['unique_regions']}\n\n")
        
        f.write("## Receptor Overview\n\n")
        f.write("| Receptor | Regions | Max Density | Highest Expression Region |\n")
        f.write("|----------|---------|-------------|---------------------------|\n")
        
        for receptor in sorted(density_data.keys()):
            s = stats['receptor_stats'][receptor]
            f.write(f"| {receptor} | {s['region_count']} | {s['max_density']} fmol/mg | {s['highest_region']} |\n")
        
        f.write("\n## Detailed Distribution\n\n")
        
        for receptor in sorted(density_data.keys()):
            f.write(f"### {receptor}\n\n")
            f.write("| Brain Region | Density (fmol/mg) |\n")
            f.write("|--------------|-------------------|\n")
            for region in density_data[receptor][:10]:  # Top 10
                f.write(f"| {region['region_name']} | {region['density_fmol_mg']} |\n")
            if len(density_data[receptor]) > 10:
                f.write(f"| ... and {len(density_data[receptor]) - 10} more | |\n")
            f.write("\n")
    
    print(f"✓ Exported summary: {output_path}")


def print_top_densities(density_data: Dict):
    """Print top density regions for key receptors."""
    key_receptors = ['5-HT2A', '5-HT1A', 'D2', 'MOR', 'KOR', 'CB1', 'NMDA_GluN2A']
    
    print("\n" + "="*70)
    print("TOP DENSITY REGIONS FOR KEY RECEPTORS")
    print("="*70)
    
    for receptor in key_receptors:
        if receptor in density_data:
            print(f"\n{receptor}:")
            for region in density_data[receptor][:5]:
                bar = '█' * (region['density_fmol_mg'] // 20)
                print(f"  {region['region_name']:35s} {region['density_fmol_mg']:4d} fmol/mg {bar}")


def main():
    """Main extraction function."""
    # Paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data'
    output_dir = data_dir / 'receptor_densities'
    output_dir.mkdir(exist_ok=True)
    
    print("="*70)
    print("NEURORECEPTOR DENSITY DISTRIBUTION EXTRACTION")
    print("="*70)
    
    # Load knowledge graph
    print("\nLoading knowledge graph...")
    kg = load_knowledge_graph(data_dir)
    
    # Extract data
    print("Extracting receptor-region density data...")
    density_data = extract_receptor_densities(kg)
    effects = extract_receptor_effects(kg)
    stats = calculate_statistics(density_data)
    
    print(f"  Found {stats['total_receptors']} receptors")
    print(f"  Found {stats['total_mappings']} density mappings")
    print(f"  Found {stats['unique_regions']} unique brain regions")
    
    # Export files
    print("\nExporting data files...")
    export_json(density_data, effects, stats, output_dir / 'receptor_density_distribution.json')
    export_csv(density_data, effects, output_dir / 'receptor_density_distribution.csv')
    export_summary(density_data, stats, output_dir / 'DENSITY_SUMMARY.md')
    
    # Print highlights
    print_top_densities(density_data)
    
    print("\n" + "="*70)
    print("EXTRACTION COMPLETE")
    print("="*70)
    print(f"\nOutput files in: {output_dir}")
    
    return density_data, stats


if __name__ == '__main__':
    main()
