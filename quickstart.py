#!/usr/bin/env python
"""
Quick Start Demo
Demonstrates the consciousness research pipeline
"""

import asyncio
import json
from pathlib import Path


def demo_pharmacology():
    """Demonstrate pharmacology analysis"""
    print("\n" + "="*60)
    print("PHARMACOLOGY DEMO: Analyzing LSD receptor binding")
    print("="*60)
    
    from consciousness_study.src.pharmacology.receptor_mapping import PharmacologyAnalyzer
    
    analyzer = PharmacologyAnalyzer()
    
    # Get LSD pharmacology
    lsd = analyzer.get_receptor_profile("LSD")
    print(f"\n{lsd.name} Pharmacology:")
    print(f"  Half-life: {lsd.half_life_hours}h")
    print(f"  Primary targets: {lsd.primary_targets}")
    print("\n  Receptor bindings (sorted by affinity):")
    for b in sorted(lsd.bindings, key=lambda x: x.ki_nm or float('inf')):
        print(f"    {b.receptor}: Ki={b.ki_nm}nM ({b.activity.value})")
    
    # Predict effects
    print("\n  Predicted effects:")
    effects = analyzer.predict_effects_from_pharmacology("LSD")
    for e in effects[:8]:
        print(f"    {e['effect']} (via {e['receptor']}, intensity={e['predicted_intensity']:.1f})")
    
    # Compare substances
    print("\n  Similar compounds to LSD:")
    comparison = analyzer.compare_substances("LSD", "Psilocin")
    print(f"    Shared receptors with Psilocin: {comparison['shared_targets']}")


def demo_knowledge_graph():
    """Demonstrate knowledge graph queries"""
    print("\n" + "="*60)
    print("KNOWLEDGE GRAPH DEMO: Querying consciousness relationships")
    print("="*60)
    
    from consciousness_study.src.knowledge_graph.graph_builder import build_base_graph
    
    kg = build_base_graph()
    
    summary = kg.summary()
    print(f"\nGraph Statistics:")
    print(f"  Nodes: {summary['total_nodes']}")
    print(f"  Edges: {summary['total_edges']}")
    print(f"  Node types: {summary['node_types']}")
    
    # Query: What causes ego dissolution?
    print("\nQuery: Brain regions involved in ego dissolution")
    regions = kg.get_brain_regions_for_effect("ego_dissolution")
    for r in regions[:5]:
        print(f"  {r['region']}: {r['score']:.2f}")
    
    # Query: Effects of LSD
    print("\nQuery: Effects of LSD (via receptor mechanisms)")
    effects = kg.get_substance_effects("LSD")
    for e in effects[:8]:
        print(f"  {e['effect']} ({e['relationship']})")
    
    # Find similar substances
    print("\nQuery: Substances pharmacologically similar to LSD")
    similar = kg.find_similar_substances("LSD")
    for name, score in similar:
        print(f"  {name.split(':')[1]}: {score:.2f} similarity")


def demo_effect_taxonomy():
    """Demonstrate the effect taxonomy from PsychonautWiki"""
    print("\n" + "="*60)
    print("EFFECT TAXONOMY: Standardized phenomenological categories")
    print("="*60)
    
    from consciousness_study.src.nlp.effect_extractor import EFFECT_TAXONOMY
    
    for category, effects in EFFECT_TAXONOMY.items():
        print(f"\n{category.upper()}:")
        for effect in effects[:5]:
            print(f"  • {effect}")
        if len(effects) > 5:
            print(f"  ... and {len(effects)-5} more")


def demo_molecular_dynamics():
    """Demonstrate MD simulation setup"""
    print("\n" + "="*60)
    print("MOLECULAR DYNAMICS: Available receptor structures")
    print("="*60)
    
    from consciousness_study.src.molecular_dynamics.md_setup import (
        RECEPTOR_STRUCTURES, LIGAND_SMILES
    )
    
    print("\nReceptor structures available for MD simulations:")
    for name, struct in RECEPTOR_STRUCTURES.items():
        print(f"  {name}: PDB {struct.pdb_id}")
        print(f"    {struct.name}")
        print(f"    Binding site residues: {struct.binding_site_residues[:5]}...")
    
    print("\nLigand SMILES available:")
    for name, smiles in LIGAND_SMILES.items():
        print(f"  {name}: {smiles[:40]}...")


async def demo_full_pipeline():
    """Demonstrate the full research pipeline (requires API access)"""
    print("\n" + "="*60)
    print("FULL PIPELINE DEMO")
    print("="*60)
    print("\nThis demo shows the pipeline structure.")
    print("For actual data collection, you need:")
    print("  1. Network access to APIs")
    print("  2. Ollama running with llama3.2 model")
    print("  3. Rate-limited scraping (~20 reports)")
    
    print("\nPipeline stages:")
    print("  1. Data Collection:")
    print("     - PsychonautWiki API (substance info, effects)")
    print("     - TripSit API (factsheets, interactions)")
    print("     - Erowid scraping (experience reports)")
    print("  2. NLP Processing:")
    print("     - LLM-based effect extraction")
    print("     - Temporal phase identification")
    print("     - Standardization to taxonomy")
    print("  3. Analysis:")
    print("     - Effect frequency computation")
    print("     - Brain region mapping")
    print("     - Prediction vs observation comparison")
    print("  4. Knowledge Graph:")
    print("     - Multi-relational graph building")
    print("     - Mechanism pathway queries")
    print("  5. Molecular Dynamics:")
    print("     - Receptor-ligand binding simulations")
    print("     - Binding site analysis")


def main():
    """Run all demos"""
    print("\n" + "#"*60)
    print("CONSCIOUSNESS RESEARCH TOOLKIT - QUICK START DEMO")
    print("#"*60)
    
    # Run demos that don't require external resources
    demo_pharmacology()
    demo_knowledge_graph()
    demo_effect_taxonomy()
    demo_molecular_dynamics()
    
    # Show pipeline structure
    asyncio.run(demo_full_pipeline())
    
    print("\n" + "="*60)
    print("NEXT STEPS")
    print("="*60)
    print("""
1. Install dependencies:
   pip install -r requirements.txt
   
2. Install Ollama and pull llama3.2:
   ollama pull llama3.2
   
3. Run the full analysis:
   python -m consciousness_study.src.analysis.integration_pipeline
   
4. For molecular dynamics:
   - Install OpenMM: conda install -c conda-forge openmm
   - Install GROMACS: sudo apt install gromacs
   - Download receptor structures from RCSB PDB
   
5. Explore the knowledge graph:
   - Export to Neo4j for visualization
   - Run custom queries on effect-receptor-region relationships
""")


if __name__ == "__main__":
    main()
