#!/usr/bin/env python3
"""
Batch Processing Script for Receptor-Ligand Pipeline
=====================================================

Run the full Chai-1 → PoseBusters → FEP pipeline on priority pairs.

Usage:
    python run_batch_pipeline.py [--max-pairs N] [--test]
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

# Import pipeline
from pose_prediction_pipeline import run_full_pipeline

# Priority pairs for consciousness study
PRIORITY_PAIRS = [
    # Serotonin receptors + psychedelics
    ("5-HT2A", "DMT", "Serotonin"),
    ("5-HT2A", "LSD", "Serotonin"),
    ("5-HT2A", "Psilocin", "Serotonin"),
    ("5-HT2A", "5-MeO-DMT", "Serotonin"),
    ("5-HT2A", "Mescaline", "Serotonin"),
    ("5-HT2B", "LSD", "Serotonin"),
    ("5-HT2C", "DMT", "Serotonin"),
    ("5-HT1A", "DMT", "Serotonin"),
    
    # Dopamine system
    ("D2", "Dopamine", None),
    ("D2", "LSD", "Dopamine"),
    ("D1", "Dopamine", None),
    
    # Opioid system
    ("MOR", "Morphine", None),
    ("MOR", "Salvinorin_A", "Morphine"),
    ("KOR", "Salvinorin_A", None),
    
    # Glutamate system  
    ("NMDA_NR1", "Ketamine", None),
    ("NMDA_NR1", "PCP", "Ketamine"),
    ("NMDA_NR2B", "Ketamine", None),
    
    # Cannabinoid system
    ("CB1", "THC", "Anandamide"),
    ("CB2", "THC", None),
    
    # GABA system
    ("GABA_A", "Diazepam", "GABA"),
    
    # Nicotinic
    ("nAChR_a4b2", "Nicotine", "Acetylcholine"),
]


def main():
    parser = argparse.ArgumentParser(description="Batch pipeline processing")
    parser.add_argument("--max-pairs", type=int, default=None, 
                       help="Maximum number of pairs to process")
    parser.add_argument("--test", action="store_true",
                       help="Run in test mode (1 pair only)")
    args = parser.parse_args()
    
    # Determine pairs to process
    pairs = PRIORITY_PAIRS[:1] if args.test else PRIORITY_PAIRS
    if args.max_pairs:
        pairs = pairs[:args.max_pairs]
    
    print("=" * 70)
    print("CONSCIOUSNESS STUDY - BATCH PIPELINE")
    print(f"Started: {datetime.now().isoformat()}")
    print(f"Processing {len(pairs)} pairs")
    print("=" * 70)
    
    results = []
    for i, (receptor, ligand, reference) in enumerate(pairs, 1):
        print(f"\n[{i}/{len(pairs)}] Processing: {receptor} + {ligand}")
        
        try:
            result = run_full_pipeline(receptor, ligand, reference_ligand=reference)
            results.append({
                "receptor": receptor,
                "ligand": ligand,
                "reference": reference,
                "success": result.get("success", False),
                "chai1": result.get("chai1"),
                "posebusters": result.get("posebusters"),
                "fep": result.get("fep")
            })
        except Exception as e:
            print(f"    ✗ Error: {e}")
            results.append({
                "receptor": receptor,
                "ligand": ligand,
                "error": str(e)
            })
    
    # Save summary
    summary_file = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/batch_results.json")
    with open(summary_file, "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_pairs": len(pairs),
            "successful": sum(1 for r in results if r.get("success")),
            "results": results
        }, f, indent=2, default=str)
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    successful = sum(1 for r in results if r.get("success"))
    print(f"Total: {len(results)}")
    print(f"Successful: {successful}")
    print(f"Failed: {len(results) - successful}")
    print(f"\nResults saved to: {summary_file}")


if __name__ == "__main__":
    main()
