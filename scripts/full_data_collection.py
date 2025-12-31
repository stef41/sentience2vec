#!/usr/bin/env python3
"""
Full Data Collection Script
Scrapes all three data sources: PsychonautWiki, TripSit, and Erowid
"""

import json
import time
import os
import sys
from pathlib import Path
from datetime import datetime

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_collection.tripsit_client import TripSitAPI
from src.data_collection.psychonautwiki_client import PsychonautWikiAPI
from src.data_collection.erowid_scraper import ErowidScraper


def collect_tripsit_data(output_dir: Path):
    """Collect all data from TripSit API"""
    print("\n" + "="*60)
    print("TRIPSIT DATA COLLECTION")
    print("="*60)
    
    api = TripSitAPI()
    
    # Get all drug names
    print("\nFetching all drug names...")
    drug_names = api.get_all_drug_names()
    print(f"Found {len(drug_names)} drugs")
    
    # Get all categories
    print("\nFetching categories...")
    categories = api.get_all_categories()
    print(f"Categories: {categories}")
    
    # Get all drug data (more efficient than individual calls)
    print("\nFetching all drug data...")
    all_drugs = api.get_all_drugs()
    print(f"Retrieved {len(all_drugs)} drug records")
    
    # Save data
    tripsit_dir = output_dir / "tripsit"
    tripsit_dir.mkdir(parents=True, exist_ok=True)
    
    with open(tripsit_dir / "all_drugs.json", "w") as f:
        json.dump(all_drugs, f, indent=2)
    
    with open(tripsit_dir / "drug_names.json", "w") as f:
        json.dump(drug_names, f, indent=2)
    
    with open(tripsit_dir / "categories.json", "w") as f:
        json.dump(categories, f, indent=2)
    
    # Get alias mapping
    print("\nFetching alias mappings...")
    aliases = api.get_all_drug_name_aliases()
    with open(tripsit_dir / "aliases.json", "w") as f:
        json.dump(aliases, f, indent=2)
    print(f"Retrieved {len(aliases)} alias mappings")
    
    # Extract combo/interaction data
    combo_data = {}
    for drug in all_drugs:
        name = drug.get("name") or drug.get("pretty_name")
        if name and drug.get("combos"):
            combo_data[name] = drug["combos"]
    
    with open(tripsit_dir / "drug_combinations.json", "w") as f:
        json.dump(combo_data, f, indent=2)
    print(f"Extracted combination data for {len(combo_data)} substances")
    
    return {
        "drugs": len(all_drugs),
        "categories": len(categories),
        "aliases": len(aliases),
        "combinations": len(combo_data)
    }


def collect_psychonautwiki_data(output_dir: Path):
    """Collect data from PsychonautWiki GraphQL API"""
    print("\n" + "="*60)
    print("PSYCHONAUTWIKI DATA COLLECTION")
    print("="*60)
    
    api = PsychonautWikiAPI()
    
    # Target substances (common psychoactives for consciousness research)
    target_substances = [
        # Classic psychedelics
        "LSD", "Psilocybin", "Psilocin", "DMT", "5-MeO-DMT", "Mescaline",
        "2C-B", "2C-E", "2C-I", "2C-C", "DOC", "DOM", "DOI",
        # Tryptamines
        "4-AcO-DMT", "4-HO-MET", "4-HO-MiPT", "DPT", "DiPT",
        # Lysergamides
        "1P-LSD", "AL-LAD", "ETH-LAD", "LSA", "ALD-52",
        # Dissociatives
        "Ketamine", "DXM", "PCP", "MXE", "3-MeO-PCP", "3-HO-PCP",
        "Nitrous", "Methoxetamine",
        # Empathogens
        "MDMA", "MDA", "6-APB", "5-MAPB", "Methylone",
        # Deliriants
        "DPH", "Scopolamine", "Datura",
        # Atypical
        "Salvinorin A", "Ibogaine", "Cannabis", "Amanita muscaria",
        # Common recreational
        "Cocaine", "Amphetamine", "Methamphetamine", 
        "GHB", "Alcohol", "Caffeine", "Nicotine",
        # Opioids (for comparison)
        "Morphine", "Heroin", "Fentanyl", "Kratom",
        # Benzodiazepines
        "Alprazolam", "Diazepam", "Clonazepam",
    ]
    
    psychonaut_dir = output_dir / "psychonautwiki"
    psychonaut_dir.mkdir(parents=True, exist_ok=True)
    
    substances_data = []
    effects_collected = set()
    failed = []
    
    print(f"\nFetching {len(target_substances)} substances...")
    
    for i, name in enumerate(target_substances):
        try:
            print(f"  [{i+1}/{len(target_substances)}] {name}...", end=" ")
            substance = api.get_substance(name)
            
            if substance:
                substances_data.append(substance)
                # Collect unique effects
                effects = substance.get("effects", [])
                for effect in effects:
                    if isinstance(effect, dict):
                        effects_collected.add(effect.get("name", ""))
                    elif isinstance(effect, str):
                        effects_collected.add(effect)
                print(f"✓ ({len(effects)} effects)")
            else:
                failed.append(name)
                print("✗ (not found)")
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            failed.append(name)
            print(f"✗ ({e})")
    
    # Save substance data
    with open(psychonaut_dir / "substances.json", "w") as f:
        json.dump(substances_data, f, indent=2, default=str)
    
    # Save effects list
    effects_list = sorted(list(effects_collected))
    with open(psychonaut_dir / "effects.json", "w") as f:
        json.dump(effects_list, f, indent=2)
    
    # Save failed queries for debugging
    with open(psychonaut_dir / "failed_queries.json", "w") as f:
        json.dump(failed, f, indent=2)
    
    # Try to get all effects from the API
    print("\nFetching complete effects taxonomy...")
    try:
        all_effects = api.get_all_effects()
        with open(psychonaut_dir / "all_effects.json", "w") as f:
            json.dump(all_effects, f, indent=2, default=str)
        print(f"Retrieved {len(all_effects)} effects from taxonomy")
    except Exception as e:
        print(f"Could not fetch effects taxonomy: {e}")
    
    print(f"\nResults: {len(substances_data)} substances, {len(effects_collected)} unique effects")
    if failed:
        print(f"Failed: {failed}")
    
    return {
        "substances": len(substances_data),
        "effects": len(effects_collected),
        "failed": len(failed)
    }


def collect_erowid_data(output_dir: Path, max_reports_per_substance: int = 20):
    """Collect experience reports from Erowid"""
    print("\n" + "="*60)
    print("EROWID DATA COLLECTION")
    print("="*60)
    
    scraper = ErowidScraper(cache_dir=str(output_dir / "erowid" / "cache"))
    
    erowid_dir = output_dir / "erowid"
    erowid_dir.mkdir(parents=True, exist_ok=True)
    
    # First get the substance list from Erowid
    print("\nFetching substance list from Erowid...")
    substance_list = scraper.get_substance_list()
    print(f"Found {len(substance_list)} substances with experience vaults")
    
    # Save full substance list
    with open(erowid_dir / "substance_list.json", "w") as f:
        json.dump(substance_list, f, indent=2)
    
    # Priority substances for consciousness research (map to Erowid names/URLs)
    priority_keywords = [
        "lsd", "mushroom", "psilocybin", "dmt", "ayahuasca", "mescaline", "peyote",
        "2c-b", "2c-e", "mdma", "mda", "ketamine", "dxm", "pcp", "nitrous",
        "salvia", "ibogaine", "cannabis", "marijuana", "datura", "diphenhydramine"
    ]
    
    # Filter substance list to priority substances
    target_substances = []
    for sub in substance_list:
        name_lower = sub["name"].lower()
        if any(kw in name_lower for kw in priority_keywords):
            target_substances.append(sub)
    
    print(f"Targeting {len(target_substances)} priority substances")
    
    all_reports = []
    substance_counts = {}
    
    print(f"\nScraping experience reports...")
    print(f"(Max {max_reports_per_substance} reports per substance)")
    print("Note: This respects rate limits (2s delay) and may take time...\n")
    
    for sub in target_substances[:15]:  # Limit to first 15 for initial collection
        name = sub["name"]
        url = sub["url"]
        
        try:
            print(f"  {name}...", end=" ", flush=True)
            
            reports = scraper.scrape_substance(
                substance_name=name,
                substance_url=url,
                max_reports=max_reports_per_substance,
                output_dir=str(erowid_dir)
            )
            
            substance_counts[name] = len(reports)
            all_reports.extend(reports)
            
            print(f"✓ ({len(reports)} reports)")
            
        except Exception as e:
            print(f"✗ ({e})")
            substance_counts[name] = 0
    
    # Save combined data
    print("\nSaving combined dataset...")
    from dataclasses import asdict
    with open(erowid_dir / "all_reports.json", "w") as f:
        json.dump([asdict(r) for r in all_reports], f, indent=2)
    
    with open(erowid_dir / "collection_stats.json", "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "substance_counts": substance_counts,
            "total_reports": len(all_reports)
        }, f, indent=2)
    
    print(f"\nTotal reports collected: {len(all_reports)}")
    
    return {
        "reports": len(all_reports),
        "substances": len([s for s, c in substance_counts.items() if c > 0])
    }


def main():
    """Run full data collection pipeline"""
    print("="*60)
    print("CONSCIOUSNESS STUDY - FULL DATA COLLECTION")
    print(f"Started: {datetime.now().isoformat()}")
    print("="*60)
    
    output_dir = Path("data/raw")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # 1. TripSit (fast, bulk API)
    try:
        results["tripsit"] = collect_tripsit_data(output_dir)
    except Exception as e:
        print(f"\nTripSit collection failed: {e}")
        results["tripsit"] = {"error": str(e)}
    
    # 2. PsychonautWiki (GraphQL, individual queries)
    try:
        results["psychonautwiki"] = collect_psychonautwiki_data(output_dir)
    except Exception as e:
        print(f"\nPsychonautWiki collection failed: {e}")
        results["psychonautwiki"] = {"error": str(e)}
    
    # 3. Erowid (scraping, slower)
    # Only collect a sample due to rate limiting
    try:
        results["erowid"] = collect_erowid_data(output_dir, max_reports_per_substance=20)
    except Exception as e:
        print(f"\nErowid collection failed: {e}")
        results["erowid"] = {"error": str(e)}
    
    # Summary
    print("\n" + "="*60)
    print("COLLECTION SUMMARY")
    print("="*60)
    
    for source, data in results.items():
        print(f"\n{source.upper()}:")
        for key, value in data.items():
            print(f"  {key}: {value}")
    
    # Save summary
    with open(output_dir / "collection_summary.json", "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "results": results
        }, f, indent=2)
    
    print(f"\n\nData saved to: {output_dir.absolute()}")
    print(f"Completed: {datetime.now().isoformat()}")


if __name__ == "__main__":
    main()
