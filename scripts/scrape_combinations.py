#!/usr/bin/env python3
"""
Comprehensive Drug Combination Scraper
Scrapes interaction data from multiple sources:
1. Mixtures.info API (35 substances, 195 interactions)
2. TripSit Combo Chart (from wiki and existing drug data)
3. PsychonautWiki interactions
4. Polydrug extraction from Erowid reports
"""

import json
import os
import re
import time
import requests
from pathlib import Path
from collections import defaultdict
from datetime import datetime

BASE_DIR = Path("/data/users/zacharie/whc/consciousness_study")
DATA_DIR = BASE_DIR / "data" / "raw"
COMBO_DIR = DATA_DIR / "combinations"
COMBO_DIR.mkdir(parents=True, exist_ok=True)

# Risk level standardization
RISK_LEVELS = {
    # Mixtures.info levels (text)
    "NEUTRAL": "low_risk_no_synergy",
    "POSITIVE": "low_risk_synergy",
    "CAUTION": "caution",
    "NEGATIVE": "unsafe",
    "DANGEROUS": "dangerous",
    "UNKNOWN": "unknown",
    # Mixtures.info levels (numeric 0-4 scale)
    0: "low_risk_synergy",       # Positive synergy
    1: "low_risk_no_synergy",    # Neutral
    2: "caution",                # Caution
    3: "unsafe",                 # Risky
    4: "dangerous",              # Dangerous
    # TripSit levels
    "Low Risk & Synergy": "low_risk_synergy",
    "Low Risk & No Synergy": "low_risk_no_synergy", 
    "Low Risk & Decrease": "low_risk_decrease",
    "Caution": "caution",
    "Unsafe": "unsafe",
    "Dangerous": "dangerous",
    "Serotonin Syndrome": "dangerous",
}


def scrape_mixtures_api():
    """Scrape all data from Mixtures.info API"""
    print("\n" + "="*60)
    print("SCRAPING MIXTURES.INFO API")
    print("="*60)
    
    base_url = "https://mixtures.info/en/api/v1"
    session = requests.Session()
    session.headers.update({"User-Agent": "ConsciousnessResearch/1.0"})
    
    data = {
        "substances": [],
        "interactions": [],
        "aliases": {},
        "source": "mixtures.info",
        "scraped_at": datetime.now().isoformat()
    }
    
    # Get all aliases for local search
    print("\n[1/4] Fetching aliases...")
    try:
        resp = session.get(f"{base_url}/aliases/", timeout=30)
        if resp.status_code == 200:
            data["aliases"] = resp.json()
            print(f"  ✓ Got {len(data['aliases'])} aliases")
    except Exception as e:
        print(f"  ✗ Error: {e}")
    
    time.sleep(1)
    
    # Get all substances
    print("\n[2/4] Fetching substance list...")
    try:
        resp = session.get(f"{base_url}/substances/", timeout=30)
        if resp.status_code == 200:
            substances_dict = resp.json()
            # Convert dict to list with slug included
            substances = []
            for slug, sub_data in substances_dict.items():
                sub_data["slug"] = slug
                substances.append(sub_data)
            data["substances"] = substances
            print(f"  ✓ Got {len(substances)} substances")
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return data
    
    time.sleep(1)
    
    # Get detailed info for each substance
    print("\n[3/4] Fetching detailed substance info...")
    for i, sub in enumerate(data["substances"]):
        slug = sub.get("slug")
        if slug:
            try:
                resp = session.get(f"{base_url}/substance/{slug}/", timeout=30)
                if resp.status_code == 200:
                    detail = resp.json()
                    sub.update(detail)
                    
                    # Extract interactions from the response
                    interactions = detail.get("interactions", {})
                    if isinstance(interactions, dict):
                        for combo_key, interaction_data in interactions.items():
                            interactants = interaction_data.get("interactants", [])
                            # Get the other substance
                            other = [x for x in interactants if x.lower() != slug.lower()]
                            if other:
                                data["interactions"].append({
                                    "substance_a": slug,
                                    "substance_b": other[0].lower(),
                                    "risk": interaction_data.get("risk"),
                                    "synergy": interaction_data.get("synergy"),
                                    "risk_description": interaction_data.get("risk_description", ""),
                                    "effects_description": interaction_data.get("effects_description", ""),
                                    "url": interaction_data.get("site_url"),
                                    "source": "mixtures.info"
                                })
                
                print(f"  [{i+1}/{len(data['substances'])}] {slug}", end="\r")
                time.sleep(0.5)  # Rate limiting
            except Exception as e:
                print(f"  Error on {slug}: {e}")
    
    print(f"\n  ✓ Got details for {len(data['substances'])} substances")
    print(f"  ✓ Got {len(data['interactions'])} interaction records")
    
    # Deduplicate interactions (A+B is same as B+A)
    seen = set()
    unique_interactions = []
    for interaction in data["interactions"]:
        pair = tuple(sorted([interaction["substance_a"], interaction["substance_b"]]))
        if pair not in seen:
            seen.add(pair)
            unique_interactions.append(interaction)
    
    data["interactions"] = unique_interactions
    print(f"  ✓ {len(unique_interactions)} unique interactions after dedup")
    
    # Save
    output_path = COMBO_DIR / "mixtures_info.json"
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\n✓ Saved to {output_path}")
    
    return data


def scrape_tripsit_combos():
    """Extract combo data from TripSit drug database (already scraped)"""
    print("\n" + "="*60)
    print("EXTRACTING TRIPSIT COMBINATION DATA")
    print("="*60)
    
    tripsit_path = DATA_DIR / "tripsit" / "all_drugs_full.json"
    if not tripsit_path.exists():
        print(f"  ✗ TripSit data not found at {tripsit_path}")
        return {}
    
    with open(tripsit_path) as f:
        tripsit_data = json.load(f)
    
    data = {
        "combinations": [],
        "substances_with_combos": 0,
        "source": "tripsit",
        "scraped_at": datetime.now().isoformat()
    }
    
    # Also parse the wiki page data we fetched
    wiki_combinations = parse_tripsit_wiki_combinations()
    
    # The data is structured as a dict with drug names as keys
    # Each drug may have a 'combos' field
    for drug_name, drug_data in tripsit_data.items():
        if not isinstance(drug_data, dict):
            continue
            
        combos = drug_data.get("combos", {})
        if combos:
            data["substances_with_combos"] += 1
            
            for other_drug, combo_info in combos.items():
                status = combo_info.get("status", "unknown")
                note = combo_info.get("note", "")
                
                data["combinations"].append({
                    "substance_a": drug_name.lower(),
                    "substance_b": other_drug.lower(),
                    "risk": RISK_LEVELS.get(status, status.lower().replace(" ", "_").replace("&", "").replace("  ", "_")),
                    "status_original": status,
                    "note": note,
                    "source": "tripsit_api"
                })
    
    # Merge wiki combinations
    data["combinations"].extend(wiki_combinations)
    
    print(f"  ✓ Found {data['substances_with_combos']} substances with combo data")
    print(f"  ✓ Extracted {len(data['combinations'])} combination records")
    
    # Save
    output_path = COMBO_DIR / "tripsit_combos.json"
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  ✓ Saved to {output_path}")
    
    return data


def parse_tripsit_wiki_combinations():
    """Parse the TripSit wiki page for additional combination data"""
    # This would parse the wiki content we fetched
    # The wiki has standardized format: "### substance_a & substance_b"
    # with Status: and Note: fields
    
    combinations = []
    
    # Wiki combination data (parsed from fetch results)
    wiki_combos = [
        ("cannabis", "lsd", "Caution", "Cannabis has an unexpectedly strong and somewhat unpredictable synergy with psychedelics."),
        ("amphetamines", "lsd", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("cocaine", "lsd", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("tramadol", "lsd", "Unsafe", "Tramadol is well known to lower seizure threshold and psychedelics also cause occasional seizures."),
        ("cannabis", "mushrooms", "Caution", "Cannabis has an unexpectedly strong and somewhat unpredictable synergy with psychedelics."),
        ("amphetamines", "mushrooms", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("cocaine", "mushrooms", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("tramadol", "mushrooms", "Unsafe", "Tramadol is well known to lower seizure threshold and psychedelics also cause occasional seizures."),
        ("cannabis", "dmt", "Caution", "Cannabis has an unexpectedly strong and somewhat unpredictable synergy with psychedelics."),
        ("amphetamines", "dmt", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("cocaine", "dmt", "Caution", "Stimulants increase anxiety levels and the risk of thought loops which can lead to negative experiences"),
        ("tramadol", "dmt", "Unsafe", "Tramadol is well known to lower seizure threshold and psychedelics also cause occasional seizures."),
        ("alcohol", "ketamine", "Dangerous", "Both substances cause ataxia and bring a very high risk of vomiting and unconsciousness."),
        ("ghb", "ketamine", "Dangerous", "Both substances cause ataxia and bring a risk of vomiting and unconsciousness."),
        ("opioids", "ketamine", "Dangerous", "Both substances bring a risk of vomiting and unconsciousness."),
        ("alcohol", "mxe", "Dangerous", "There is a high risk of memory loss, vomiting and severe ataxia from this combination."),
        ("ghb", "mxe", "Dangerous", "Both substances cause ataxia and bring a risk of vomiting and unconsciousness."),
        ("opioids", "mxe", "Dangerous", "This combination can potentiate the effects of the opioid"),
        ("maois", "dxm", "Dangerous", "High risk of serotonin syndrome"),
        ("ssris", "dxm", "Dangerous", "High risk of serotonin syndrome."),
        ("alcohol", "dxm", "Dangerous", "Both substances potentiate the ataxia and sedation, CNS depression can lead to difficulty breathing."),
        ("ghb", "dxm", "Dangerous", "Both substances cause ataxia and bring a risk of vomiting and unconsciousness."),
        ("opioids", "dxm", "Dangerous", "CNS depression, difficult breathing, heart issues, hepatoxic."),
        ("tramadol", "amphetamines", "Dangerous", "Tramadol and stimulants both increase the risk of seizures."),
        ("tramadol", "mdma", "Dangerous", "Tramadol and stimulants both increase the risk of seizures."),
        ("tramadol", "cocaine", "Dangerous", "Tramadol and stimulants both increase the risk of seizures."),
        ("maois", "amphetamines", "Dangerous", "MAO-B inhibitors can increase potency unpredictably. MAO-A inhibitors can lead to hypertensive crises."),
        ("maois", "mdma", "Dangerous", "MAO-B inhibitors can increase potency unpredictably. MAO-A inhibitors with MDMA will lead to hypertensive crises."),
        ("ghb", "alcohol", "Dangerous", "Even in very low doses this combination rapidly leads to memory loss, severe ataxia and unconsciousness."),
        ("opioids", "alcohol", "Dangerous", "Both substances potentiate each other, leading to unconsciousness. Risk of vomit aspiration."),
        ("tramadol", "alcohol", "Dangerous", "Heavy CNS depressants, risk of seizures. Both substances potentiate each other."),
        ("benzodiazepines", "alcohol", "Dangerous", "Potentiate each other strongly, very rapidly leading to unconsciousness. Blacking out almost certain."),
        ("opioids", "ghb", "Dangerous", "Potentiate each other strongly and unpredictably, rapidly leading to unconsciousness."),
        ("tramadol", "ghb", "Dangerous", "The sedative effects can lead to dangerous respiratory depression."),
        ("benzodiazepines", "ghb", "Dangerous", "Potentiate each other strongly, rapidly leading to unconsciousness."),
        ("tramadol", "opioids", "Dangerous", "Increases seizure risk, CNS/respiratory depression may be synergistic."),
        ("benzodiazepines", "opioids", "Dangerous", "CNS/respiratory depression synergistic, rapidly leads to unconsciousness."),
        ("benzodiazepines", "tramadol", "Dangerous", "CNS/respiratory depression may be synergistic. Vomit aspiration risk."),
        ("maois", "ssris", "Dangerous", "High risk of serotonin syndrome."),
        ("opioids", "cocaine", "Dangerous", "Stimulants increase respiration allowing higher opioid dose. If stimulant wears off first, respiratory arrest."),
        ("maois", "cocaine", "Dangerous", "This combination is poorly explored"),
        ("maois", "pcp", "Dangerous", "This combination is very poorly explored"),
        # Low risk synergy combos
        ("mdma", "cannabis", "Low Risk & Synergy", "Can be intense as psychedelics. Consider set and setting."),
        ("mdma", "ketamine", "Low Risk & Synergy", "No unexpected interactions, may increase blood pressure."),
        ("ketamine", "dox", "Low Risk & Synergy", "Ketamine and psychedelics tend to potentiate each other - go slowly."),
        ("mdma", "amphetamines", "Low Risk & Synergy", "Amphetamines increase the neurotoxic effects of MDMA"),
        ("alcohol", "cannabis", "Low Risk & Synergy", "In excess, this combination can cause nausea."),
    ]
    
    for sub_a, sub_b, status, note in wiki_combos:
        combinations.append({
            "substance_a": sub_a,
            "substance_b": sub_b,
            "risk": RISK_LEVELS.get(status, status.lower().replace(" ", "_")),
            "status_original": status,
            "note": note,
            "source": "tripsit_wiki"
        })
    
    return combinations


def extract_polydrug_reports():
    """Extract multi-substance experiences from Erowid reports"""
    print("\n" + "="*60)
    print("EXTRACTING POLYDRUG REPORTS FROM EROWID")
    print("="*60)
    
    reports_dir = DATA_DIR / "erowid" / "reports"
    if not reports_dir.exists():
        print(f"  ✗ Erowid reports not found at {reports_dir}")
        return {}
    
    data = {
        "polydrug_reports": [],
        "combination_frequencies": defaultdict(int),
        "source": "erowid",
        "scraped_at": datetime.now().isoformat()
    }
    
    # Common substance patterns
    substance_patterns = [
        r'\b(LSD|acid)\b',
        r'\b(psilocybin|mushrooms?|shrooms?)\b',
        r'\b(DMT|ayahuasca)\b',
        r'\b(MDMA|ecstasy|molly)\b',
        r'\b(ketamine|k-hole)\b',
        r'\b(cocaine|coke)\b',
        r'\b(cannabis|marijuana|weed|THC)\b',
        r'\b(alcohol|drunk|beer|wine)\b',
        r'\b(amphetamine|adderall|speed)\b',
        r'\b(heroin|oxycodone|morphine|opioid)\b',
        r'\b(benzodiazepine|xanax|valium|benzo)\b',
        r'\b(DXM|dextromethorphan)\b',
        r'\b(mescaline|peyote|san pedro)\b',
        r'\b(2C-B|2C-E|2C-I|2C-x)\b',
        r'\b(nitrous|N2O|laughing gas)\b',
        r'\b(GHB|GBL)\b',
        r'\b(MXE|methoxetamine)\b',
        r'\b(salvia)\b',
        r'\b(caffeine|coffee)\b',
        r'\b(nicotine|cigarette|tobacco)\b',
    ]
    
    # Normalize substance names
    substance_normalize = {
        "acid": "lsd",
        "mushrooms": "psilocybin",
        "shrooms": "psilocybin",
        "ecstasy": "mdma",
        "molly": "mdma",
        "k-hole": "ketamine",
        "coke": "cocaine",
        "marijuana": "cannabis",
        "weed": "cannabis",
        "thc": "cannabis",
        "drunk": "alcohol",
        "beer": "alcohol",
        "wine": "alcohol",
        "adderall": "amphetamine",
        "speed": "amphetamine",
        "oxycodone": "opioids",
        "morphine": "opioids",
        "heroin": "opioids",
        "xanax": "benzodiazepines",
        "valium": "benzodiazepines",
        "benzo": "benzodiazepines",
        "dextromethorphan": "dxm",
        "peyote": "mescaline",
        "san pedro": "mescaline",
        "laughing gas": "nitrous",
        "n2o": "nitrous",
        "methoxetamine": "mxe",
        "coffee": "caffeine",
        "cigarette": "nicotine",
        "tobacco": "nicotine",
        "ayahuasca": "dmt_maoi",  # Special case - ayahuasca is DMT + MAOI
    }
    
    report_files = list(reports_dir.glob("*.json"))
    print(f"  Analyzing {len(report_files)} reports...")
    
    multi_substance_count = 0
    
    for i, report_file in enumerate(report_files):
        try:
            with open(report_file) as f:
                report = json.load(f)
            
            # Get substances from metadata
            substances_listed = set()
            if "substances" in report:
                for sub in report["substances"]:
                    name = sub.get("name", "").lower()
                    # Normalize
                    name = substance_normalize.get(name, name)
                    if name:
                        substances_listed.add(name)
            
            # Also scan text for substance mentions
            text = report.get("text", "").lower()
            for pattern in substance_patterns:
                matches = re.findall(pattern, text, re.IGNORECASE)
                for match in matches:
                    name = match.lower()
                    name = substance_normalize.get(name, name)
                    substances_listed.add(name)
            
            # If multiple substances, record as polydrug
            if len(substances_listed) >= 2:
                multi_substance_count += 1
                
                # Create sorted list for consistent combo keys
                substances_sorted = sorted(substances_listed)
                
                data["polydrug_reports"].append({
                    "report_id": report.get("id"),
                    "title": report.get("title"),
                    "substances": substances_sorted,
                    "substance_count": len(substances_sorted),
                    "url": report.get("url"),
                })
                
                # Count pairwise combinations
                for j, sub_a in enumerate(substances_sorted):
                    for sub_b in substances_sorted[j+1:]:
                        combo_key = f"{sub_a}+{sub_b}"
                        data["combination_frequencies"][combo_key] += 1
            
            if (i + 1) % 100 == 0:
                print(f"  [{i+1}/{len(report_files)}] Found {multi_substance_count} polydrug reports...", end="\r")
                
        except Exception as e:
            continue
    
    # Convert defaultdict to regular dict for JSON
    data["combination_frequencies"] = dict(data["combination_frequencies"])
    
    # Sort combinations by frequency
    sorted_combos = sorted(
        data["combination_frequencies"].items(),
        key=lambda x: x[1],
        reverse=True
    )
    data["top_combinations"] = sorted_combos[:50]
    
    print(f"\n  ✓ Found {multi_substance_count} polydrug reports")
    print(f"  ✓ Identified {len(data['combination_frequencies'])} unique combinations")
    
    # Show top 10
    print("\n  Top 10 most frequent combinations:")
    for combo, count in sorted_combos[:10]:
        print(f"    {combo}: {count} reports")
    
    # Save
    output_path = COMBO_DIR / "erowid_polydrug.json"
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\n  ✓ Saved to {output_path}")
    
    return data


def extract_psychonautwiki_interactions():
    """Extract interaction data from PsychonautWiki substances"""
    print("\n" + "="*60)
    print("EXTRACTING PSYCHONAUTWIKI INTERACTIONS")
    print("="*60)
    
    pw_path = DATA_DIR / "psychonautwiki" / "substances.json"
    if not pw_path.exists():
        print(f"  ✗ PsychonautWiki data not found")
        return {}
    
    with open(pw_path) as f:
        substances = json.load(f)
    
    data = {
        "interactions": [],
        "substances_with_interactions": 0,
        "source": "psychonautwiki",
        "scraped_at": datetime.now().isoformat()
    }
    
    for substance in substances:
        name = substance.get("name", "").lower()
        
        # Check for dangerous/unsafe interactions
        dangerous = substance.get("dangerousInteractions") or []
        unsafe = substance.get("unsafeInteractions") or []
        uncertain = substance.get("uncertainInteractions") or []
        
        if dangerous or unsafe or uncertain:
            data["substances_with_interactions"] += 1
        
        for other in dangerous:
            other_name = other.get("name", "").lower() if isinstance(other, dict) else str(other).lower()
            data["interactions"].append({
                "substance_a": name,
                "substance_b": other_name,
                "risk": "dangerous",
                "source": "psychonautwiki"
            })
        
        for other in unsafe:
            other_name = other.get("name", "").lower() if isinstance(other, dict) else str(other).lower()
            data["interactions"].append({
                "substance_a": name,
                "substance_b": other_name,
                "risk": "unsafe",
                "source": "psychonautwiki"
            })
        
        for other in uncertain:
            other_name = other.get("name", "").lower() if isinstance(other, dict) else str(other).lower()
            data["interactions"].append({
                "substance_a": name,
                "substance_b": other_name,
                "risk": "caution",
                "source": "psychonautwiki"
            })
    
    print(f"  ✓ Found {data['substances_with_interactions']} substances with interaction data")
    print(f"  ✓ Extracted {len(data['interactions'])} interaction records")
    
    # Save
    output_path = COMBO_DIR / "psychonautwiki_interactions.json"
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  ✓ Saved to {output_path}")
    
    return data


def build_unified_combination_database():
    """Merge all combination data into a unified database"""
    print("\n" + "="*60)
    print("BUILDING UNIFIED COMBINATION DATABASE")
    print("="*60)
    
    unified = {
        "combinations": {},  # key: "sub_a+sub_b" (sorted)
        "sources": ["mixtures.info", "tripsit", "psychonautwiki", "erowid"],
        "risk_levels": list(set(RISK_LEVELS.values())),
        "created_at": datetime.now().isoformat()
    }
    
    # Load all sources
    sources = [
        (COMBO_DIR / "mixtures_info.json", "mixtures.info"),
        (COMBO_DIR / "tripsit_combos.json", "tripsit"),
        (COMBO_DIR / "psychonautwiki_interactions.json", "psychonautwiki"),
    ]
    
    for path, source in sources:
        if path.exists():
            with open(path) as f:
                data = json.load(f)
            
            # Extract interactions based on source structure
            interactions = []
            if source == "mixtures.info":
                interactions = data.get("interactions", []) + data.get("combo_details", [])
            elif source == "tripsit":
                interactions = data.get("combinations", [])
            elif source == "psychonautwiki":
                interactions = data.get("interactions", [])
            
            for interaction in interactions:
                sub_a = interaction.get("substance_a", "").lower().strip()
                sub_b = interaction.get("substance_b", "").lower().strip()
                
                if not sub_a or not sub_b:
                    continue
                
                # Create sorted key
                key = "+".join(sorted([sub_a, sub_b]))
                
                if key not in unified["combinations"]:
                    unified["combinations"][key] = {
                        "substances": sorted([sub_a, sub_b]),
                        "sources": {},
                        "notes": []
                    }
                
                # Add source-specific data
                risk = interaction.get("risk", "unknown")
                note = interaction.get("note", "")
                
                unified["combinations"][key]["sources"][source] = risk
                if note and note not in unified["combinations"][key]["notes"]:
                    unified["combinations"][key]["notes"].append(note)
    
    # Add polydrug frequency data
    polydrug_path = COMBO_DIR / "erowid_polydrug.json"
    if polydrug_path.exists():
        with open(polydrug_path) as f:
            polydrug = json.load(f)
        
        for combo_key, freq in polydrug.get("combination_frequencies", {}).items():
            if combo_key in unified["combinations"]:
                unified["combinations"][combo_key]["erowid_report_count"] = freq
            else:
                subs = combo_key.split("+")
                if len(subs) == 2:
                    unified["combinations"][combo_key] = {
                        "substances": sorted(subs),
                        "sources": {},
                        "notes": [],
                        "erowid_report_count": freq
                    }
    
    # Calculate consensus risk level
    for key, combo in unified["combinations"].items():
        risks = list(combo["sources"].values())
        if risks:
            # Priority: dangerous > unsafe > caution > low_risk > unknown
            priority = ["dangerous", "unsafe", "caution", "low_risk_decrease", 
                       "low_risk_no_synergy", "low_risk_synergy", "unknown"]
            for level in priority:
                if level in risks:
                    combo["consensus_risk"] = level
                    break
            else:
                combo["consensus_risk"] = risks[0]
    
    # Summary statistics
    unified["statistics"] = {
        "total_combinations": len(unified["combinations"]),
        "by_risk": defaultdict(int),
        "with_erowid_reports": 0,
        "multi_source": 0
    }
    
    for combo in unified["combinations"].values():
        risk = combo.get("consensus_risk", "unknown")
        unified["statistics"]["by_risk"][risk] += 1
        
        if combo.get("erowid_report_count", 0) > 0:
            unified["statistics"]["with_erowid_reports"] += 1
        
        if len(combo.get("sources", {})) > 1:
            unified["statistics"]["multi_source"] += 1
    
    unified["statistics"]["by_risk"] = dict(unified["statistics"]["by_risk"])
    
    print(f"\n  ✓ Unified {unified['statistics']['total_combinations']} unique combinations")
    print(f"  ✓ {unified['statistics']['multi_source']} confirmed by multiple sources")
    print(f"  ✓ {unified['statistics']['with_erowid_reports']} with Erowid report data")
    
    print("\n  Risk level distribution:")
    for level, count in sorted(unified["statistics"]["by_risk"].items(), key=lambda x: -x[1]):
        print(f"    {level}: {count}")
    
    # Save
    output_path = COMBO_DIR / "unified_combinations.json"
    with open(output_path, "w") as f:
        json.dump(unified, f, indent=2)
    print(f"\n  ✓ Saved to {output_path}")
    
    return unified


def main():
    print("\n" + "="*70)
    print("COMPREHENSIVE DRUG COMBINATION DATA SCRAPER")
    print("="*70)
    
    # 1. Mixtures.info API
    mixtures_data = scrape_mixtures_api()
    
    # 2. TripSit combinations
    tripsit_data = scrape_tripsit_combos()
    
    # 3. PsychonautWiki interactions
    pw_data = extract_psychonautwiki_interactions()
    
    # 4. Erowid polydrug reports
    erowid_data = extract_polydrug_reports()
    
    # 5. Build unified database
    unified = build_unified_combination_database()
    
    # Final summary
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    print(f"""
Data Collected:
  • Mixtures.info: {len(mixtures_data.get('substances', []))} substances, {len(mixtures_data.get('interactions', []))} interactions
  • TripSit: {tripsit_data.get('substances_with_combos', 0)} substances with combos, {len(tripsit_data.get('combinations', []))} records
  • PsychonautWiki: {pw_data.get('substances_with_interactions', 0)} substances with interactions
  • Erowid: {len(erowid_data.get('polydrug_reports', []))} polydrug reports

Unified Database:
  • {unified['statistics']['total_combinations']} unique substance combinations
  • {unified['statistics']['multi_source']} confirmed by multiple sources
  • {unified['statistics']['with_erowid_reports']} with real-world experience reports

Files saved to: {COMBO_DIR}
    """)


if __name__ == "__main__":
    main()
