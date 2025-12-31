#!/usr/bin/env python3
"""
Process Ferrari/Berkeley Erowid Trip Report Dataset
====================================================
4,910 trip narratives from UC Berkeley Center for the Science of Psychedelics
Source: github.com/jonathanferrari/scraper

This script:
1. Loads all narrative text files with metadata
2. Extracts polydrug combinations from reports
3. Analyzes effect descriptions
4. Builds substance frequency statistics
5. Exports unified format for consciousness mapping
"""

import json
import re
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import pandas as pd

# Paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "raw" / "ferrari_narratives"
OUTPUT_DIR = BASE_DIR / "data" / "raw" / "ferrari_processed"

# Substance name normalization
SUBSTANCE_ALIASES = {
    # Psychedelics
    'lsd': 'lsd', 'acid': 'lsd', 'lsd-25': 'lsd',
    'mushrooms': 'psilocybin', 'mushrooms - p. cubensis': 'psilocybin',
    'mushrooms - p. semilanceata': 'psilocybin', 'mushrooms - p. mexicana': 'psilocybin',
    'mushrooms (p. cubensis)': 'psilocybin', 'mushrooms - p. cubensis': 'psilocybin',
    'psilocybin mushrooms': 'psilocybin', 'magic mushrooms': 'psilocybin',
    'shrooms': 'psilocybin',
    'dmt': 'dmt', 'n,n-dmt': 'dmt', 'nn-dmt': 'dmt',
    '5-meo-dmt': '5-meo-dmt', '5meodmt': '5-meo-dmt',
    'mescaline': 'mescaline',
    'cacti - t. pachanoi': 'mescaline', 'cacti - t. peruvianus': 'mescaline',
    'cacti - t. bridgesii': 'mescaline', 'peyote': 'mescaline',
    'san pedro': 'mescaline',
    '2c-b': '2c-b', '2cb': '2c-b',
    '2c-e': '2c-e', '2c-i': '2c-i',
    'ayahuasca': 'ayahuasca', 'aya': 'ayahuasca',
    'salvia': 'salvia', 'salvia divinorum': 'salvia',
    
    # Empathogens
    'mdma': 'mdma', 'mdma (ecstasy)': 'mdma', 'ecstasy': 'mdma',
    'molly': 'mdma', 'e': 'mdma', 'x': 'mdma',
    'mda': 'mda',
    
    # Dissociatives
    'ketamine': 'ketamine', 'ket': 'ketamine', 'k': 'ketamine',
    'dxm': 'dxm', 'dextromethorphan': 'dxm', 'robitussin': 'dxm',
    'pcp': 'pcp', 'phencyclidine': 'pcp',
    'nitrous oxide': 'nitrous', 'nitrous': 'nitrous', 'n2o': 'nitrous', 'nos': 'nitrous',
    
    # Cannabis
    'cannabis': 'cannabis', 'marijuana': 'cannabis', 'weed': 'cannabis',
    'cannabis (edible)': 'cannabis', 'cannabis - hash': 'cannabis',
    'thc': 'cannabis', 'pot': 'cannabis', 'herb': 'cannabis',
    
    # Opioids
    'heroin': 'heroin', 'morphine': 'morphine', 'codeine': 'codeine',
    'oxycodone': 'oxycodone', 'hydrocodone': 'hydrocodone',
    'kratom': 'kratom', 'opium': 'opium',
    
    # Stimulants
    'cocaine': 'cocaine', 'coke': 'cocaine',
    'amphetamine': 'amphetamine', 'adderall': 'amphetamine',
    'methamphetamine': 'methamphetamine', 'meth': 'methamphetamine',
    'caffeine': 'caffeine', 'coffee': 'caffeine',
    
    # Depressants
    'alcohol': 'alcohol', 'beer': 'alcohol', 'wine': 'alcohol',
    'ghb': 'ghb', 'gbl': 'ghb',
    'alprazolam': 'alprazolam', 'xanax': 'alprazolam',
    'diazepam': 'diazepam', 'valium': 'diazepam',
    
    # Others
    'ibogaine': 'ibogaine', 'iboga': 'ibogaine',
    'dph': 'diphenhydramine', 'diphenhydramine': 'diphenhydramine',
    'benadryl': 'diphenhydramine',
}

# Effect keywords to search for
EFFECT_KEYWORDS = {
    'visual': [
        'visual', 'hallucination', 'pattern', 'geometric', 'fractal', 'color',
        'morphing', 'breathing', 'trail', 'tracers', 'distort', 'warp', 'shimmer'
    ],
    'auditory': [
        'sound', 'music', 'hear', 'auditory', 'echo', 'ringing', 'voice'
    ],
    'cognitive': [
        'think', 'thought', 'understand', 'insight', 'realize', 'consciousness',
        'awareness', 'mind', 'ego', 'identity', 'meaning', 'profound'
    ],
    'emotional': [
        'feel', 'emotion', 'love', 'fear', 'anxiety', 'euphoria', 'happy',
        'sad', 'peace', 'terror', 'joy', 'bliss', 'panic', 'empathy'
    ],
    'physical': [
        'body', 'nausea', 'heart', 'tingle', 'numb', 'warm', 'cold',
        'shake', 'tremble', 'sweat', 'dizzy', 'heavy', 'light'
    ],
    'temporal': [
        'time', 'eternal', 'loop', 'forever', 'hour', 'minute', 'slow', 'fast'
    ],
    'mystical': [
        'god', 'spirit', 'universe', 'cosmic', 'transcend', 'divine',
        'sacred', 'infinite', 'eternal', 'entity', 'being', 'death', 'rebirth'
    ],
    'dissociative': [
        'dissociate', 'disconnect', 'separate', 'float', 'detach', 'hole',
        'void', 'nothing', 'ego death', 'dissolve'
    ],
    'synesthesia': [
        'synesthesia', 'taste color', 'hear color', 'see sound', 'smell music'
    ]
}


def normalize_substance(name: str) -> str:
    """Normalize substance name to canonical form."""
    name_lower = name.lower().strip()
    return SUBSTANCE_ALIASES.get(name_lower, name_lower)


def load_metadata() -> pd.DataFrame:
    """Load the info.csv metadata file."""
    info_path = DATA_DIR / "narratives" / "csv" / "info.csv"
    return pd.read_csv(info_path)


def load_narrative(report_id: int) -> str | None:
    """Load narrative text for a given report ID."""
    txt_path = DATA_DIR / "narratives" / "txt" / f"{report_id}.txt"
    if txt_path.exists():
        try:
            return txt_path.read_text(encoding='utf-8', errors='ignore')
        except Exception:
            return None
    return None


def extract_mentioned_substances(text: str) -> list[str]:
    """Extract substance names mentioned in the narrative text."""
    text_lower = text.lower()
    mentioned = set()
    
    # Check for all known substance aliases
    for alias, canonical in SUBSTANCE_ALIASES.items():
        # Use word boundaries
        pattern = r'\b' + re.escape(alias) + r'\b'
        if re.search(pattern, text_lower):
            mentioned.add(canonical)
    
    return list(mentioned)


def extract_effect_categories(text: str) -> dict[str, int]:
    """Count effect keyword occurrences by category."""
    text_lower = text.lower()
    effect_counts = {}
    
    for category, keywords in EFFECT_KEYWORDS.items():
        count = 0
        for kw in keywords:
            count += len(re.findall(r'\b' + re.escape(kw) + r'\b', text_lower))
        effect_counts[category] = count
    
    return effect_counts


def process_all_narratives():
    """Process all narratives and build unified dataset."""
    print("=" * 70)
    print("PROCESSING FERRARI/BERKELEY TRIP REPORTS")
    print("=" * 70)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load metadata
    print("\n📊 Loading metadata...")
    metadata = load_metadata()
    print(f"   Found {len(metadata)} reports with metadata")
    
    # Process each narrative
    processed = []
    polydrug_reports = []
    substance_effects = defaultdict(lambda: defaultdict(int))
    combination_frequencies = Counter()
    total_text_size = 0
    
    print(f"\n📖 Processing narratives...")
    for idx, row in metadata.iterrows():
        if idx % 500 == 0:
            print(f"   Processing {idx}/{len(metadata)}...")
        
        report_id = row['id']
        text = load_narrative(report_id)
        
        if not text:
            continue
        
        total_text_size += len(text)
        
        # Primary substance from metadata
        primary_substance = normalize_substance(row['substance_name'])
        
        # Find all mentioned substances in text
        mentioned = extract_mentioned_substances(text)
        if primary_substance not in mentioned:
            mentioned.append(primary_substance)
        
        # Extract effect categories
        effects = extract_effect_categories(text)
        
        # Accumulate effects for this substance
        for category, count in effects.items():
            substance_effects[primary_substance][category] += count
        
        # Build processed record
        record = {
            'id': report_id,
            'primary_substance': primary_substance,
            'original_substance': row['substance_name'],
            'trip_type': row['trip_type'],
            'intensity': row.get('intensity', 3),
            'rating': row.get('rating', ''),
            'gender': row.get('gender', ''),
            'experience_year': row.get('experience-year', ''),
            'title': row.get('title', ''),
            'text_length': len(text),
            'mentioned_substances': mentioned,
            'effect_scores': effects,
            'is_polydrug': len(mentioned) > 1
        }
        processed.append(record)
        
        # Track polydrug reports
        if len(mentioned) > 1:
            combo_key = '+'.join(sorted(mentioned))
            combination_frequencies[combo_key] += 1
            polydrug_reports.append({
                'id': report_id,
                'substances': mentioned,
                'primary': primary_substance,
                'combo_key': combo_key,
                'trip_type': row['trip_type'],
                'intensity': row.get('intensity', 3),
                'text_length': len(text)
            })
    
    print(f"   Processed {len(processed)} reports")
    print(f"   Total text corpus: {total_text_size:,} chars ({total_text_size // (1024*1024):.1f} MB)")
    
    # Statistics
    print("\n" + "=" * 70)
    print("STATISTICS")
    print("=" * 70)
    
    # Substance distribution
    substance_counts = Counter(r['primary_substance'] for r in processed)
    print(f"\n📊 Top 15 Substances:")
    for sub, count in substance_counts.most_common(15):
        print(f"   {sub:25s}: {count:4d} reports")
    
    # Trip type distribution
    trip_types = Counter(r['trip_type'] for r in processed)
    print(f"\n📊 Trip Types:")
    for tt, count in trip_types.items():
        print(f"   {tt:15s}: {count:4d} reports")
    
    # Polydrug statistics
    print(f"\n📊 Polydrug Reports:")
    print(f"   Total polydrug reports: {len(polydrug_reports)}")
    print(f"   Unique combinations: {len(combination_frequencies)}")
    print(f"\n   Top 15 combinations:")
    for combo, count in combination_frequencies.most_common(15):
        print(f"   {combo:50s}: {count:4d}")
    
    # Effect category by substance
    print(f"\n📊 Effect Categories by Substance (top 5):")
    for substance in list(substance_counts.keys())[:5]:
        effects = substance_effects[substance]
        top_effects = sorted(effects.items(), key=lambda x: -x[1])[:3]
        print(f"   {substance}: {', '.join(f'{e}({c})' for e, c in top_effects)}")
    
    # Save outputs
    print("\n" + "=" * 70)
    print("SAVING OUTPUTS")
    print("=" * 70)
    
    # 1. All processed reports
    output_file = OUTPUT_DIR / "all_reports.json"
    with open(output_file, 'w') as f:
        json.dump(processed, f, indent=2, default=str)
    print(f"✅ Saved {len(processed)} reports to {output_file.name}")
    
    # 2. Polydrug reports
    polydrug_file = OUTPUT_DIR / "polydrug_reports.json"
    with open(polydrug_file, 'w') as f:
        json.dump({
            'reports': polydrug_reports,
            'combination_frequencies': dict(combination_frequencies.most_common()),
            'total_polydrug': len(polydrug_reports),
            'unique_combinations': len(combination_frequencies)
        }, f, indent=2)
    print(f"✅ Saved {len(polydrug_reports)} polydrug reports to {polydrug_file.name}")
    
    # 3. Substance effect profiles
    effect_profiles = {}
    for substance, effects in substance_effects.items():
        total = sum(effects.values())
        if total > 0:
            effect_profiles[substance] = {
                'raw_counts': dict(effects),
                'total_mentions': total,
                'report_count': substance_counts[substance],
                'normalized': {k: round(v / total, 3) for k, v in effects.items()}
            }
    
    profiles_file = OUTPUT_DIR / "substance_effect_profiles.json"
    with open(profiles_file, 'w') as f:
        json.dump(effect_profiles, f, indent=2)
    print(f"✅ Saved {len(effect_profiles)} effect profiles to {profiles_file.name}")
    
    # 4. Summary statistics
    summary = {
        'source': 'jonathanferrari/scraper (UC Berkeley)',
        'processed_at': datetime.now().isoformat(),
        'total_reports': len(processed),
        'total_text_chars': total_text_size,
        'total_text_mb': round(total_text_size / (1024 * 1024), 2),
        'unique_substances': len(substance_counts),
        'polydrug_reports': len(polydrug_reports),
        'unique_combinations': len(combination_frequencies),
        'trip_types': dict(trip_types),
        'substance_counts': dict(substance_counts.most_common()),
        'top_combinations': dict(combination_frequencies.most_common(50))
    }
    
    summary_file = OUTPUT_DIR / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"✅ Saved summary to {summary_file.name}")
    
    # Final summary
    print("\n" + "=" * 70)
    print("FERRARI/BERKELEY DATASET PROCESSING COMPLETE")
    print("=" * 70)
    print(f"""
📊 SUMMARY:
   • Total reports: {len(processed):,}
   • Text corpus: {total_text_size // (1024*1024):.1f} MB
   • Unique substances: {len(substance_counts)}
   • Polydrug reports: {len(polydrug_reports):,}
   • Unique combinations: {len(combination_frequencies)}
   
📂 OUTPUT FILES:
   • {OUTPUT_DIR / 'all_reports.json'}
   • {OUTPUT_DIR / 'polydrug_reports.json'}
   • {OUTPUT_DIR / 'substance_effect_profiles.json'}
   • {OUTPUT_DIR / 'summary.json'}
""")


if __name__ == "__main__":
    process_all_narratives()
