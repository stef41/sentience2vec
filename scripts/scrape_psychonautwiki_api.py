#!/usr/bin/env python3
"""
Scrape drug -> effect relationships from PsychonautWiki
This creates a mapping of which drugs produce which subjective effects.
"""

import requests
import json
import time
import re
from pathlib import Path

# PsychonautWiki API endpoint
API_URL = "https://api.psychonautwiki.org/graphql"

# First get all substances with their effects
SUBSTANCES_QUERY = """
query {
    substances(limit: 500) {
        name
        class {
            chemical
            psychoactive
        }
        effects {
            name
            url
        }
    }
}
"""

OUTPUT_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

print("Fetching drug data from PsychonautWiki API...")

try:
    response = requests.post(
        API_URL,
        json={'query': SUBSTANCES_QUERY},
        headers={'Content-Type': 'application/json'},
        timeout=60
    )
    response.raise_for_status()
    data = response.json()
    
    substances = data.get('data', {}).get('substances', [])
    print(f"Retrieved {len(substances)} substances from PsychonautWiki API")
    
    # Process into drug -> effects mapping
    drug_effects = {}
    effect_drugs = {}  # Also build reverse mapping: effect -> drugs
    
    for substance in substances:
        if not substance:
            continue
            
        name = substance.get('name', '')
        if not name:
            continue
            
        effects = substance.get('effects', []) or []
        drug_class = substance.get('class', {})
        
        psychoactive_class = None
        if drug_class:
            psychoactive = drug_class.get('psychoactive', []) or []
            if psychoactive:
                psychoactive_class = psychoactive[0] if isinstance(psychoactive, list) else psychoactive
        
        effect_names = []
        for effect in effects:
            if effect and effect.get('name'):
                effect_name = effect['name']
                effect_names.append(effect_name)
                
                # Build reverse mapping
                if effect_name not in effect_drugs:
                    effect_drugs[effect_name] = []
                if name not in effect_drugs[effect_name]:
                    effect_drugs[effect_name].append(name)
        
        drug_effects[name] = {
            'effects': effect_names,
            'class': psychoactive_class,
            'source': 'PsychonautWiki'
        }
    
    print(f"\nProcessed {len(drug_effects)} drugs with effects")
    print(f"Found {len(effect_drugs)} unique effects")
    
    # Save drug -> effects mapping
    output_file = OUTPUT_PATH / 'psychonautwiki_drug_effects.json'
    with open(output_file, 'w') as f:
        json.dump(drug_effects, f, indent=2)
    print(f"\nSaved drug effects to {output_file}")
    
    # Save effect -> drugs mapping
    output_file2 = OUTPUT_PATH / 'psychonautwiki_effect_drugs.json'
    with open(output_file2, 'w') as f:
        json.dump(effect_drugs, f, indent=2)
    print(f"Saved effect-drug mapping to {output_file2}")
    
    # Print some stats
    print("\n=== Effect Statistics ===")
    sorted_effects = sorted(effect_drugs.items(), key=lambda x: len(x[1]), reverse=True)
    for effect, drugs in sorted_effects[:30]:
        print(f"  {effect}: {len(drugs)} drugs")
    
    # Print sample drug data
    print("\n=== Sample Drug Data ===")
    for drug, data in list(drug_effects.items())[:10]:
        print(f"\n{drug}:")
        print(f"  Class: {data['class']}")
        print(f"  Effects ({len(data['effects'])}): {', '.join(data['effects'][:5])}...")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
