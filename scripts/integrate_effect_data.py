#!/usr/bin/env python3
"""
Integrate scraped drug-effect data from PsychonautWiki and TripSit into the website.
Uses the direct effect->drug mappings from effect pages.
"""

import json
import re
from pathlib import Path
from datetime import datetime

DATA_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')

# Load the direct effect->drug mappings (scraped from effect pages)
print("Loading scraped data...")
with open(DATA_PATH / 'effect_drugs_direct.json') as f:
    effect_drug_data = json.load(f)
print(f"  Effect->Drug mappings: {len(effect_drug_data)} effects")

# Load TripSit data for drug classifications
with open(DATA_PATH / 'tripsit_data.json') as f:
    tripsit_data = json.load(f)
print(f"  TripSit drugs: {len(tripsit_data)}")

# Map TripSit categories to our class names
CATEGORY_MAP = {
    'psychedelic': 'psychedelic',
    'dissociative': 'dissociative', 
    'stimulant': 'stimulant',
    'depressant': 'depressant',
    'opioid': 'opioid',
    'empathogen': 'empathogen',
    'entactogen': 'empathogen',
    'benzodiazepine': 'depressant',
    'barbiturate': 'depressant',
    'deliriant': 'deliriant',
    'nootropic': 'nootropic',
    'cannabinoid': 'cannabinoid',
}

def get_drug_class(tripsit_entry):
    """Get standardized drug class from TripSit data."""
    categories = tripsit_entry.get('categories', [])
    priority = ['psychedelic', 'dissociative', 'empathogen', 'stimulant', 
                'opioid', 'depressant', 'deliriant', 'nootropic', 'cannabinoid']
    
    for p in priority:
        if p in categories:
            return CATEGORY_MAP.get(p, p)
    for cat in categories:
        if cat in CATEGORY_MAP:
            return CATEGORY_MAP[cat]
    return 'other'

# Generate JavaScript code
js_code = f"""
// ========================================
// Effect-Drug Relationships from PsychonautWiki
// Scraped directly from effect pages on {datetime.now().strftime('%Y-%m-%d')}
// Source: https://psychonautwiki.org/wiki/Subjective_effect_index
// ========================================

const EFFECT_DRUG_RELATIONSHIPS = {{
"""

for effect_id, data in sorted(effect_drug_data.items()):
    drugs_json = json.dumps(data['drugs'])
    js_code += f'  "{effect_id}": {drugs_json},\n'

js_code += """}};

// Effect metadata (names and URLs)
const EFFECT_METADATA = {
"""

for effect_id, data in sorted(effect_drug_data.items()):
    name = data['name'].replace('"', '\\"')
    url = data['url']
    js_code += f'  "{effect_id}": {{"name": "{name}", "url": "{url}"}},\n'

js_code += """};

// Function to find drugs that produce a given effect
function getDrugsForEffectFromPsychonautWiki(effectId) {
    // Normalize the effect ID
    var normalized = effectId.toLowerCase()
        .replace(/ /g, '_')
        .replace(/-/g, '_')
        .replace(/[^a-z0-9_]/g, '');
    
    // Direct lookup
    if (EFFECT_DRUG_RELATIONSHIPS[normalized]) {
        return EFFECT_DRUG_RELATIONSHIPS[normalized];
    }
    
    // Try partial matching
    var matches = [];
    Object.keys(EFFECT_DRUG_RELATIONSHIPS).forEach(function(key) {
        if (key.indexOf(normalized) !== -1 || normalized.indexOf(key) !== -1) {
            matches = matches.concat(EFFECT_DRUG_RELATIONSHIPS[key]);
        }
    });
    
    // Remove duplicates
    return matches.filter(function(item, pos) {
        return matches.indexOf(item) === pos;
    });
}

// ========================================
// Drug Classifications from TripSit
// Source: https://tripsit.me
// ========================================

const TRIPSIT_DRUG_DATA = {
"""

for drug_id, data in sorted(tripsit_data.items()):
    name = data['name'].replace('"', '\\"')
    categories = data.get('categories', [])
    drug_class = get_drug_class(data)
    aliases = data.get('aliases', [])[:5]
    
    js_code += f'  "{name}": {{\n'
    js_code += f'    "class": "{drug_class}",\n'
    js_code += f'    "categories": {json.dumps(categories)},\n'
    js_code += f'    "aliases": {json.dumps(aliases)},\n'
    js_code += f'    "source": "TripSit"\n'
    js_code += '  },\n'

js_code += "};\n"

# Save
output_path = DATA_PATH / 'scraped_relationships.js'
with open(output_path, 'w') as f:
    f.write(js_code)

print(f"\nSaved to {output_path}")
print(f"\nStats:")
print(f"  Total effects: {len(effect_drug_data)}")
print(f"  Total drug-effect relationships: {sum(len(d['drugs']) for d in effect_drug_data.values())}")

# Show sample
print("\n=== Verification ===")
test_effects = ['visual_processing_deceleration', 'euphoria', 'sedation', 'time_distortion', 'ego_death']
for effect in test_effects:
    if effect in effect_drug_data:
        data = effect_drug_data[effect]
        print(f"\n{data['name']}:")
        print(f"  {', '.join(data['drugs'][:10])}" + (f"... (+{len(data['drugs'])-10})" if len(data['drugs']) > 10 else ""))
