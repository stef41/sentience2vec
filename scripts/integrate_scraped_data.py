#!/usr/bin/env python3
"""
Integrate scraped drug data from TripSit and PsychonautWiki into data.js.
Creates proper drug -> effect mappings from real sources.
"""

import json
import re
from pathlib import Path

DATA_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')
OUTPUT_FILE = Path('/data/users/zacharie/whc/consciousness_study/website/data.js')

# Load scraped data
print("Loading scraped data...")

with open(DATA_PATH / 'tripsit_data.json') as f:
    tripsit_data = json.load(f)
print(f"  TripSit: {len(tripsit_data)} drugs")

with open(DATA_PATH / 'psychonautwiki_drug_effects.json') as f:
    pwiki_drugs = json.load(f)
print(f"  PsychonautWiki: {len(pwiki_drugs)} drugs with effects")

with open(DATA_PATH / 'psychonautwiki_effect_drugs.json') as f:
    pwiki_effects = json.load(f)
print(f"  PsychonautWiki: {len(pwiki_effects)} unique effects")

# Clean up effect names - remove non-effect items
EXCLUDE_EFFECTS = {
    'contributors', 'Subjective Effect Index', 'experience index',
    'PsychonautWiki', 'edit', 'source', 'citation', 'reference'
}

for effect in list(pwiki_effects.keys()):
    if any(x in effect.lower() for x in ['contributor', 'index', 'wiki', 'edit', 'source']):
        del pwiki_effects[effect]

print(f"  After cleanup: {len(pwiki_effects)} effects")

# Create normalized name lookup for drugs
def normalize(name):
    return re.sub(r'[^a-z0-9]', '', name.lower())

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

# Build comprehensive drug database
print("\nBuilding integrated drug database...")

# Create lookup indices
tripsit_lookup = {}
for drug_id, data in tripsit_data.items():
    norm = normalize(data['name'])
    tripsit_lookup[norm] = data
    for alias in data.get('aliases', []):
        tripsit_lookup[normalize(alias)] = data

pwiki_lookup = {}
for drug_name, data in pwiki_drugs.items():
    pwiki_lookup[normalize(drug_name)] = data

# Build effect_id -> drugs mapping
# Clean effect names and create IDs
def effect_to_id(effect_name):
    """Convert effect name to ID format."""
    # Remove quotes and special characters
    clean = effect_name.replace('"', '').replace("'", '').strip()
    return clean.lower().replace(' ', '_').replace('-', '_')

effect_drug_mapping = {}
for effect_name, drugs in pwiki_effects.items():
    # Skip non-effects
    if len(effect_name) < 3 or len(effect_name) > 60:
        continue
    if any(x in effect_name.lower() for x in ['http', '[', ']', 'edit']):
        continue
    
    effect_id = effect_to_id(effect_name)
    effect_drug_mapping[effect_id] = {
        'name': effect_name,
        'drugs': drugs,
        'source': 'PsychonautWiki'
    }

print(f"Built {len(effect_drug_mapping)} effect -> drugs mappings")

# Generate JavaScript code for effect-drug relationships
js_code = """
// ========================================
// Effect-Drug Relationships from PsychonautWiki
// Scraped on """ + str(__import__('datetime').datetime.now().strftime('%Y-%m-%d')) + """
// Source: https://psychonautwiki.org
// ========================================

const EFFECT_DRUG_RELATIONSHIPS = {
"""

for effect_id, data in sorted(effect_drug_mapping.items()):
    drugs_list = json.dumps(data['drugs'])
    js_code += f'  "{effect_id}": {drugs_list},\n'

js_code += """};

// Function to find drugs that produce a given effect
function getDrugsForEffectFromPsychonautWiki(effectId) {
    const normalizedId = effectId.toLowerCase().replace(/ /g, '_').replace(/-/g, '_');
    return EFFECT_DRUG_RELATIONSHIPS[normalizedId] || [];
}

// ========================================
// Drug Classifications from TripSit
// Source: https://tripsit.me
// ========================================

const TRIPSIT_DRUG_DATA = {
"""

for drug_id, data in sorted(tripsit_data.items()):
    name = data['name']
    categories = data.get('categories', [])
    drug_class = get_drug_class(data)
    aliases = data.get('aliases', [])[:5]  # Limit aliases
    
    js_code += f'  "{name}": {{\n'
    js_code += f'    "class": "{drug_class}",\n'
    js_code += f'    "categories": {json.dumps(categories)},\n'
    js_code += f'    "aliases": {json.dumps(aliases)},\n'
    js_code += f'    "source": "TripSit"\n'
    js_code += '  },\n'

js_code += "};\n"

# Save to separate file
output_path = DATA_PATH / 'scraped_relationships.js'
with open(output_path, 'w') as f:
    f.write(js_code)

print(f"\nSaved relationships to {output_path}")

# Now update findDrugsForEffect in brain3d.js to use this data
print("\nStats:")
print(f"  Total effect-drug relationships: {sum(len(d['drugs']) for d in effect_drug_mapping.values())}")
print(f"  Drugs in TripSit: {len(tripsit_data)}")

# Show sample mappings
print("\n=== Sample Effect -> Drugs ===")
samples = ['euphoria', 'time_distortion', 'visual_hallucination', 'sedation', 'stimulation', 'anxiety']
for s in samples:
    for eid, data in effect_drug_mapping.items():
        if s in eid:
            print(f"  {data['name']}: {data['drugs'][:5]}...")
            break
