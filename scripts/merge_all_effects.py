#!/usr/bin/env python3
"""
Merge effect-drug relationships from both PsychonautWiki and TripSit.
Creates the final scraped_relationships.js with data from both sources.
"""

import json
import re
from pathlib import Path
from datetime import datetime
from collections import defaultdict

DATA_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')

print("Loading data from both sources...")

# Load PsychonautWiki effect->drug data (scraped from effect pages)
with open(DATA_PATH / 'effect_drugs_direct.json') as f:
    pwiki_effects = json.load(f)
print(f"  PsychonautWiki: {len(pwiki_effects)} effects")

# Load TripSit effect data
with open(DATA_PATH / 'tripsit_effects.json') as f:
    tripsit_data = json.load(f)
    tripsit_effect_drugs = tripsit_data['effect_drugs']
print(f"  TripSit: {len(tripsit_effect_drugs)} effects")

# Load TripSit drug classifications
with open(DATA_PATH / 'tripsit_data.json') as f:
    tripsit_drugs = json.load(f)
print(f"  TripSit drugs: {len(tripsit_drugs)}")

# Merge effect->drug mappings
# PsychonautWiki is primary, TripSit supplements
merged_effects = {}

# First add all PsychonautWiki effects
for effect_id, data in pwiki_effects.items():
    merged_effects[effect_id] = {
        'name': data['name'],
        'drugs': list(data['drugs']),
        'url': data['url'],
        'source': 'PsychonautWiki'
    }

# Then add/merge TripSit effects
for effect_id, drugs in tripsit_effect_drugs.items():
    if effect_id in merged_effects:
        # Merge drugs (add new ones from TripSit)
        existing_drugs = set(d.lower() for d in merged_effects[effect_id]['drugs'])
        for drug in drugs:
            if drug.lower() not in existing_drugs:
                merged_effects[effect_id]['drugs'].append(drug)
        merged_effects[effect_id]['source'] = 'PsychonautWiki + TripSit'
    else:
        # New effect from TripSit
        # Format the name nicely
        nice_name = effect_id.replace('_', ' ').title()
        merged_effects[effect_id] = {
            'name': nice_name,
            'drugs': drugs,
            'url': f'https://tripsit.me/wiki/{effect_id}',
            'source': 'TripSit'
        }

print(f"\nMerged: {len(merged_effects)} total effects")
print(f"Total relationships: {sum(len(e['drugs']) for e in merged_effects.values())}")

# Map TripSit categories to drug classes
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

def get_drug_class(entry):
    categories = entry.get('categories', [])
    priority = ['psychedelic', 'dissociative', 'empathogen', 'stimulant', 
                'opioid', 'depressant', 'deliriant', 'nootropic', 'cannabinoid']
    for p in priority:
        if p in categories:
            return CATEGORY_MAP.get(p, p)
    for cat in categories:
        if cat in CATEGORY_MAP:
            return CATEGORY_MAP[cat]
    return 'other'

# Generate JavaScript
js_code = f"""
// ========================================
// Effect-Drug Relationships
// Merged from PsychonautWiki and TripSit
// Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}
// ========================================
// Sources:
//   - PsychonautWiki: https://psychonautwiki.org/wiki/Subjective_effect_index
//   - TripSit: https://tripsit.me
// ========================================

const EFFECT_DRUG_RELATIONSHIPS = {{
"""

for effect_id, data in sorted(merged_effects.items()):
    drugs_json = json.dumps(data['drugs'])
    js_code += f'  "{effect_id}": {drugs_json},\n'

js_code += """}};

// Effect metadata (names, URLs, sources)
const EFFECT_METADATA = {
"""

for effect_id, data in sorted(merged_effects.items()):
    name = data['name'].replace('"', '\\"')
    url = data.get('url', '')
    source = data.get('source', 'Unknown')
    js_code += f'  "{effect_id}": {{"name": "{name}", "url": "{url}", "source": "{source}"}},\n'

js_code += """};

// Function to find drugs that produce a given effect
function getDrugsForEffectFromScrapedData(effectId) {
    var normalized = effectId.toLowerCase()
        .replace(/ /g, '_')
        .replace(/-/g, '_')
        .replace(/[^a-z0-9_]/g, '');
    
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
    
    return matches.filter(function(item, pos) {
        return matches.indexOf(item) === pos;
    });
}

// Get effect metadata
function getEffectMetadata(effectId) {
    var normalized = effectId.toLowerCase()
        .replace(/ /g, '_')
        .replace(/-/g, '_')
        .replace(/[^a-z0-9_]/g, '');
    return EFFECT_METADATA[normalized] || null;
}

// ========================================
// Drug Classifications from TripSit
// Source: https://tripsit.me
// ========================================

const TRIPSIT_DRUG_DATA = {
"""

for drug_id, data in sorted(tripsit_drugs.items()):
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

# Stats
pwiki_only = sum(1 for e in merged_effects.values() if e['source'] == 'PsychonautWiki')
tripsit_only = sum(1 for e in merged_effects.values() if e['source'] == 'TripSit')
both = sum(1 for e in merged_effects.values() if e['source'] == 'PsychonautWiki + TripSit')

print(f"\n=== Source Breakdown ===")
print(f"  PsychonautWiki only: {pwiki_only} effects")
print(f"  TripSit only: {tripsit_only} effects")
print(f"  Both sources: {both} effects")

# Verification
print("\n=== Verification ===")
test_effects = ['visual_processing_deceleration', 'euphoria', 'sedation', 'time_distortion', 'anxiety_suppression']
for effect in test_effects:
    if effect in merged_effects:
        data = merged_effects[effect]
        print(f"\n{data['name']} ({data['source']}):")
        print(f"  {', '.join(data['drugs'][:8])}" + (f"... (+{len(data['drugs'])-8})" if len(data['drugs']) > 8 else ""))
