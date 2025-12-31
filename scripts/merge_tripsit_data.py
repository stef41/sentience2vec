#!/usr/bin/env python3
"""
Merge TripSit drug classifications into the website data.js
Uses real data from TripSit API instead of guessing.
"""

import json
import re

# Load TripSit data
with open('/data/users/zacharie/whc/consciousness_study/website/data/tripsit_data.json', 'r') as f:
    tripsit_data = json.load(f)

print(f"Loaded {len(tripsit_data)} drugs from TripSit")

# Create name lookup (lowercase, no special chars)
def normalize_name(name):
    return re.sub(r'[^a-z0-9]', '', name.lower())

tripsit_lookup = {}
for drug_id, data in tripsit_data.items():
    # Add main name
    norm = normalize_name(data['name'])
    tripsit_lookup[norm] = data
    
    # Add aliases
    for alias in data.get('aliases', []):
        norm_alias = normalize_name(alias)
        if norm_alias not in tripsit_lookup:
            tripsit_lookup[norm_alias] = data

print(f"Built lookup with {len(tripsit_lookup)} entries")

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

def get_class_from_tripsit(data):
    """Get the primary drug class from TripSit categories."""
    categories = data.get('categories', [])
    
    # Priority order for class assignment
    priority = ['psychedelic', 'dissociative', 'empathogen', 'stimulant', 
                'opioid', 'depressant', 'deliriant', 'nootropic', 'cannabinoid']
    
    for p in priority:
        if p in categories:
            return CATEGORY_MAP.get(p, p)
    
    # Check for mapped categories
    for cat in categories:
        if cat in CATEGORY_MAP:
            return CATEGORY_MAP[cat]
    
    return None

# Read data.js
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'r') as f:
    content = f.read()

# Find drugs section
drugs_match = re.search(r'drugs:\s*\{', content)
if not drugs_match:
    print("Could not find drugs section!")
    exit(1)

# Extract drug names and update classes
updated = 0
not_found = []

# Find all drug entries
drug_pattern = re.compile(r'"([^"]+)":\s*\{\s*"receptors"')
matches = list(drug_pattern.finditer(content, drugs_match.end()))

for match in matches:
    drug_name = match.group(1)
    
    # Skip entries that are clearly not drug names
    if any(x in drug_name.lower() for x in ['alpha', 'beta', 'glun', 'gluk', 'glua', 'muscle', 'sigma']):
        continue
    
    norm_name = normalize_name(drug_name)
    
    # Look up in TripSit
    tripsit_info = tripsit_lookup.get(norm_name)
    
    if tripsit_info:
        new_class = get_class_from_tripsit(tripsit_info)
        
        if new_class:
            # Find the class field for this drug
            drug_start = match.start()
            class_pattern = re.compile(r'"class":\s*"([^"]*)"')
            class_match = class_pattern.search(content, drug_start)
            
            if class_match and class_match.start() < drug_start + 10000:
                old_class = class_match.group(1)
                if old_class != new_class:
                    content = content[:class_match.start()] + f'"class": "{new_class}"' + content[class_match.end():]
                    updated += 1
                    print(f"  {drug_name}: {old_class} -> {new_class} (TripSit: {tripsit_info['categories']})")
    else:
        not_found.append(drug_name)

# Write back
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'w') as f:
    f.write(content)

print(f"\nUpdated {updated} drug classifications from TripSit data")
print(f"Not found in TripSit: {len(not_found)} drugs")

# Show some not found
if not_found:
    print("\nSample drugs not in TripSit:")
    for d in not_found[:20]:
        print(f"  - {d}")
