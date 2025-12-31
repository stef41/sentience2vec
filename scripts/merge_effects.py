#!/usr/bin/env python3
"""
Merge PsychonautWiki effects into website data.js
"""

import json
import re

# Load the new effects
with open('/data/users/zacharie/whc/consciousness_study/website/data/effects_full.json', 'r') as f:
    new_data = json.load(f)

# Read the existing data.js
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'r') as f:
    data_js = f.read()

# Extract existing effects section for reference
# We'll replace the effects object entirely

# Build the new effects object 
new_effects = {}
effect_descriptions = {}

for key, effect in new_data['effects'].items():
    new_effects[key] = {
        "name": effect['name'],
        "receptors": effect['receptors'],
        "regions": effect['regions'],
        "category": effect['category'],
        "subcategory": effect.get('subcategory', ''),
        "description": effect.get('description', effect['name'])
    }
    effect_descriptions[key] = effect.get('description', effect['name'])

# Convert to JS format
effects_js = json.dumps(new_effects, indent=2)
descriptions_js = json.dumps(effect_descriptions, indent=2)

# Find and replace the effects section
# Pattern to match effects: { ... },
effects_pattern = r'(\s*effects:\s*)\{[^}]*(?:\{[^}]*\}[^}]*)*\}'

# Count braces to find the full effects object
def find_effects_section(js_content):
    match = re.search(r'(\s*effects:\s*)\{', js_content)
    if not match:
        return None, None
    
    start = match.start()
    brace_start = match.end() - 1
    
    # Count braces to find matching close
    depth = 0
    i = brace_start
    while i < len(js_content):
        if js_content[i] == '{':
            depth += 1
        elif js_content[i] == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    return None, None

# Also find and replace effectDescriptions
def find_effect_descriptions_section(js_content):
    match = re.search(r'(\s*effectDescriptions:\s*)\{', js_content)
    if not match:
        return None, None
    
    start = match.start()
    brace_start = match.end() - 1
    
    depth = 0
    i = brace_start
    while i < len(js_content):
        if js_content[i] == '{':
            depth += 1
        elif js_content[i] == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    return None, None


start, end = find_effects_section(data_js)
if start and end:
    print(f"Found effects section at {start}-{end}")
    # Replace
    new_data_js = data_js[:start] + f"\n    effects: {effects_js}" + data_js[end:]
    
    # Now update effectDescriptions
    start2, end2 = find_effect_descriptions_section(new_data_js)
    if start2 and end2:
        print(f"Found effectDescriptions section at {start2}-{end2}")
        new_data_js = new_data_js[:start2] + f"\n    effectDescriptions: {descriptions_js}" + new_data_js[end2:]
    
    # Write back
    with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'w') as f:
        f.write(new_data_js)
    
    print(f"Updated data.js with {len(new_effects)} effects!")
else:
    print("Could not find effects section!")
    
# Also create effectCategories for the dropdown
categories_js = json.dumps(new_data['effectCategories'], indent=2)
print(f"\nEffect categories: {list(new_data['effectCategories'].keys())}")
