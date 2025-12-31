#!/usr/bin/env python3
"""
Merge PsychonautWiki effects with source URLs into website data.js
"""

import json
import re

# Load the scraped effects with URLs
with open('/data/users/zacharie/whc/consciousness_study/website/data/psychonautwiki_effects.json', 'r') as f:
    scraped_data = json.load(f)

# Read the existing data.js
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'r') as f:
    data_js = f.read()

# Build the new effects object with source URLs
new_effects = {}
effect_descriptions = {}

for key, effect in scraped_data['effects'].items():
    new_effects[key] = {
        "name": effect['name'],
        "receptors": effect['receptors'],
        "regions": effect['regions'],
        "category": effect['category'],
        "subcategory": effect.get('subcategory', ''),
        "description": effect.get('description', effect['name']),
        "source": effect.get('url', f"https://psychonautwiki.org/wiki/{effect['name'].replace(' ', '_')}"),
        "sourceLabel": "PsychonautWiki"
    }
    effect_descriptions[key] = effect.get('description', effect['name'])

# Convert to JS format
effects_js = json.dumps(new_effects, indent=2)
descriptions_js = json.dumps(effect_descriptions, indent=2)

# Find and replace the effects section
def find_section(js_content, section_name):
    pattern = rf'(\s*{section_name}:\s*)\{{'
    match = re.search(pattern, js_content)
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

# Replace effects section
start, end = find_section(data_js, 'effects')
if start and end:
    print(f"Found effects section at {start}-{end}")
    new_data_js = data_js[:start] + f"\n    effects: {effects_js}" + data_js[end:]
    
    # Replace effectDescriptions section
    start2, end2 = find_section(new_data_js, 'effectDescriptions')
    if start2 and end2:
        print(f"Found effectDescriptions section at {start2}-{end2}")
        new_data_js = new_data_js[:start2] + f"\n    effectDescriptions: {descriptions_js}" + new_data_js[end2:]
    
    # Write back
    with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'w') as f:
        f.write(new_data_js)
    
    print(f"Updated data.js with {len(new_effects)} effects including source URLs!")
else:
    print("Could not find effects section!")

# Verify one effect
sample = list(new_effects.items())[10]
print(f"\nSample effect: {sample[0]}")
print(json.dumps(sample[1], indent=2))
