#!/usr/bin/env python3
"""
Scrape drug information and their effects from PsychonautWiki.
"""

import requests
from bs4 import BeautifulSoup
import json
import time
import re

BASE_URL = "https://psychonautwiki.org"

# Get the list of substances from the substance index
def get_substance_list():
    """Get list of all substances from PsychonautWiki."""
    url = f"{BASE_URL}/wiki/Psychoactive_substance_index"
    print(f"Fetching substance index...")
    
    response = requests.get(url, timeout=30)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    substances = []
    content = soup.find('div', {'class': 'mw-parser-output'})
    
    if content:
        for link in content.find_all('a'):
            href = link.get('href', '')
            if href.startswith('/wiki/') and ':' not in href:
                name = link.get_text().strip()
                if name and len(name) > 1:
                    # Skip meta pages
                    if any(skip in name.lower() for skip in ['index', 'category', 'template', 'help', 'special']):
                        continue
                    substances.append({
                        'name': name,
                        'url': f"{BASE_URL}{href}"
                    })
    
    # Remove duplicates
    seen = set()
    unique = []
    for s in substances:
        if s['name'] not in seen:
            seen.add(s['name'])
            unique.append(s)
    
    return unique

def scrape_substance_page(url, name):
    """Scrape a single substance page for effects and class."""
    try:
        response = requests.get(url, timeout=15)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        data = {
            'name': name,
            'url': url,
            'class': None,
            'effects': [],
            'pharmacology': {}
        }
        
        content = soup.find('div', {'class': 'mw-parser-output'})
        if not content:
            return data
        
        # Find drug class from infobox or text
        infobox = soup.find('table', {'class': 'infobox'})
        if infobox:
            for row in infobox.find_all('tr'):
                header = row.find('th')
                cell = row.find('td')
                if header and cell:
                    header_text = header.get_text().strip().lower()
                    cell_text = cell.get_text().strip()
                    
                    if 'class' in header_text or 'chemical class' in header_text:
                        data['class'] = cell_text
                    elif 'pharmacology' in header_text:
                        data['pharmacology']['mechanism'] = cell_text
        
        # Find effects section
        effects_section = None
        for header in content.find_all(['h2', 'h3']):
            if 'subjective effects' in header.get_text().lower():
                effects_section = header
                break
        
        if effects_section:
            # Get all links in the effects section
            current = effects_section.find_next_sibling()
            while current and current.name not in ['h2']:
                if current.name in ['ul', 'div', 'p']:
                    for link in current.find_all('a'):
                        href = link.get('href', '')
                        effect_name = link.get_text().strip()
                        if effect_name and '/wiki/' in href:
                            # Clean effect name
                            effect_key = re.sub(r'[^a-z0-9]+', '_', effect_name.lower()).strip('_')
                            if effect_key and len(effect_key) > 2:
                                data['effects'].append({
                                    'name': effect_name,
                                    'key': effect_key
                                })
                current = current.find_next_sibling() if current else None
        
        # Deduplicate effects
        seen_effects = set()
        unique_effects = []
        for e in data['effects']:
            if e['key'] not in seen_effects:
                seen_effects.add(e['key'])
                unique_effects.append(e)
        data['effects'] = unique_effects
        
        return data
        
    except Exception as e:
        print(f"  Error scraping {name}: {e}")
        return {'name': name, 'url': url, 'effects': [], 'class': None}

def main():
    # Get substance list
    substances = get_substance_list()
    print(f"Found {len(substances)} substances")
    
    # Limit for testing - remove for full scrape
    # substances = substances[:50]
    
    all_data = {}
    
    for i, substance in enumerate(substances):
        print(f"[{i+1}/{len(substances)}] Scraping {substance['name']}...")
        
        data = scrape_substance_page(substance['url'], substance['name'])
        
        if data['effects'] or data['class']:
            all_data[substance['name']] = data
            print(f"  Found {len(data['effects'])} effects, class: {data['class']}")
        
        # Rate limit
        time.sleep(0.3)
    
    # Save results
    output_path = '/data/users/zacharie/whc/consciousness_study/website/data/psychonautwiki_drugs.json'
    with open(output_path, 'w') as f:
        json.dump(all_data, f, indent=2)
    
    print(f"\nSaved {len(all_data)} substances to {output_path}")
    
    # Summary
    classes = {}
    for name, data in all_data.items():
        cls = data.get('class') or 'unknown'
        classes[cls] = classes.get(cls, 0) + 1
    
    print("\nClass distribution:")
    for cls, count in sorted(classes.items(), key=lambda x: -x[1]):
        print(f"  {cls}: {count}")

if __name__ == "__main__":
    main()
