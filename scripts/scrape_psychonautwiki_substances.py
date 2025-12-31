#!/usr/bin/env python3
"""
Scrape drug -> effect relationships from PsychonautWiki substance pages.
Uses cached category pages to get substance lists, then scrapes individual pages.
"""

import requests
from bs4 import BeautifulSoup
import json
import time
import re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

OUTPUT_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://psychonautwiki.org"
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (compatible; research project)',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
}

# Categories that contain substances
SUBSTANCE_CATEGORIES = [
    '/wiki/Psychedelics',
    '/wiki/Dissociatives', 
    '/wiki/Entactogens',
    '/wiki/Stimulants',
    '/wiki/Depressants',
    '/wiki/Cannabinoids',
    '/wiki/Opioids',
    '/wiki/Deliriants',
    '/wiki/Nootropics'
]

def get_page(url, retries=3):
    """Fetch a page with retries."""
    for i in range(retries):
        try:
            response = requests.get(url, headers=HEADERS, timeout=30)
            response.raise_for_status()
            return response.text
        except Exception as e:
            if i < retries - 1:
                time.sleep(1)
            else:
                print(f"  Failed to fetch {url}: {e}")
                return None

def extract_substances_from_category(category_url):
    """Extract substance names and links from a category page."""
    html = get_page(BASE_URL + category_url)
    if not html:
        return []
    
    soup = BeautifulSoup(html, 'html.parser')
    substances = []
    
    # Look for substance links in the main content
    content = soup.find('div', {'id': 'mw-content-text'})
    if content:
        # Find all internal wiki links
        for link in content.find_all('a', href=True):
            href = link.get('href', '')
            # Filter for substance pages (they typically have /wiki/ prefix and are capitalized)
            if href.startswith('/wiki/') and not any(x in href.lower() for x in 
                ['category:', 'file:', 'special:', 'help:', 'talk:', 'user:', 
                 'template:', 'subjective', 'effect', 'index', 'list']):
                name = link.get_text().strip()
                if name and len(name) > 1 and not name.startswith('['):
                    substances.append({
                        'name': name,
                        'url': href
                    })
    
    return substances

def extract_effects_from_substance(substance_url, substance_name):
    """Extract subjective effects from a substance page."""
    html = get_page(BASE_URL + substance_url)
    if not html:
        return None
    
    soup = BeautifulSoup(html, 'html.parser')
    
    effects = []
    drug_class = None
    
    # Try to find drug class
    class_row = soup.find(lambda tag: tag.name == 'th' and 'Class' in tag.get_text())
    if class_row:
        class_cell = class_row.find_next_sibling('td')
        if class_cell:
            drug_class = class_cell.get_text().strip()
    
    # Find subjective effects section
    effects_heading = soup.find(lambda tag: tag.name in ['h2', 'h3'] and 
                                'Subjective effects' in tag.get_text())
    
    if effects_heading:
        # Get all effect links in the effects section
        current = effects_heading.find_next_sibling()
        while current and current.name not in ['h2', 'h3']:
            if current.name in ['ul', 'div', 'p']:
                for link in current.find_all('a', href=True):
                    href = link.get('href', '')
                    if '/wiki/' in href and 'effect' not in href.lower() and 'subjective' not in href.lower():
                        # This might be an effect link
                        effect_name = link.get_text().strip()
                        if effect_name and len(effect_name) > 1:
                            effects.append(effect_name)
            current = current.find_next_sibling()
    
    # Also look for effects in info boxes or tables
    infobox = soup.find('table', class_='infobox')
    if infobox:
        effects_row = infobox.find(lambda tag: tag.name == 'th' and 'effect' in tag.get_text().lower())
        if effects_row:
            effects_cell = effects_row.find_next_sibling('td')
            if effects_cell:
                for link in effects_cell.find_all('a'):
                    effect_name = link.get_text().strip()
                    if effect_name and effect_name not in effects:
                        effects.append(effect_name)
    
    return {
        'name': substance_name,
        'url': substance_url,
        'class': drug_class,
        'effects': list(set(effects)),  # Remove duplicates
        'source': 'PsychonautWiki'
    }

def main():
    print("=== PsychonautWiki Drug-Effect Scraper ===\n")
    
    # Collect all substance URLs from category pages
    all_substances = {}
    
    print("Step 1: Collecting substances from category pages...")
    for category in SUBSTANCE_CATEGORIES:
        print(f"  Fetching {category}...")
        substances = extract_substances_from_category(category)
        for sub in substances:
            if sub['name'] not in all_substances:
                all_substances[sub['name']] = sub
        print(f"    Found {len(substances)} substances")
        time.sleep(0.5)
    
    print(f"\nTotal unique substances found: {len(all_substances)}")
    
    # Scrape individual substance pages for effects
    print("\nStep 2: Scraping substance pages for effects...")
    
    drug_effects = {}
    effect_drugs = {}
    
    # Process substances (limit to first 100 for speed)
    substances_list = list(all_substances.values())[:150]
    
    for i, sub in enumerate(substances_list):
        print(f"  [{i+1}/{len(substances_list)}] {sub['name']}...", end=' ')
        
        data = extract_effects_from_substance(sub['url'], sub['name'])
        if data and data.get('effects'):
            drug_effects[sub['name']] = data
            
            # Build reverse mapping
            for effect in data['effects']:
                if effect not in effect_drugs:
                    effect_drugs[effect] = []
                if sub['name'] not in effect_drugs[effect]:
                    effect_drugs[effect].append(sub['name'])
            
            print(f"{len(data['effects'])} effects")
        else:
            print("no effects found")
        
        time.sleep(0.3)  # Be nice to the server
    
    print(f"\nProcessed {len(drug_effects)} drugs with effects")
    print(f"Found {len(effect_drugs)} unique effects")
    
    # Save results
    output_file = OUTPUT_PATH / 'psychonautwiki_drug_effects.json'
    with open(output_file, 'w') as f:
        json.dump(drug_effects, f, indent=2)
    print(f"\nSaved drug effects to {output_file}")
    
    output_file2 = OUTPUT_PATH / 'psychonautwiki_effect_drugs.json'
    with open(output_file2, 'w') as f:
        json.dump(effect_drugs, f, indent=2)
    print(f"Saved effect-drug mapping to {output_file2}")
    
    # Stats
    print("\n=== Top Effects by Drug Count ===")
    sorted_effects = sorted(effect_drugs.items(), key=lambda x: len(x[1]), reverse=True)
    for effect, drugs in sorted_effects[:20]:
        print(f"  {effect}: {len(drugs)} drugs")

if __name__ == '__main__':
    main()
