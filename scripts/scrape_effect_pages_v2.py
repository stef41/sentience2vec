#!/usr/bin/env python3
"""
Scrape effect -> drug relationships directly from PsychonautWiki effect pages.
Uses the actual drug lists from each effect page.
"""

import requests
from bs4 import BeautifulSoup
import json
import time
import re
from pathlib import Path

OUTPUT_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://psychonautwiki.org"
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (compatible; research project)',
    'Accept': 'text/html'
}

def get_page(url, retries=3):
    """Fetch a page with retries."""
    for i in range(retries):
        try:
            response = requests.get(url, headers=HEADERS, timeout=30)
            if response.status_code == 404:
                return None
            response.raise_for_status()
            return response.text
        except Exception as e:
            if i < retries - 1:
                time.sleep(1)
            else:
                return None

def get_all_effect_urls():
    """Get all effect URLs from the Subjective Effect Index."""
    print("Fetching Subjective Effect Index...")
    html = get_page(f"{BASE_URL}/wiki/Subjective_effect_index")
    if not html:
        return []
    
    soup = BeautifulSoup(html, 'html.parser')
    effect_urls = []
    seen = set()
    
    # Find all links in the main content
    content = soup.find('div', {'id': 'mw-content-text'})
    if content:
        for link in content.find_all('a', href=True):
            href = link.get('href', '')
            text = link.get_text().strip()
            
            # Filter for effect pages
            if href.startswith('/wiki/') and text and len(text) > 2:
                # Skip non-effect pages
                skip_patterns = ['Category:', 'File:', 'Special:', 'Help:', 'Talk:', 
                                'Template:', 'User:', 'Project:', '#']
                if not any(s in href for s in skip_patterns):
                    if href not in seen:
                        seen.add(href)
                        effect_urls.append({
                            'name': text,
                            'url': href
                        })
    
    return effect_urls

def extract_drugs_from_effect_page(effect_url):
    """Extract the list of drugs from an effect page's 'Psychoactive substances' section."""
    html = get_page(BASE_URL + effect_url)
    if not html:
        return []
    
    soup = BeautifulSoup(html, 'html.parser')
    drugs = []
    
    # Find the SMW format list (this is what PsychonautWiki uses for drug lists)
    smw_list = soup.find('ul', class_='smw-format')
    if smw_list:
        for item in smw_list.find_all('li', class_='smw-row'):
            link = item.find('a')
            if link:
                drug_name = link.get_text().strip()
                if drug_name and drug_name not in drugs:
                    drugs.append(drug_name)
    
    return drugs

def main():
    print("=== PsychonautWiki Effect->Drug Scraper ===\n")
    
    # Get all effect URLs
    effect_urls = get_all_effect_urls()
    print(f"Found {len(effect_urls)} potential effect pages\n")
    
    effect_drug_mapping = {}
    effects_with_drugs = 0
    
    for i, effect in enumerate(effect_urls):
        effect_name = effect['name']
        effect_url = effect['url']
        
        # Skip category/index pages
        skip_names = ['index', 'Index', 'psychonaut', 'wiki', 'guideline', 'tutorial']
        if any(s.lower() in effect_name.lower() for s in skip_names):
            continue
        
        print(f"[{i+1}/{len(effect_urls)}] {effect_name}...", end=' ', flush=True)
        
        drugs = extract_drugs_from_effect_page(effect_url)
        
        if drugs:
            # Create normalized key
            key = effect_name.lower().replace(' ', '_').replace('-', '_')
            key = re.sub(r'[^a-z0-9_]', '', key)
            
            effect_drug_mapping[key] = {
                'name': effect_name,
                'drugs': drugs,
                'url': BASE_URL + effect_url
            }
            effects_with_drugs += 1
            print(f"{len(drugs)} drugs")
        else:
            print("-")
        
        time.sleep(0.15)
    
    print(f"\n=== Results ===")
    print(f"Effects with drug lists: {effects_with_drugs}")
    print(f"Total relationships: {sum(len(e['drugs']) for e in effect_drug_mapping.values())}")
    
    # Save
    output_file = OUTPUT_PATH / 'effect_drugs_direct.json'
    with open(output_file, 'w') as f:
        json.dump(effect_drug_mapping, f, indent=2)
    print(f"\nSaved to {output_file}")
    
    # Show examples
    print("\n=== Sample Mappings ===")
    for key in ['visual_processing_deceleration', 'euphoria', 'time_distortion', 'ego_death', 'sedation']:
        if key in effect_drug_mapping:
            data = effect_drug_mapping[key]
            print(f"\n{data['name']}: {', '.join(data['drugs'][:8])}" + 
                  (f"... (+{len(data['drugs'])-8})" if len(data['drugs']) > 8 else ""))

if __name__ == '__main__':
    main()
