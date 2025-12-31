#!/usr/bin/env python3
"""
Scrape effect -> drug relationships directly from PsychonautWiki effect pages.
This is more accurate than scraping drug pages, as it gets the actual listed substances.
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
                print(f"  Failed: {e}")
                return None

def get_all_effect_urls():
    """Get all effect URLs from the Subjective Effect Index."""
    print("Fetching Subjective Effect Index...")
    html = get_page(f"{BASE_URL}/wiki/Subjective_effect_index")
    if not html:
        return []
    
    soup = BeautifulSoup(html, 'html.parser')
    effect_urls = []
    
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
                                'Template:', 'index', 'Index', 'http', 'User:', 'Project:']
                if not any(s in href for s in skip_patterns):
                    full_url = href
                    if full_url not in [e['url'] for e in effect_urls]:
                        effect_urls.append({
                            'name': text,
                            'url': full_url
                        })
    
    return effect_urls

def extract_drugs_from_effect_page(effect_url, effect_name):
    """Extract the list of drugs from an effect page's 'Psychoactive substances' section."""
    html = get_page(BASE_URL + effect_url)
    if not html:
        return None
    
    soup = BeautifulSoup(html, 'html.parser')
    drugs = []
    
    # Find the "Psychoactive substances" heading
    substances_heading = None
    for heading in soup.find_all(['h2', 'h3', 'span']):
        text = heading.get_text().strip()
        if 'Psychoactive substance' in text:
            substances_heading = heading
            break
    
    if substances_heading:
        # Get the parent element and find the list after it
        parent = substances_heading.find_parent(['h2', 'h3', 'div'])
        if not parent:
            parent = substances_heading
        
        # Look for lists (ul) after the heading
        current = parent
        while current:
            current = current.find_next_sibling()
            if not current:
                break
            
            # Stop at next heading
            if current.name in ['h2', 'h3']:
                break
            
            # Extract drug links from lists
            if current.name == 'ul' or current.name == 'div':
                for link in current.find_all('a', href=True):
                    href = link.get('href', '')
                    drug_name = link.get_text().strip()
                    
                    # Filter for drug pages (skip effect pages, categories, etc.)
                    if href.startswith('/wiki/') and drug_name and len(drug_name) > 1:
                        skip = ['Category:', 'File:', 'Special:', 'effect', 'Effect',
                               'index', 'Index', '#', 'Template:', 'User:']
                        if not any(s in href for s in skip):
                            if drug_name not in drugs:
                                drugs.append(drug_name)
    
    return drugs

def main():
    print("=== PsychonautWiki Effect->Drug Scraper ===\n")
    
    # Get all effect URLs
    effect_urls = get_all_effect_urls()
    print(f"Found {len(effect_urls)} potential effect pages\n")
    
    # Filter to only actual effect pages (those with drugs listed)
    effect_drug_mapping = {}
    
    for i, effect in enumerate(effect_urls):
        effect_name = effect['name']
        effect_url = effect['url']
        
        # Skip obviously non-effect pages
        skip_names = ['psychonaut', 'wiki', 'substance', 'experience', 'index', 
                      'responsible', 'harm', 'guideline', 'tutorial', 'media']
        if any(s in effect_name.lower() for s in skip_names):
            continue
        
        print(f"[{i+1}/{len(effect_urls)}] {effect_name}...", end=' ', flush=True)
        
        drugs = extract_drugs_from_effect_page(effect_url, effect_name)
        
        if drugs and len(drugs) > 0:
            # Create normalized key
            key = effect_name.lower().replace(' ', '_').replace('-', '_')
            key = re.sub(r'[^a-z0-9_]', '', key)
            
            effect_drug_mapping[key] = {
                'name': effect_name,
                'drugs': drugs,
                'url': BASE_URL + effect_url,
                'source': 'PsychonautWiki'
            }
            print(f"{len(drugs)} drugs")
        else:
            print("no drugs listed")
        
        time.sleep(0.2)  # Be nice to server
    
    print(f"\n=== Results ===")
    print(f"Found {len(effect_drug_mapping)} effects with drug lists")
    print(f"Total drug-effect relationships: {sum(len(e['drugs']) for e in effect_drug_mapping.values())}")
    
    # Save the mapping
    output_file = OUTPUT_PATH / 'effect_drugs_direct.json'
    with open(output_file, 'w') as f:
        json.dump(effect_drug_mapping, f, indent=2)
    print(f"\nSaved to {output_file}")
    
    # Show some examples
    print("\n=== Sample Effect->Drugs Mappings ===")
    samples = ['visual_processing_deceleration', 'euphoria', 'time_distortion', 
               'ego_death', 'hallucination', 'sedation', 'stimulation']
    for sample in samples:
        for key, data in effect_drug_mapping.items():
            if sample in key:
                print(f"\n{data['name']}:")
                print(f"  Drugs: {', '.join(data['drugs'][:10])}" + 
                      (f"... (+{len(data['drugs'])-10})" if len(data['drugs']) > 10 else ""))
                break

if __name__ == '__main__':
    main()
