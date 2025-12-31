#!/usr/bin/env python3
"""
Comprehensive substance scraper - extracts ALL substances from:
1. PsychonautWiki Substance Index (300+ substances)
2. TripSit Wiki Categories (150+ ethnobotanicals, 69 main drugs)
3. TripSit Factsheets API (551 entries)
"""

import json
import time
import re
import requests
from pathlib import Path
from bs4 import BeautifulSoup
from urllib.parse import urljoin

# Base paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "raw"

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (compatible; ConsciousnessResearchBot/1.0; academic research)'
}

def scrape_psychonautwiki_index():
    """Scrape ALL substances from PsychonautWiki substance index page"""
    print("\n" + "="*70)
    print("PSYCHONAUTWIKI - SCRAPING FULL INDEX")
    print("="*70)
    
    url = "https://psychonautwiki.org/wiki/Psychoactive_substance_index"
    
    try:
        response = requests.get(url, headers=HEADERS, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        substances = set()
        
        # Find all links to substance pages
        for link in soup.find_all('a', href=True):
            href = link['href']
            text = link.get_text(strip=True)
            
            # Filter for substance pages (wiki internal links)
            if href.startswith('/wiki/') and not any(x in href for x in [
                'Category:', 'Special:', 'Template:', 'File:', 'Talk:', 
                'Help:', 'Project:', 'index', 'Index', 'User:', 'Main_Page'
            ]):
                # Extract substance name
                name = href.replace('/wiki/', '').replace('_', ' ')
                # Skip meta pages
                if name and len(name) > 1 and not name.startswith('#'):
                    # Clean up names with anchors
                    if '#' in name:
                        name = name.split('#')[0]
                    if name:
                        substances.add(name)
        
        # Also look for specific substance patterns in the HTML
        content = soup.find('div', {'class': 'mw-parser-output'})
        if content:
            for link in content.find_all('a'):
                href = link.get('href', '')
                if '/wiki/' in href:
                    name = href.split('/wiki/')[-1].replace('_', ' ')
                    if '#' in name:
                        name = name.split('#')[0]
                    if name and len(name) > 1:
                        substances.add(name)
        
        print(f"Found {len(substances)} unique substance names")
        return sorted(list(substances))
        
    except Exception as e:
        print(f"Error scraping index: {e}")
        return []


def scrape_tripsit_wiki_category(category_url):
    """Scrape a TripSit wiki category page"""
    try:
        response = requests.get(category_url, headers=HEADERS, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        substances = []
        
        # Find category members
        category_div = soup.find('div', {'class': 'mw-category'})
        if category_div:
            for link in category_div.find_all('a'):
                name = link.get_text(strip=True)
                if name and not name.startswith('Category:'):
                    substances.append(name)
        
        # Also check subcategories
        subcats = []
        subcat_div = soup.find('div', {'class': 'mw-subcategories'})
        if subcat_div:
            for link in subcat_div.find_all('a'):
                href = link.get('href', '')
                if 'Category:' in href:
                    subcats.append(urljoin(category_url, href))
        
        return substances, subcats
        
    except Exception as e:
        print(f"  Error: {e}")
        return [], []


def scrape_tripsit_wiki_all():
    """Scrape ALL TripSit wiki categories"""
    print("\n" + "="*70)
    print("TRIPSIT WIKI - SCRAPING ALL CATEGORIES")
    print("="*70)
    
    base_categories = [
        "https://wiki.tripsit.me/wiki/Category:Drugs",
        "https://wiki.tripsit.me/wiki/Category:Ethnobotanical",
        "https://wiki.tripsit.me/wiki/Category:Psychedelic",
        "https://wiki.tripsit.me/wiki/Category:Dissociative",
        "https://wiki.tripsit.me/wiki/Category:Stimulant",
        "https://wiki.tripsit.me/wiki/Category:Depressant",
        "https://wiki.tripsit.me/wiki/Category:Opioid",
        "https://wiki.tripsit.me/wiki/Category:Benzodiazepine",
        "https://wiki.tripsit.me/wiki/Category:Deliriant",
        "https://wiki.tripsit.me/wiki/Category:Nootropic",
        "https://wiki.tripsit.me/wiki/Category:Research_Chemical",
    ]
    
    all_substances = set()
    visited = set()
    to_visit = list(base_categories)
    
    while to_visit:
        url = to_visit.pop(0)
        if url in visited:
            continue
        visited.add(url)
        
        cat_name = url.split('Category:')[-1] if 'Category:' in url else url
        print(f"  {cat_name}...", end=' ', flush=True)
        
        substances, subcats = scrape_tripsit_wiki_category(url)
        all_substances.update(substances)
        print(f"{len(substances)} items")
        
        for subcat in subcats:
            if subcat not in visited:
                to_visit.append(subcat)
        
        time.sleep(1)  # Rate limiting
    
    print(f"\nTotal unique substances from wiki: {len(all_substances)}")
    return sorted(list(all_substances))


def scrape_tripsit_factsheets():
    """Scrape TripSit factsheets API (551 drugs)"""
    print("\n" + "="*70)
    print("TRIPSIT FACTSHEETS API")
    print("="*70)
    
    # The factsheets use the tripsit API
    api_url = "https://tripbot.tripsit.me/api/tripsit/getAllDrugs"
    
    try:
        response = requests.get(api_url, headers=HEADERS, timeout=60)
        response.raise_for_status()
        data = response.json()
        
        if 'data' in data and isinstance(data['data'], dict):
            drugs = list(data['data'].keys())
            print(f"Found {len(drugs)} drugs in factsheets API")
            
            # Save full data
            output_dir = DATA_DIR / "tripsit"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            with open(output_dir / "all_drugs_api.json", 'w') as f:
                json.dump(data['data'], f, indent=2, default=str)
            
            return drugs, data['data']
        else:
            print(f"Unexpected response format")
            return [], {}
            
    except Exception as e:
        print(f"Error: {e}")
        return [], {}


def fetch_psychonautwiki_substance(name, api):
    """Fetch detailed substance data from PsychonautWiki GraphQL API"""
    try:
        substance = api.get_substance(name)
        if substance:
            return substance
    except Exception as e:
        pass
    return None


def main():
    print("="*70)
    print("COMPREHENSIVE SUBSTANCE SCRAPER")
    print("Target: ALL substances from PsychonautWiki + TripSit")
    print("="*70)
    
    # 1. Get TripSit factsheets first (this is the richest source)
    tripsit_names, tripsit_full_data = scrape_tripsit_factsheets()
    
    # 2. Get TripSit wiki categories
    wiki_substances = scrape_tripsit_wiki_all()
    
    # 3. Get PsychonautWiki index
    pw_substances = scrape_psychonautwiki_index()
    
    # Combine all substance names
    all_names = set()
    all_names.update(tripsit_names)
    all_names.update(wiki_substances)
    all_names.update(pw_substances)
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"TripSit Factsheets API: {len(tripsit_names)} drugs")
    print(f"TripSit Wiki: {len(wiki_substances)} entries")
    print(f"PsychonautWiki Index: {len(pw_substances)} entries")
    print(f"TOTAL UNIQUE: {len(all_names)} substance names")
    
    # Save master list
    output_dir = DATA_DIR / "master"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(output_dir / "all_substance_names.json", 'w') as f:
        json.dump({
            'tripsit_factsheets': sorted(tripsit_names),
            'tripsit_wiki': sorted(wiki_substances),
            'psychonautwiki_index': sorted(pw_substances),
            'all_unique': sorted(list(all_names)),
            'counts': {
                'tripsit_factsheets': len(tripsit_names),
                'tripsit_wiki': len(wiki_substances),
                'psychonautwiki_index': len(pw_substances),
                'total_unique': len(all_names)
            }
        }, f, indent=2)
    
    print(f"\nSaved master list to {output_dir / 'all_substance_names.json'}")
    
    # Save TripSit full data separately
    tripsit_dir = DATA_DIR / "tripsit"
    tripsit_dir.mkdir(parents=True, exist_ok=True)
    
    with open(tripsit_dir / "factsheets_full.json", 'w') as f:
        json.dump(tripsit_full_data, f, indent=2, default=str)
    
    print(f"Saved TripSit full data to {tripsit_dir / 'factsheets_full.json'}")
    
    return all_names, tripsit_full_data


if __name__ == "__main__":
    main()
