#!/usr/bin/env python3
"""
Scrape specific drug pages from PsychonautWiki for their subjective effects.
Uses known drug names rather than crawling category pages.
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

# Known drug names to scrape (these have dedicated pages with effects)
DRUGS = [
    # Psychedelics
    "LSD", "Psilocybin_mushrooms", "DMT", "Mescaline", "2C-B", "2C-E", "2C-I", "2C-C", "2C-P",
    "4-AcO-DMT", "4-AcO-MET", "4-HO-MET", "4-HO-MiPT", "5-MeO-DMT", "5-MeO-MiPT",
    "AL-LAD", "1P-LSD", "ETH-LAD", "LSA", "LSZ", "ALD-52",
    "DOM", "DOC", "DOI", "DOB", "TMA-2", "TMA-6",
    "25I-NBOMe", "25C-NBOMe", "25B-NBOMe",
    "Ayahuasca", "Ibogaine", "Salvia_divinorum", "DPT", "MET", "DiPT", "DIPT",
    
    # Dissociatives
    "Ketamine", "PCP", "DXM", "Methoxetamine", "3-MeO-PCP", "3-MeO-PCE",
    "Nitrous", "Diphenidine", "Ephenidine", "O-PCE", "3-HO-PCP",
    
    # Entactogens/Empathogens
    "MDMA", "MDA", "MDEA", "6-APB", "5-APB", "5-MAPB", "6-APDB",
    "Mephedrone", "Methylone", "Butylone", "4-FA", "4-FMA", "Bk-MDMA",
    
    # Stimulants
    "Amphetamine", "Methamphetamine", "Cocaine", "Caffeine", "Nicotine",
    "Methylphenidate", "Modafinil", "Adrafinil", "MDPV", "A-PVP", "A-PHP",
    "2-FMA", "3-FPM", "NEP", "Hexen", "4F-MPH", "Ephedrine", "Pseudoephedrine",
    
    # Depressants
    "Alcohol", "GHB", "GBL", "Phenibut", "Pregabalin", "Gabapentin",
    "Alprazolam", "Diazepam", "Clonazepam", "Lorazepam", "Etizolam",
    "Clonazolam", "Flualprazolam", "Barbiturates", "Methaqualone",
    "Zopiclone", "Zolpidem",
    
    # Opioids
    "Heroin", "Morphine", "Codeine", "Oxycodone", "Hydrocodone", "Fentanyl",
    "Tramadol", "Kratom", "O-Desmethyltramadol", "Buprenorphine", "Methadone",
    "U-47700", "Tianeptine",
    
    # Cannabinoids
    "Cannabis", "THC", "CBD", "JWH-018", "Synthetic_cannabinoids",
    
    # Deliriants
    "Diphenhydramine", "Datura", "DPH", "Scopolamine", "Benzydamine",
    
    # Nootropics
    "Piracetam", "Aniracetam", "Phenylpiracetam", "Noopept", "Modafinil"
]

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

def extract_effects(html, drug_name):
    """Extract subjective effects from a PsychonautWiki substance page."""
    soup = BeautifulSoup(html, 'html.parser')
    
    effects = set()
    drug_class = None
    
    # Find the Subjective effects section
    for heading in soup.find_all(['h2', 'h3', 'span']):
        text = heading.get_text()
        if 'Subjective effect' in text or 'subjective effect' in text:
            # Found effects section - look for all links in subsequent content
            parent = heading.find_parent(['h2', 'h3', 'div'])
            if not parent:
                parent = heading
            
            # Get next siblings until next h2
            current = parent
            while current:
                current = current.find_next_sibling()
                if not current:
                    break
                if current.name == 'h2':
                    break
                    
                # Extract effect links
                for link in current.find_all('a', href=True):
                    href = link.get('href', '')
                    text = link.get_text().strip()
                    
                    # Check if it's an effect link (contains effect-related keywords or links to effect pages)
                    if text and len(text) > 2 and len(text) < 60:
                        # Filter out navigation and non-effect links
                        if not any(x in text.lower() for x in 
                                   ['edit', 'source', 'wiki', 'cite', 'reference', 
                                    'http', 'category', '↑', '[', ']']):
                            effects.add(text)
    
    # Try to find drug class from infobox
    infobox = soup.find('table', class_='infobox')
    if infobox:
        for row in infobox.find_all('tr'):
            th = row.find('th')
            td = row.find('td')
            if th and td:
                header = th.get_text().strip().lower()
                if 'class' in header or 'type' in header:
                    drug_class = td.get_text().strip()
                    break
    
    return {
        'name': drug_name.replace('_', ' '),
        'effects': list(effects),
        'class': drug_class,
        'source': 'PsychonautWiki',
        'url': f'{BASE_URL}/wiki/{drug_name}'
    }

def main():
    print("=== PsychonautWiki Targeted Drug Scraper ===\n")
    
    drug_effects = {}
    effect_drugs = {}
    
    for i, drug in enumerate(DRUGS):
        url = f"{BASE_URL}/wiki/{drug}"
        print(f"[{i+1}/{len(DRUGS)}] {drug}...", end=' ', flush=True)
        
        html = get_page(url)
        if html:
            data = extract_effects(html, drug)
            if data['effects']:
                drug_effects[data['name']] = data
                
                # Build reverse mapping
                for effect in data['effects']:
                    if effect not in effect_drugs:
                        effect_drugs[effect] = []
                    if data['name'] not in effect_drugs[effect]:
                        effect_drugs[effect].append(data['name'])
                
                print(f"{len(data['effects'])} effects")
            else:
                print("no effects")
        else:
            print("page not found")
        
        time.sleep(0.2)
    
    print(f"\n=== Results ===")
    print(f"Scraped {len(drug_effects)} drugs with effects")
    print(f"Found {len(effect_drugs)} unique effects")
    
    # Save
    with open(OUTPUT_PATH / 'psychonautwiki_drug_effects.json', 'w') as f:
        json.dump(drug_effects, f, indent=2)
    print(f"\nSaved to psychonautwiki_drug_effects.json")
    
    with open(OUTPUT_PATH / 'psychonautwiki_effect_drugs.json', 'w') as f:
        json.dump(effect_drugs, f, indent=2)
    print(f"Saved to psychonautwiki_effect_drugs.json")
    
    # Top effects
    print("\n=== Top 20 Effects ===")
    sorted_effects = sorted(effect_drugs.items(), key=lambda x: len(x[1]), reverse=True)
    for effect, drugs in sorted_effects[:20]:
        print(f"  {effect}: {', '.join(drugs[:5])}" + ("..." if len(drugs) > 5 else ""))

if __name__ == '__main__':
    main()
