#!/usr/bin/env python3
"""
Comprehensive data collection script - run all sources
"""
import json
import time
import os
from pathlib import Path
from dataclasses import asdict
import sys

sys.path.insert(0, '/data/users/zacharie/whc/consciousness_study')

from src.data_collection.erowid_scraper import ErowidScraper
from src.data_collection.tripsit_client import TripSitAPI
from src.data_collection.psychonautwiki_client import PsychonautWikiAPI

BASE = Path('/data/users/zacharie/whc/consciousness_study/data/raw')


def scrape_erowid_comprehensive():
    """Scrape all key substances from Erowid"""
    print("\n" + "="*60)
    print("EROWID COMPREHENSIVE SCRAPE")
    print("="*60)
    
    scraper = ErowidScraper(cache_dir=str(BASE / 'erowid/cache'))
    reports_dir = BASE / 'erowid/reports'
    reports_dir.mkdir(parents=True, exist_ok=True)
    
    # All substances of interest
    substances = {
        # Classic Psychedelics
        'LSD': 'https://www.erowid.org/experiences/subs/exp_LSD.shtml',
        'Mushrooms': 'https://www.erowid.org/experiences/subs/exp_Mushrooms.shtml',
        'DMT': 'https://www.erowid.org/experiences/subs/exp_DMT.shtml',
        'Ayahuasca': 'https://www.erowid.org/experiences/subs/exp_Ayahuasca.shtml',
        'Mescaline': 'https://www.erowid.org/experiences/subs/exp_Mescaline.shtml',
        'Peyote': 'https://www.erowid.org/experiences/subs/exp_Peyote.shtml',
        '5-MeO-DMT': 'https://www.erowid.org/experiences/subs/exp_5-MeO-DMT.shtml',
        'San Pedro': 'https://www.erowid.org/experiences/subs/exp_Cacti_-_T._pachanoi.shtml',
        # Phenethylamines
        '2C-B': 'https://www.erowid.org/experiences/subs/exp_2C-B.shtml',
        '2C-E': 'https://www.erowid.org/experiences/subs/exp_2C-E.shtml',
        '2C-I': 'https://www.erowid.org/experiences/subs/exp_2C-I.shtml',
        '2C-T-7': 'https://www.erowid.org/experiences/subs/exp_2C-T-7.shtml',
        '2C-T-2': 'https://www.erowid.org/experiences/subs/exp_2C-T-2.shtml',
        'DOC': 'https://www.erowid.org/experiences/subs/exp_DOC.shtml',
        'DOB': 'https://www.erowid.org/experiences/subs/exp_DOB.shtml',
        # Tryptamines
        '4-AcO-DMT': 'https://www.erowid.org/experiences/subs/exp_4-AcO-DMT.shtml',
        'AMT': 'https://www.erowid.org/experiences/subs/exp_AMT.shtml',
        '5-MeO-DiPT': 'https://www.erowid.org/experiences/subs/exp_5-MeO-DiPT.shtml',
        'DiPT': 'https://www.erowid.org/experiences/subs/exp_DiPT.shtml',
        'DPT': 'https://www.erowid.org/experiences/subs/exp_DPT.shtml',
        # Lysergamides
        'LSA': 'https://www.erowid.org/experiences/subs/exp_Morning_Glory.shtml',
        # Empathogens
        'MDMA': 'https://www.erowid.org/experiences/subs/exp_MDMA.shtml',
        'MDA': 'https://www.erowid.org/experiences/subs/exp_MDA.shtml',
        # Dissociatives
        'Ketamine': 'https://www.erowid.org/experiences/subs/exp_Ketamine.shtml',
        'DXM': 'https://www.erowid.org/experiences/subs/exp_DXM.shtml',
        'PCP': 'https://www.erowid.org/experiences/subs/exp_PCP.shtml',
        'Nitrous': 'https://www.erowid.org/experiences/subs/exp_Nitrous_Oxide.shtml',
        'MXE': 'https://www.erowid.org/experiences/subs/exp_Methoxetamine.shtml',
        # Atypical
        'Salvia': 'https://www.erowid.org/experiences/subs/exp_Salvia_divinorum.shtml',
        'Ibogaine': 'https://www.erowid.org/experiences/subs/exp_Ibogaine.shtml',
        'Cannabis': 'https://www.erowid.org/experiences/subs/exp_Cannabis.shtml',
        'Amanita': 'https://www.erowid.org/experiences/subs/exp_Amanitas.shtml',
        # Deliriants
        'DPH': 'https://www.erowid.org/experiences/subs/exp_Diphenhydramine.shtml',
        'Datura': 'https://www.erowid.org/experiences/subs/exp_Datura.shtml',
        # Stimulants
        'Amphetamines': 'https://www.erowid.org/experiences/subs/exp_Amphetamines.shtml',
        'Cocaine': 'https://www.erowid.org/experiences/subs/exp_Cocaine.shtml',
        'Methamphetamine': 'https://www.erowid.org/experiences/subs/exp_Methamphetamine.shtml',
        # Depressants
        'GHB': 'https://www.erowid.org/experiences/subs/exp_GHB.shtml',
        'Alcohol': 'https://www.erowid.org/experiences/subs/exp_Alcohol.shtml',
        # Opioids
        'Opioids': 'https://www.erowid.org/experiences/subs/exp_Opioids.shtml',
        'Heroin': 'https://www.erowid.org/experiences/subs/exp_Heroin.shtml',
        'Kratom': 'https://www.erowid.org/experiences/subs/exp_Kratom.shtml',
        # MAOIs
        'Harmaline': 'https://www.erowid.org/experiences/subs/exp_Harmaline.shtml',
        'Syrian Rue': 'https://www.erowid.org/experiences/subs/exp_Syrian_Rue.shtml',
    }
    
    existing_ids = {f.stem for f in reports_dir.glob('*.json')}
    print(f"Already have {len(existing_ids)} reports cached")
    
    MAX_REPORTS = 100  # Per substance
    MAX_PAGES = 5
    
    report_index = {}
    
    for name, url in substances.items():
        print(f"\n{name}:")
        try:
            reports_list = scraper.get_experience_list_for_substance(url, max_pages=MAX_PAGES)
            report_index[name] = len(reports_list)
            print(f"  {len(reports_list)} indexed")
            
            scraped = 0
            for info in reports_list[:MAX_REPORTS]:
                if info['id'] in existing_ids:
                    continue
                try:
                    report = scraper.parse_experience_report(info['url'])
                    if report and len(report.text) > 100:
                        d = asdict(report)
                        with open(reports_dir / f'{report.id}.json', 'w') as f:
                            json.dump(d, f, indent=2)
                        existing_ids.add(report.id)
                        scraped += 1
                except:
                    pass
            print(f"  {scraped} new scraped")
        except Exception as e:
            print(f"  error: {e}")
    
    # Save index
    with open(BASE / 'erowid/report_index.json', 'w') as f:
        json.dump(report_index, f, indent=2)
    
    print(f"\nTotal reports: {len(existing_ids)}")
    return len(existing_ids)


def scrape_tripsit():
    """Scrape TripSit with careful rate limiting"""
    print("\n" + "="*60)
    print("TRIPSIT DATA COLLECTION")
    print("="*60)
    
    api = TripSitAPI()
    output_dir = BASE / 'tripsit'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    substances = [
        'LSD', 'DMT', 'Psilocybin', 'Mescaline', 'MDMA', 'Ketamine',
        'Cannabis', 'Salvia', 'DXM', 'Ibogaine', '2C-B', 'PCP',
        '5-MeO-DMT', 'Ayahuasca', 'DOC', 'MDA', 'MXE', '2C-E',
        'Cocaine', 'Amphetamine', 'GHB', 'Nitrous', 'DPH',
        '4-AcO-DMT', 'AL-LAD', '1P-LSD', 'Heroin', 'Morphine',
        'Alprazolam', 'Diazepam', 'GBL', 'Kratom', 'Caffeine',
        'Methamphetamine', 'Alcohol', 'Mushrooms', 'Oxycodone',
    ]
    
    # Load existing
    existing_file = output_dir / 'key_substances.json'
    if existing_file.exists():
        with open(existing_file) as f:
            collected = json.load(f)
    else:
        collected = []
    
    existing_names = {(s.get('name') or s.get('pretty_name', '')).lower() for s in collected}
    print(f"Already have {len(collected)} cached")
    
    for name in substances:
        if name.lower() in existing_names:
            continue
        try:
            print(f"  {name}...", end=' ', flush=True)
            data = api.get_drug(name)
            if data:
                collected.append(data)
                existing_names.add(name.lower())
                print(f"✓ ({len(data.get('combos', {}))} combos)")
            else:
                print("✗")
            time.sleep(5)  # Conservative rate limit
        except Exception as e:
            if '429' in str(e):
                print("rate limited - waiting 60s")
                time.sleep(60)
            else:
                print(f"error: {str(e)[:30]}")
            time.sleep(5)
    
    with open(output_dir / 'key_substances.json', 'w') as f:
        json.dump(collected, f, indent=2)
    
    print(f"\nTotal: {len(collected)} substances")
    return len(collected)


def print_status():
    """Print current collection status"""
    print("\n" + "="*60)
    print("CURRENT DATA STATUS")
    print("="*60)
    
    # PsychonautWiki
    pw = BASE / 'psychonautwiki/substances.json'
    if pw.exists():
        with open(pw) as f:
            data = json.load(f)
        print(f"\n📊 PsychonautWiki: {len(data)} substances")
    
    # TripSit
    ts = BASE / 'tripsit/key_substances.json'
    if ts.exists():
        with open(ts) as f:
            data = json.load(f)
        print(f"📊 TripSit: {len(data)} substances")
    
    # Erowid
    er = BASE / 'erowid/reports'
    if er.exists():
        reports = list(er.glob('*.json'))
        total = sum(len(json.load(open(f)).get('text', '')) for f in reports)
        print(f"📊 Erowid: {len(reports)} reports ({total//1000000}MB text)")


if __name__ == '__main__':
    print_status()
    
    # Run Erowid first (most important for consciousness research)
    erowid_count = scrape_erowid_comprehensive()
    
    # Then TripSit if time allows
    # tripsit_count = scrape_tripsit()
    
    print_status()
