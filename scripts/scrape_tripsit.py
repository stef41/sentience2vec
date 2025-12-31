#!/usr/bin/env python3
"""
Fetch drug data from TripSit's factsheets API.
TripSit has a well-structured JSON API with drug information.
"""

import requests
import json

def fetch_tripsit_data():
    """Fetch all drug data from TripSit's API."""
    
    # TripSit factsheets API
    url = "https://tripbot.tripsit.me/api/tripsit/getAllDrugs"
    
    print("Fetching TripSit drug database...")
    response = requests.get(url, timeout=30)
    data = response.json()
    
    if 'data' not in data:
        print("Error: No data in response")
        return {}
    
    drugs = data['data'][0]  # All drugs are in first element
    
    print(f"Found {len(drugs)} drugs")
    
    # Process each drug
    processed = {}
    
    for drug_name, drug_data in drugs.items():
        processed_drug = {
            'name': drug_data.get('pretty_name', drug_name),
            'aliases': drug_data.get('aliases', []),
            'categories': drug_data.get('categories', []),
            'summary': drug_data.get('summary', ''),
            'effects': drug_data.get('effects', ''),
            'dose': drug_data.get('dose', {}),
            'duration': drug_data.get('duration', {}),
            'onset': drug_data.get('onset', {}),
            'properties': drug_data.get('properties', {}),
            'links': drug_data.get('links', {}),
            'sources': drug_data.get('sources', {})
        }
        
        # Extract class from categories
        if processed_drug['categories']:
            processed_drug['class'] = processed_drug['categories'][0] if processed_drug['categories'] else 'other'
        else:
            processed_drug['class'] = 'other'
        
        processed[drug_name] = processed_drug
    
    return processed

def main():
    drugs = fetch_tripsit_data()
    
    # Save to file
    output_path = '/data/users/zacharie/whc/consciousness_study/website/data/tripsit_data.json'
    with open(output_path, 'w') as f:
        json.dump(drugs, f, indent=2)
    
    print(f"Saved {len(drugs)} drugs to {output_path}")
    
    # Show category distribution
    categories = {}
    for name, data in drugs.items():
        for cat in data.get('categories', ['unknown']):
            categories[cat] = categories.get(cat, 0) + 1
    
    print("\nCategory distribution:")
    for cat, count in sorted(categories.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")
    
    # Show sample
    print("\nSample drugs:")
    for name, data in list(drugs.items())[:5]:
        print(f"  {data['name']}: {data['categories']} - {data['effects'][:100] if data['effects'] else 'N/A'}...")

if __name__ == "__main__":
    main()
