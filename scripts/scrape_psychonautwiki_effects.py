#!/usr/bin/env python3
"""
Scrape the full PsychonautWiki Subjective Effect Index
and map effects to receptors and brain regions.
"""

import requests
from bs4 import BeautifulSoup
import json
import re
import time

BASE_URL = "https://psychonautwiki.org"

# Main categories with their subcategories
EFFECT_CATEGORIES = {
    "Visual": {
        "url": "/wiki/Visual_effects",
        "subcategories": ["Amplifications", "Suppressions", "Distortions", "Geometry", "Hallucinatory states"]
    },
    "Auditory": {
        "url": "/wiki/Auditory_effects",
        "subcategories": []
    },
    "Tactile": {
        "url": "/wiki/Tactile_effects",
        "subcategories": []
    },
    "Disconnective": {
        "url": "/wiki/Disconnective_effects",
        "subcategories": []
    },
    "Smell & Taste": {
        "url": "/wiki/Smell_and_taste_effects",
        "subcategories": []
    },
    "Multisensory": {
        "url": "/wiki/Multisensory_effects",
        "subcategories": []
    },
    "Cognitive": {
        "url": "/wiki/Cognitive_effects",
        "subcategories": ["Enhancements", "Depressions", "Intensifications", "Suppressions", "Novel", "Psychological", "Transpersonal"]
    },
    "Physical": {
        "url": "/wiki/Physical_effects",
        "subcategories": ["Amplifications", "Suppressions", "Alterations"]
    },
    "Uncomfortable Physical": {
        "url": "/wiki/Uncomfortable_physical_effects",
        "subcategories": ["Cardiovascular", "Cerebrovascular", "Bodily"]
    }
}

# Known receptor associations for effect types
RECEPTOR_MAPPINGS = {
    # Visual/Hallucinatory effects - primarily 5-HT2A mediated
    "geometry": ["5-HT2A"],
    "visual": ["5-HT2A", "5-HT2C"],
    "hallucination": ["5-HT2A", "5-HT2C", "D2"],
    "color": ["5-HT2A"],
    "drifting": ["5-HT2A"],
    "morphing": ["5-HT2A"],
    "pattern": ["5-HT2A"],
    "tracers": ["5-HT2A"],
    
    # Cognitive enhancements
    "creativity": ["5-HT2A", "D2"],
    "analysis": ["5-HT2A", "D1"],
    "memory": ["M1", "D1", "NMDA"],
    "focus": ["D1", "D2", "alpha2A"],
    "thought": ["5-HT2A", "D2"],
    
    # Emotional/Mood
    "euphoria": ["MOR", "D2", "5-HT2A", "CB1"],
    "dysphoria": ["KOR", "D2"],
    "empathy": ["5-HT2A", "SERT"],
    "anxiety": ["5-HT1A", "GABA-A", "5-HT2C"],
    "depression": ["SERT", "NET", "5-HT1A"],
    
    # Physical effects
    "stimulation": ["D1", "D2", "NET", "DAT", "alpha2A"],
    "sedation": ["GABA-A", "H1", "alpha2A", "MOR"],
    "pain": ["MOR", "KOR", "CB1"],
    "nausea": ["5-HT3", "D2", "MOR"],
    "appetite": ["CB1", "5-HT2C", "H1"],
    
    # Dissociative effects
    "dissociation": ["NMDA", "KOR", "Sigma1"],
    "derealization": ["NMDA", "5-HT2A"],
    "depersonalization": ["NMDA", "5-HT2A"],
    
    # Transpersonal
    "ego": ["5-HT2A", "NMDA"],
    "mystical": ["5-HT2A"],
    "unity": ["5-HT2A"],
    "spiritual": ["5-HT2A"],
    
    # Physical alterations
    "pupil": ["alpha1A", "5-HT2A", "MOR"],
    "heart": ["beta1", "alpha1A", "M2"],
    "blood pressure": ["alpha1A", "alpha2A"],
    "temperature": ["5-HT2A", "D2"],
    "muscle": ["GABA-A", "alpha2A"],
    
    # Sleep-related
    "dream": ["5-HT2A", "M1"],
    "sleep": ["H1", "GABA-A", "melatonin"],
    "wakefulness": ["D1", "D2", "H3", "OX1"],
}

# Brain region associations
REGION_MAPPINGS = {
    "visual": ["visual_cortex", "V1", "V2", "V4"],
    "auditory": ["auditory_cortex", "temporal_lobe"],
    "cognitive": ["prefrontal_cortex", "anterior_cingulate"],
    "emotion": ["amygdala", "insula", "anterior_cingulate"],
    "memory": ["hippocampus", "entorhinal_cortex"],
    "motor": ["motor_cortex", "basal_ganglia", "cerebellum"],
    "reward": ["nucleus_accumbens", "ventral_tegmental_area", "striatum"],
    "pain": ["anterior_cingulate", "insula", "thalamus"],
    "self": ["default_mode_network", "medial_prefrontal_cortex", "posterior_cingulate"],
    "dissociation": ["claustrum", "insula", "temporal_parietal_junction"],
    "physical": ["hypothalamus", "brainstem", "spinal_cord"],
}


def get_receptors_for_effect(effect_name: str) -> list:
    """Determine likely receptor involvement based on effect name."""
    effect_lower = effect_name.lower()
    receptors = set()
    
    for keyword, receptor_list in RECEPTOR_MAPPINGS.items():
        if keyword in effect_lower:
            receptors.update(receptor_list)
    
    # Default to 5-HT2A for psychedelic-type effects
    if not receptors and any(word in effect_lower for word in ["hallucin", "visual", "geometr", "perception"]):
        receptors.add("5-HT2A")
    
    return list(receptors) if receptors else ["5-HT2A"]  # Default


def get_regions_for_effect(effect_name: str, category: str) -> list:
    """Determine likely brain regions based on effect name and category."""
    effect_lower = effect_name.lower()
    category_lower = category.lower()
    regions = set()
    
    # Category-based mapping
    if "visual" in category_lower:
        regions.update(["visual_cortex", "V1", "V2"])
    elif "auditory" in category_lower:
        regions.update(["auditory_cortex", "temporal_lobe"])
    elif "cognitive" in category_lower:
        regions.update(["prefrontal_cortex", "anterior_cingulate"])
    elif "physical" in category_lower:
        regions.update(["hypothalamus", "brainstem"])
    elif "disconnective" in category_lower:
        regions.update(["claustrum", "insula"])
    
    # Effect-specific mapping
    for keyword, region_list in REGION_MAPPINGS.items():
        if keyword in effect_lower:
            regions.update(region_list)
    
    return list(regions) if regions else ["prefrontal_cortex"]


def scrape_effect_page(url: str) -> dict:
    """Scrape an individual effect page for details."""
    try:
        response = requests.get(url, timeout=10)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Get description from first paragraph
        content = soup.find('div', {'class': 'mw-parser-output'})
        if content:
            paragraphs = content.find_all('p', limit=3)
            description = ' '.join(p.get_text().strip() for p in paragraphs if p.get_text().strip())
            # Clean up
            description = re.sub(r'\[\d+\]', '', description)  # Remove citations
            description = description[:500] if len(description) > 500 else description
            return {"description": description}
    except Exception as e:
        print(f"  Error scraping {url}: {e}")
    
    return {"description": ""}


def scrape_effect_index():
    """Scrape the main Subjective Effect Index page."""
    print("Fetching Subjective Effect Index...")
    
    url = f"{BASE_URL}/wiki/Subjective_effect_index"
    response = requests.get(url, timeout=10)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    effects = {}
    effect_descriptions = {}
    current_category = "Uncategorized"
    current_subcategory = ""
    
    # Find all list items and headers
    content = soup.find('div', {'class': 'mw-parser-output'})
    
    if not content:
        print("Could not find content!")
        return {}, {}
    
    for element in content.find_all(['h3', 'h4', 'li', 'a']):
        # Track category from h3
        if element.name == 'h3':
            span = element.find('span', {'class': 'mw-headline'})
            if span:
                current_category = span.get_text().strip()
                current_subcategory = ""
                print(f"\n=== {current_category} ===")
        
        # Track subcategory from h4
        elif element.name == 'h4':
            span = element.find('span', {'class': 'mw-headline'})
            if span:
                current_subcategory = span.get_text().strip()
                print(f"  -- {current_subcategory} --")
        
        # Get effects from list items
        elif element.name == 'li':
            link = element.find('a')
            if link and link.get('href', '').startswith('/wiki/'):
                href = link.get('href')
                # Skip meta/category links
                if any(skip in href for skip in ['Category:', 'Special:', 'Project:', 'Talk:', 'File:']):
                    continue
                
                effect_name = link.get_text().strip()
                if not effect_name or len(effect_name) < 2:
                    continue
                
                # Create effect key (snake_case)
                effect_key = re.sub(r'[^a-z0-9]+', '_', effect_name.lower()).strip('_')
                
                if effect_key and effect_key not in effects:
                    full_category = f"{current_category} > {current_subcategory}" if current_subcategory else current_category
                    
                    effects[effect_key] = {
                        "name": effect_name,
                        "category": current_category,
                        "subcategory": current_subcategory,
                        "receptors": get_receptors_for_effect(effect_name),
                        "regions": get_regions_for_effect(effect_name, current_category),
                        "url": f"{BASE_URL}{href}"
                    }
                    
                    effect_descriptions[effect_key] = effect_name
                    print(f"    + {effect_name}")
    
    return effects, effect_descriptions


def fetch_effect_descriptions(effects: dict, limit: int = 50):
    """Fetch detailed descriptions for effects (rate-limited)."""
    print(f"\nFetching descriptions for up to {limit} effects...")
    
    count = 0
    for key, effect in effects.items():
        if count >= limit:
            break
        
        if 'url' in effect:
            print(f"  Fetching: {effect['name']}")
            details = scrape_effect_page(effect['url'])
            if details.get('description'):
                effect['description'] = details['description']
            count += 1
            time.sleep(0.5)  # Rate limit
    
    return effects


def generate_js_data(effects: dict, effect_descriptions: dict):
    """Generate JavaScript data file content."""
    
    # Group effects by category for the effect selector
    categories = {}
    for key, effect in effects.items():
        cat = effect.get('category', 'Other')
        if cat not in categories:
            categories[cat] = []
        categories[cat].append({
            "key": key,
            "name": effect['name'],
            "subcategory": effect.get('subcategory', '')
        })
    
    # Create effects object for receptor/region mapping
    effects_data = {}
    for key, effect in effects.items():
        effects_data[key] = {
            "name": effect['name'],
            "receptors": effect.get('receptors', []),
            "regions": effect.get('regions', []),
            "category": effect.get('category', ''),
            "subcategory": effect.get('subcategory', ''),
            "description": effect.get('description', effect['name'])
        }
    
    return {
        "effectCategories": categories,
        "effects": effects_data,
        "effectDescriptions": effect_descriptions
    }


def main():
    # Scrape the effect index
    effects, effect_descriptions = scrape_effect_index()
    
    print(f"\nFound {len(effects)} effects!")
    
    # Optionally fetch detailed descriptions (slow, rate-limited)
    # effects = fetch_effect_descriptions(effects, limit=30)
    
    # Save raw data
    with open('/data/users/zacharie/whc/consciousness_study/website/data/psychonautwiki_effects.json', 'w') as f:
        json.dump({
            "effects": effects,
            "descriptions": effect_descriptions
        }, f, indent=2)
    
    print(f"\nSaved to data/psychonautwiki_effects.json")
    
    # Generate JS-ready data
    js_data = generate_js_data(effects, effect_descriptions)
    
    with open('/data/users/zacharie/whc/consciousness_study/website/data/effects_full.json', 'w') as f:
        json.dump(js_data, f, indent=2)
    
    print(f"Saved JS-ready data to data/effects_full.json")
    
    # Print summary
    print("\n=== Summary by Category ===")
    for cat, items in js_data['effectCategories'].items():
        print(f"  {cat}: {len(items)} effects")


if __name__ == "__main__":
    main()
