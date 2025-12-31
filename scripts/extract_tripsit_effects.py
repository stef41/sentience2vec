#!/usr/bin/env python3
"""
Extract drug-effect relationships from TripSit data.
TripSit provides effect descriptions for many drugs.
"""

import json
import re
from pathlib import Path
from collections import defaultdict

DATA_PATH = Path('/data/users/zacharie/whc/consciousness_study/website/data')

# Load TripSit data
with open(DATA_PATH / 'tripsit_data.json') as f:
    tripsit_data = json.load(f)

print(f"Loading TripSit data: {len(tripsit_data)} drugs")

# Extract effects from each drug
drug_effects = {}
effect_drugs = defaultdict(list)

# Common effect names to normalize
EFFECT_NORMALIZATIONS = {
    'open eye visuals': 'open_eye_visuals',
    'closed eye visuals': 'closed_eye_visuals',
    'open/closed eye visuals': 'visual_hallucinations',
    'closed/open eye visuals': 'visual_hallucinations',
    'cev': 'closed_eye_visuals',
    'oev': 'open_eye_visuals',
    'mood lift': 'mood_enhancement',
    'mood enhancement': 'mood_enhancement',
    'increased energy': 'stimulation',
    'energy': 'stimulation',
    'alertness': 'wakefulness',
    'increased alertness': 'wakefulness',
    'decreased need for sleep': 'wakefulness',
    'increased sociability': 'sociability_enhancement',
    'sociability': 'sociability_enhancement',
    'enhanced tactile sensation': 'tactile_enhancement',
    'tactile enhancement': 'tactile_enhancement',
    'brightened color': 'color_enhancement',
    'brightened colours': 'color_enhancement',
    'brightened colour': 'color_enhancement',
    'color enhancement': 'color_enhancement',
    'colour enhancement': 'color_enhancement',
    'visual enhancement': 'visual_enhancement',
    'mental stimulation': 'thought_acceleration',
    'physical stimulation': 'stimulation',
    'cognitive enhancement': 'analysis_enhancement',
    'increased focus': 'focus_enhancement',
    'focus': 'focus_enhancement',
    'increased confidence': 'disinhibition',
    'confidence': 'disinhibition',
    'relaxation': 'sedation',
    'sedative': 'sedation',
    'muscle relaxant': 'muscle_relaxation',
    'hypnotic': 'sedation',
    'anxiolytic': 'anxiety_suppression',
    'amnesia': 'amnesia',
    'nausea': 'nausea',
    'sweating': 'increased_perspiration',
    'increased sweating': 'increased_perspiration',
    'jaw clenching': 'teeth_grinding',
    'teeth grinding': 'teeth_grinding',
    'bruxism': 'teeth_grinding',
    'appetite suppression': 'appetite_suppression',
    'appetite loss': 'appetite_suppression',
    'analgesia': 'pain_relief',
    'pain relief': 'pain_relief',
    'itchiness': 'itchiness',
    'itching': 'itchiness',
    'respiratory depression': 'respiratory_depression',
    'constipation': 'constipation',
    'dry mouth': 'dry_mouth',
    'pupil dilation': 'pupil_dilation',
    'dilated pupils': 'pupil_dilation',
    'pupil constriction': 'pupil_constriction',
    'constricted pupils': 'pupil_constriction',
    'miosis': 'pupil_constriction',
    'mydriasis': 'pupil_dilation',
    'increased heart rate': 'increased_heart_rate',
    'tachycardia': 'increased_heart_rate',
    'increased blood pressure': 'increased_blood_pressure',
    'hypertension': 'increased_blood_pressure',
    'vasoconstriction': 'vasoconstriction',
    'time dilation': 'time_distortion',
    'time distortion': 'time_distortion',
    'altered time perception': 'time_distortion',
    'ego death': 'ego_death',
    'ego dissolution': 'ego_death',
    'ego loss': 'ego_death',
    'depersonalization': 'depersonalization',
    'derealization': 'derealization',
    'dissociation': 'cognitive_disconnection',
    'auditory hallucinations': 'auditory_hallucination',
    'auditory enhancement': 'auditory_enhancement',
    'music enhancement': 'increased_music_appreciation',
    'enhanced music appreciation': 'increased_music_appreciation',
    'giggling': 'laughter_fits',
    'laughter': 'laughter_fits',
    'paranoia': 'paranoia',
    'anxiety': 'anxiety',
    'confusion': 'confusion',
    'memory impairment': 'memory_suppression',
    'short term memory loss': 'memory_suppression',
    'thought loops': 'thought_loop',
    'thought loop': 'thought_loop',
    'introspection': 'increased_introspection',
    'increased introspection': 'increased_introspection',
    'insight': 'increased_introspection',
    'creativity': 'creativity_enhancement',
    'enhanced creativity': 'creativity_enhancement',
    'increased libido': 'increased_libido',
    'sexual arousal': 'increased_libido',
    'tactile sensations': 'tactile_enhancement',
    'body high': 'physical_euphoria',
    'body load': 'bodily_pressures',
    'nystagmus': 'double_vision',
}

def normalize_effect(effect_text):
    """Normalize effect text to standard ID."""
    effect_lower = effect_text.lower().strip()
    
    # Check direct mapping
    if effect_lower in EFFECT_NORMALIZATIONS:
        return EFFECT_NORMALIZATIONS[effect_lower]
    
    # Convert to snake_case
    normalized = re.sub(r'[^a-z0-9\s]', '', effect_lower)
    normalized = re.sub(r'\s+', '_', normalized)
    return normalized

def extract_effects_from_text(effects_text):
    """Extract individual effects from a comma/semicolon separated list."""
    if not effects_text:
        return []
    
    # Split by comma, semicolon, or period
    parts = re.split(r'[,;.]', str(effects_text))
    
    effects = []
    for part in parts:
        part = part.strip()
        if not part or len(part) < 3:
            continue
        
        # Skip dosage info, duration, notes
        skip_patterns = ['mg', 'ug', 'hour', 'minute', 'duration', 'onset', 'dose', 'note:', 
                        'warning', 'caution', 'avoid', 'contraindicated', 'similar to', 
                        'comparable to', 'described as', '...']
        if any(s in part.lower() for s in skip_patterns):
            continue
        
        # Normalize and add
        normalized = normalize_effect(part)
        if normalized and len(normalized) > 2:
            effects.append(normalized)
    
    return effects

# Process each drug
for drug_id, drug in tripsit_data.items():
    drug_name = drug['name']
    
    # Get effects from properties
    effects_text = drug.get('effects', '') or drug.get('properties', {}).get('effects', '')
    
    if effects_text:
        effects = extract_effects_from_text(effects_text)
        
        if effects:
            drug_effects[drug_name] = {
                'effects': effects,
                'raw_effects': effects_text,
                'source': 'TripSit'
            }
            
            # Build reverse mapping
            for effect in effects:
                if drug_name not in effect_drugs[effect]:
                    effect_drugs[effect].append(drug_name)

print(f"\nExtracted effects from {len(drug_effects)} drugs")
print(f"Found {len(effect_drugs)} unique effects")

# Save TripSit effect data
with open(DATA_PATH / 'tripsit_effects.json', 'w') as f:
    json.dump({
        'drug_effects': drug_effects,
        'effect_drugs': dict(effect_drugs)
    }, f, indent=2)

print(f"Saved to tripsit_effects.json")

# Show top effects
print("\n=== Top TripSit Effects (by drug count) ===")
sorted_effects = sorted(effect_drugs.items(), key=lambda x: len(x[1]), reverse=True)
for effect, drugs in sorted_effects[:20]:
    print(f"  {effect}: {len(drugs)} drugs")

# Show sample drug effects
print("\n=== Sample Drug Effects ===")
for drug_name in list(drug_effects.keys())[:5]:
    data = drug_effects[drug_name]
    print(f"\n{drug_name}:")
    print(f"  Effects: {', '.join(data['effects'][:8])}")
