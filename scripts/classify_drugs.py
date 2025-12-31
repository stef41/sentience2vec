#!/usr/bin/env python3
"""
Add proper drug classifications based on known pharmacology.
"""

import json
import re

# Known drug classifications
DRUG_CLASSES = {
    # Classic psychedelics (5-HT2A agonists)
    'psychedelic': [
        'LSD', 'Psilocybin', 'Psilocin', 'DMT', 'Mescaline', '5-MeO-DMT',
        '2C-B', '2C-E', '2C-I', '2C-C', '2C-D', '2C-T-2', '2C-T-7',
        'DOB', 'DOC', 'DOI', 'DOM',
        '25I-NBOMe', '25C-NBOMe', '25B-NBOMe',
        'AL-LAD', 'ETH-LAD', '1P-LSD', 'ALD-52', 'LSA',
        '4-AcO-DMT', '4-HO-MET', '4-HO-MiPT', '4-AcO-MET',
        'DPT', 'DiPT', '5-MeO-DiPT', '5-MeO-MiPT',
        'TMA', 'TMA-2', 'Escaline', 'Allylescaline',
        '2C-B-FLY', 'Bromo-DragonFLY',
        'Ibogaine', 'Salvinorin A', 'Ayahuasca', 'Bufotenin',
    ],
    
    # Dissociatives (NMDA antagonists, sigma agonists)
    'dissociative': [
        'Ketamine', 'PCP', 'DXM', 'Dextromethorphan', 'MXE', 'Methoxetamine',
        '3-MeO-PCP', '3-HO-PCP', '3-MeO-PCE', 'O-PCE',
        'DCK', 'Deschloroketamine', '2-FDCK', '2F-DCK',
        'Nitrous Oxide', 'N2O', 'Memantine', 'Dizocilpine', 'MK-801',
        'Ephenidine', 'Diphenidine', '3-HO-PCE',
    ],
    
    # Empathogens/Entactogens (SERT releasers)
    'empathogen': [
        'MDMA', 'MDA', 'MDEA', 'MBDB', '6-APB', '5-APB', '5-MAPB', '6-MAPB',
        'Methylone', 'Butylone', 'Ethylone', 'Mephedrone',
        '4-FA', '4-FMA', '5-EAPB', '5-APDB', '6-APDB',
        'BK-MDMA', 'Eutylone',
    ],
    
    # Stimulants (DAT/NET releasers/inhibitors)
    'stimulant': [
        'Amphetamine', 'Methamphetamine', 'Cocaine', 'Methylphenidate',
        'Lisdexamfetamine', 'Dextroamphetamine', 'Adderall',
        'Modafinil', 'Armodafinil', 'Caffeine', 'Nicotine',
        'MDPV', 'A-PVP', 'A-PHP', 'Pyrovalerone',
        '2-FMA', '3-FMA', '3-FA', '4-FMA',
        'Ephedrine', 'Pseudoephedrine', 'Phenylephrine',
        'Cathinone', 'Methcathinone', 'Flephedrone',
        'NEP', 'Hexen', 'NEH', 'Pentedrone', 'Pentylone',
        'Propylhexedrine', 'Benzedrex',
    ],
    
    # Opioids (MOR/KOR/DOR agonists)
    'opioid': [
        'Morphine', 'Heroin', 'Codeine', 'Oxycodone', 'Hydrocodone',
        'Fentanyl', 'Carfentanil', 'Sufentanil', 'Alfentanil',
        'Tramadol', 'Tapentadol', 'Methadone', 'Buprenorphine',
        'Kratom', 'Mitragynine', '7-Hydroxymitragynine',
        'O-DSMT', 'O-Desmethyltramadol',
        'Hydromorphone', 'Oxymorphone', 'Meperidine',
        'U-47700', 'Isotonitazene', 'Etonitazene',
    ],
    
    # Cannabinoids (CB1/CB2 agonists)
    'cannabinoid': [
        'THC', 'Delta-9-THC', 'Delta-8-THC', 'THCA', 'THCV',
        'CBD', 'CBN', 'CBG', 'CBC',
        'JWH-018', 'JWH-073', 'JWH-250', 'AM-2201',
        'Synthetic Cannabinoid', 'Spice', 'K2',
        'WIN 55,212-2', 'CP 55,940', 'HU-210',
        'AB-FUBINACA', 'ADB-BUTINACA', 'MDMB-4en-PINACA',
        'Cannabis', 'Marijuana', 'Hashish',
    ],
    
    # Depressants (GABA modulators)
    'depressant': [
        'Alcohol', 'Ethanol',
        'Diazepam', 'Alprazolam', 'Clonazepam', 'Lorazepam',
        'Temazepam', 'Triazolam', 'Midazolam', 'Flurazepam',
        'Etizolam', 'Clonazolam', 'Flualprazolam', 'Flunitrazolam',
        'GHB', 'GBL', '1,4-Butanediol', 'Phenibut',
        'Barbiturate', 'Phenobarbital', 'Secobarbital',
        'Zolpidem', 'Zopiclone', 'Eszopiclone', 'Zaleplon',
        'Carisoprodol', 'Meprobamate', 'Methaqualone',
        'Baclofen', 'Pregabalin', 'Gabapentin',
    ],
    
    # Deliriants (mAChR antagonists)
    'deliriant': [
        'Diphenhydramine', 'DPH', 'Benadryl',
        'Datura', 'Scopolamine', 'Atropine', 'Hyoscyamine',
        'Dimenhydrinate', 'Dramamine', 'Cyclizine',
        'Trihexyphenidyl', 'Benztropine', 'Biperiden',
        'Amanita muscaria', 'Muscimol', 'Ibotenic acid',
    ],
    
    # Nootropics
    'nootropic': [
        'Piracetam', 'Aniracetam', 'Oxiracetam', 'Pramiracetam',
        'Phenylpiracetam', 'Noopept', 'Sunifiram', 'Unifiram',
        'Adrafinil', 'Fasoracetam', 'Coluracetam',
        'Semax', 'Selank', 'NSI-189',
        'Alpha-GPC', 'CDP-Choline', 'DMAE',
        'L-Theanine', 'Bacopa', 'Lion\'s Mane',
    ],
}

# Create reverse lookup
DRUG_TO_CLASS = {}
for drug_class, drugs in DRUG_CLASSES.items():
    for drug in drugs:
        DRUG_TO_CLASS[drug.lower()] = drug_class
        # Also add without spaces/hyphens
        DRUG_TO_CLASS[drug.lower().replace(' ', '').replace('-', '')] = drug_class

def classify_drug(drug_name):
    """Classify a drug by name."""
    name_lower = drug_name.lower()
    
    # Direct match
    if name_lower in DRUG_TO_CLASS:
        return DRUG_TO_CLASS[name_lower]
    
    # Try without hyphens/spaces
    simplified = name_lower.replace(' ', '').replace('-', '').replace('_', '')
    if simplified in DRUG_TO_CLASS:
        return DRUG_TO_CLASS[simplified]
    
    # Pattern matching
    if any(x in name_lower for x in ['lsd', 'lyserg']):
        return 'psychedelic'
    if any(x in name_lower for x in ['2c-', '25i', '25c', '25b', 'nbome']):
        return 'psychedelic'
    if any(x in name_lower for x in ['dmt', 'tryptamine', 'psilocybin', 'psilocin']):
        return 'psychedelic'
    if any(x in name_lower for x in ['mescaline', 'phenethylamine']):
        return 'psychedelic'
    if any(x in name_lower for x in ['ketamine', 'pcp', 'pce', 'dxm', 'mxe']):
        return 'dissociative'
    if any(x in name_lower for x in ['mdma', 'mda', 'apb', 'mapb', 'entactogen']):
        return 'empathogen'
    if any(x in name_lower for x in ['amphetamine', 'methamphet', 'cocaine', 'methylphen']):
        return 'stimulant'
    if any(x in name_lower for x in ['cathinone', 'pvp', 'php', 'fma']):
        return 'stimulant'
    if any(x in name_lower for x in ['morphine', 'codeine', 'fentanyl', 'oxycodone', 'hydrocodone']):
        return 'opioid'
    if any(x in name_lower for x in ['tramadol', 'methadone', 'buprenorphine', 'kratom']):
        return 'opioid'
    if any(x in name_lower for x in ['thc', 'cannabinoid', 'jwh', 'cannabis']):
        return 'cannabinoid'
    if any(x in name_lower for x in ['diazepam', 'alprazolam', 'clonazepam', 'lorazepam', 'azolam', 'azepam']):
        return 'depressant'
    if any(x in name_lower for x in ['ghb', 'phenibut', 'barbital', 'barbiturate']):
        return 'depressant'
    if any(x in name_lower for x in ['diphenhydramine', 'dph', 'scopolamine', 'atropine', 'datura']):
        return 'deliriant'
    if any(x in name_lower for x in ['racetam', 'noopept', 'modafinil']):
        return 'nootropic'
    
    return 'other'


# Read data.js
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'r') as f:
    content = f.read()

# Find drugs section and update classes
# Pattern: "DrugName": { ... "class": "other" ... }
def update_drug_class(match):
    drug_block = match.group(0)
    # Extract drug name from the match
    name_match = re.search(r'"([^"]+)":\s*\{', drug_block)
    if name_match:
        drug_name = name_match.group(1)
        new_class = classify_drug(drug_name)
        if new_class != 'other':
            # Replace the class
            drug_block = re.sub(r'"class":\s*"[^"]*"', f'"class": "{new_class}"', drug_block)
            print(f"  {drug_name} -> {new_class}")
    return drug_block

# This is complex because the structure is nested. Let's do a simpler approach
# Find all drug entries and update them

# First, let's count and update
updated_count = 0
lines = content.split('\n')
new_lines = []

in_drugs_section = False
current_drug = None

for i, line in enumerate(lines):
    # Check if we're entering drugs section
    if 'drugs:' in line and '{' in line:
        in_drugs_section = True
    
    # Track current drug name
    if in_drugs_section:
        drug_match = re.match(r'\s*"([^"]+)":\s*\{', line)
        if drug_match and '"class"' not in line:
            current_drug = drug_match.group(1)
        
        # Update class if found
        if '"class":' in line and current_drug:
            new_class = classify_drug(current_drug)
            if new_class != 'other' and '"other"' in line:
                line = line.replace('"other"', f'"{new_class}"')
                updated_count += 1
                print(f"  {current_drug} -> {new_class}")
    
    new_lines.append(line)

# Write back
with open('/data/users/zacharie/whc/consciousness_study/website/data.js', 'w') as f:
    f.write('\n'.join(new_lines))

print(f"\nUpdated {updated_count} drug classifications")
