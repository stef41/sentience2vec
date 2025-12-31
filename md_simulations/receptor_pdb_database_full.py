#!/usr/bin/env python3
"""
COMPREHENSIVE Neuroreceptor PDB Database - ALL receptor subtypes

Maps ALL known neuroreceptor subtypes relevant to psychoactive compounds
to their best available crystal/cryo-EM structures from PDB.
"""

# Complete neuroreceptor database - 77 receptor subtypes
RECEPTOR_PDB_DATABASE = {
    # ==================== SEROTONIN (14 subtypes) ====================
    '5-HT1A': {'pdb_id': '7E2Y', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    '5-HT1B': {'pdb_id': '4IAR', 'chain': 'A', 'resolution': 2.7, 'species': 'human'},
    '5-HT1D': {'pdb_id': '7E32', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
    '5-HT1E': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None, 'notes': 'No structure available'},
    '5-HT1F': {'pdb_id': '7EXD', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    '5-HT2A': {'pdb_id': '6WHA', 'chain': 'A', 'resolution': 3.4, 'species': 'human'},
    '5-HT2B': {'pdb_id': '5TVN', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    '5-HT2C': {'pdb_id': '6BQG', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    '5-HT3A': {'pdb_id': '6BE1', 'chain': 'A', 'resolution': 3.5, 'species': 'mouse'},  # Ion channel
    '5-HT4': {'pdb_id': '7XT8', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    '5-HT5A': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None, 'notes': 'No structure'},
    '5-HT6': {'pdb_id': '7XTC', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    '5-HT7': {'pdb_id': '7XTB', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    
    # ==================== DOPAMINE (5 subtypes) ====================
    'D1': {'pdb_id': '7JVP', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'D2': {'pdb_id': '6CM4', 'chain': 'A', 'resolution': 2.87, 'species': 'human'},
    'D3': {'pdb_id': '3PBL', 'chain': 'A', 'resolution': 2.89, 'species': 'human'},
    'D4': {'pdb_id': '5WIU', 'chain': 'A', 'resolution': 2.14, 'species': 'human'},
    'D5': {'pdb_id': '7LJC', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    
    # ==================== OPIOID (4 subtypes) ====================
    'MOR': {'pdb_id': '5C1M', 'chain': 'A', 'resolution': 2.1, 'species': 'mouse'},
    'KOR': {'pdb_id': '4DJH', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    'DOR': {'pdb_id': '4N6H', 'chain': 'A', 'resolution': 1.8, 'species': 'human'},
    'NOP': {'pdb_id': '5DHG', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},  # Nociceptin
    
    # ==================== ADRENERGIC (9 subtypes) ====================
    'alpha1A': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None},
    'alpha1B': {'pdb_id': '7B6W', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
    'alpha1D': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None},
    'alpha2A': {'pdb_id': '6KUX', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    'alpha2B': {'pdb_id': '6K41', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    'alpha2C': {'pdb_id': '6KUW', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    'beta1': {'pdb_id': '7BTS', 'chain': 'A', 'resolution': 2.7, 'species': 'human'},
    'beta2': {'pdb_id': '3SN6', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'beta3': {'pdb_id': '7DHI', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    
    # ==================== MUSCARINIC (5 subtypes) ====================
    'M1': {'pdb_id': '6OIJ', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    'M2': {'pdb_id': '5ZKC', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    'M3': {'pdb_id': '4U15', 'chain': 'A', 'resolution': 2.8, 'species': 'rat'},
    'M4': {'pdb_id': '5DSG', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    'M5': {'pdb_id': '6OL9', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    
    # ==================== HISTAMINE (4 subtypes) ====================
    'H1': {'pdb_id': '3RZE', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
    'H2': {'pdb_id': '7UL3', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    'H3': {'pdb_id': '7F61', 'chain': 'A', 'resolution': 3.5, 'species': 'human'},
    'H4': {'pdb_id': '8JE2', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    
    # ==================== CANNABINOID (2 subtypes) ====================
    'CB1': {'pdb_id': '5TGZ', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    'CB2': {'pdb_id': '5ZTY', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    
    # ==================== GLUTAMATE - Ionotropic (6) ====================
    'NMDA_GluN1': {'pdb_id': '7EU7', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'NMDA_GluN2A': {'pdb_id': '7EU7', 'chain': 'B', 'resolution': 3.2, 'species': 'rat'},
    'NMDA_GluN2B': {'pdb_id': '7SAD', 'chain': 'B', 'resolution': 2.8, 'species': 'rat'},
    'AMPA_GluA1': {'pdb_id': '5WEO', 'chain': 'A', 'resolution': 3.4, 'species': 'rat'},
    'AMPA_GluA2': {'pdb_id': '5WEK', 'chain': 'A', 'resolution': 3.5, 'species': 'rat'},
    'Kainate_GluK2': {'pdb_id': '5KUF', 'chain': 'A', 'resolution': 3.8, 'species': 'rat'},
    
    # ==================== GLUTAMATE - Metabotropic (8) ====================
    'mGluR1': {'pdb_id': '4OR2', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    'mGluR2': {'pdb_id': '7MTS', 'chain': 'A', 'resolution': 2.5, 'species': 'human'},
    'mGluR3': {'pdb_id': '7WUE', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'mGluR4': {'pdb_id': '7E9G', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    'mGluR5': {'pdb_id': '6N52', 'chain': 'A', 'resolution': 2.6, 'species': 'human'},
    'mGluR6': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None},
    'mGluR7': {'pdb_id': '7EPA', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    'mGluR8': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None},
    
    # ==================== GABA (5 subtypes) ====================
    'GABA_A_alpha1': {'pdb_id': '6HUP', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
    'GABA_A_alpha2': {'pdb_id': '6X3T', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'GABA_A_alpha3': {'pdb_id': '7QNE', 'chain': 'A', 'resolution': 2.9, 'species': 'human'},
    'GABA_A_alpha5': {'pdb_id': '7QM5', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    'GABA_B': {'pdb_id': '7C7Q', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    
    # ==================== NICOTINIC (4 subtypes) ====================
    'nAChR_alpha4beta2': {'pdb_id': '6CNJ', 'chain': 'A', 'resolution': 3.9, 'species': 'human'},
    'nAChR_alpha7': {'pdb_id': '7KOO', 'chain': 'A', 'resolution': 2.7, 'species': 'human'},
    'nAChR_alpha3beta4': {'pdb_id': '6PV7', 'chain': 'A', 'resolution': 3.3, 'species': 'human'},
    'nAChR_muscle': {'pdb_id': '7SMM', 'chain': 'A', 'resolution': 2.5, 'species': 'human'},
    
    # ==================== SIGMA (2 subtypes) ====================
    'Sigma1': {'pdb_id': '5HK1', 'chain': 'A', 'resolution': 2.5, 'species': 'human'},
    'Sigma2': {'pdb_id': '7M95', 'chain': 'A', 'resolution': 2.2, 'species': 'human'},
    
    # ==================== ADENOSINE (4 subtypes) ====================
    'A1': {'pdb_id': '5UEN', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    'A2A': {'pdb_id': '5IU4', 'chain': 'A', 'resolution': 1.7, 'species': 'human'},
    'A2B': {'pdb_id': '6LPJ', 'chain': 'A', 'resolution': 3.4, 'species': 'human'},
    'A3': {'pdb_id': '8HN1', 'chain': 'A', 'resolution': 2.2, 'species': 'human'},
    
    # ==================== MELATONIN (2 subtypes) ====================
    'MT1': {'pdb_id': '6ME2', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    'MT2': {'pdb_id': '6ME4', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    
    # ==================== OREXIN (2 subtypes) ====================
    'OX1': {'pdb_id': '6TOD', 'chain': 'A', 'resolution': 2.8, 'species': 'human'},
    'OX2': {'pdb_id': '5WQC', 'chain': 'A', 'resolution': 2.5, 'species': 'human'},
    
    # ==================== TRACE AMINE ====================
    'TAAR1': {'pdb_id': '8JLN', 'chain': 'A', 'resolution': 3.2, 'species': 'human'},
    
    # ==================== TRP CHANNELS ====================
    'TRPV1': {'pdb_id': '7RQW', 'chain': 'A', 'resolution': 2.7, 'species': 'human'},
    'TRPM8': {'pdb_id': '6BPQ', 'chain': 'A', 'resolution': 4.1, 'species': 'human'},
    
    # ==================== GLYCINE ====================
    'GlyR_alpha1': {'pdb_id': '6PM6', 'chain': 'A', 'resolution': 3.0, 'species': 'human'},
    
    # ==================== IMIDAZOLINE ====================
    'I1': {'pdb_id': None, 'chain': None, 'resolution': None, 'species': None, 'notes': 'Binding site on alpha2A'},
    
    # ==================== TRANSPORTER TARGETS ====================
    'SERT': {'pdb_id': '6VRH', 'chain': 'A', 'resolution': 3.3, 'species': 'human'},
    'DAT': {'pdb_id': '4XP4', 'chain': 'A', 'resolution': 2.95, 'species': 'drosophila'},
    'NET': {'pdb_id': '7LI8', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
    'VMAT2': {'pdb_id': '7YOC', 'chain': 'A', 'resolution': 3.1, 'species': 'human'},
}

def count_available():
    available = sum(1 for v in RECEPTOR_PDB_DATABASE.values() if v.get('pdb_id'))
    return available, len(RECEPTOR_PDB_DATABASE)

if __name__ == '__main__':
    avail, total = count_available()
    print(f"Receptors with PDB structures: {avail}/{total}")
    
    # By family
    families = {}
    for name in RECEPTOR_PDB_DATABASE:
        if name.startswith('5-HT'): fam = 'Serotonin'
        elif name.startswith('D') and name[1:].isdigit(): fam = 'Dopamine'
        elif name in ['MOR', 'KOR', 'DOR', 'NOP']: fam = 'Opioid'
        elif name.startswith('alpha') or name.startswith('beta'): fam = 'Adrenergic'
        elif name.startswith('M') and len(name) == 2: fam = 'Muscarinic'
        elif name.startswith('H') and len(name) == 2: fam = 'Histamine'
        elif name.startswith('CB'): fam = 'Cannabinoid'
        elif 'Glu' in name or 'mGlu' in name or 'NMDA' in name or 'AMPA' in name or 'Kainate' in name: fam = 'Glutamate'
        elif 'GABA' in name: fam = 'GABA'
        elif name.startswith('nAChR'): fam = 'Nicotinic'
        elif name.startswith('Sigma'): fam = 'Sigma'
        elif name.startswith('A') and len(name) <= 3: fam = 'Adenosine'
        elif name.startswith('MT'): fam = 'Melatonin'
        elif name.startswith('OX'): fam = 'Orexin'
        elif name.startswith('TRPV') or name.startswith('TRPM'): fam = 'TRP'
        elif name in ['SERT', 'DAT', 'NET', 'VMAT2']: fam = 'Transporter'
        else: fam = 'Other'
        families[fam] = families.get(fam, 0) + 1
    
    print("\nBy family:")
    for fam, count in sorted(families.items(), key=lambda x: -x[1]):
        print(f"  {fam}: {count}")
