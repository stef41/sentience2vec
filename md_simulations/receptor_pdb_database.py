#!/usr/bin/env python3
"""
Comprehensive Neuroreceptor PDB Database

Maps all 20 neuroreceptors from the consciousness study to their best available
crystal/cryo-EM structures from the Protein Data Bank (PDB).

Selection criteria:
- High resolution (<3.0 Å preferred)
- Human receptor preferred, then close homologs
- Active or inactive state as available
- With or without ligand bound
"""

# Complete neuroreceptor to PDB mapping
# Format: receptor_name -> {pdb_id, chain, resolution, state, species, notes}

RECEPTOR_PDB_DATABASE = {
    # ============== SEROTONIN RECEPTORS ==============
    '5-HT1A': {
        'pdb_id': '7E2Y',  # Human 5-HT1A with aripiprazole
        'chain': 'A',
        'resolution': 3.0,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P08908',
        'notes': '5-HT1A with aripiprazole, Gi coupled'
    },
    '5-HT2A': {
        'pdb_id': '6WHA',  # Human 5-HT2A with LSD (already have)
        'chain': 'A',
        'resolution': 3.4,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P28223',
        'notes': '5-HT2A with LSD, key psychedelic target'
    },
    '5-HT2B': {
        'pdb_id': '5TVN',  # Human 5-HT2B with ergotamine (already have)
        'chain': 'A',
        'resolution': 2.9,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P41595',
        'notes': '5-HT2B with ergotamine'
    },
    '5-HT2C': {
        'pdb_id': '6BQG',  # Human 5-HT2C with ergotamine
        'chain': 'A',
        'resolution': 3.0,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P28335',
        'notes': '5-HT2C with ergotamine'
    },
    
    # ============== DOPAMINE RECEPTORS ==============
    'D1': {
        'pdb_id': '7JVP',  # Human D1 with SKF-81297
        'chain': 'A',
        'resolution': 3.2,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P21728',
        'notes': 'D1 with agonist, Gs coupled'
    },
    'D2': {
        'pdb_id': '6CM4',  # Human D2 with risperidone (already have)
        'chain': 'A',
        'resolution': 2.87,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'P14416',
        'notes': 'D2 with risperidone, key antipsychotic target'
    },
    'D3': {
        'pdb_id': '3PBL',  # Human D3 with eticlopride
        'chain': 'A',
        'resolution': 2.89,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'P35462',
        'notes': 'D3 with eticlopride'
    },
    
    # ============== OPIOID RECEPTORS ==============
    'MOR': {
        'pdb_id': '5C1M',  # Mouse MOR with BU72 (already have)
        'chain': 'A',
        'resolution': 2.1,
        'state': 'active',
        'species': 'mouse',
        'uniprot': 'P42866',
        'notes': 'MOR active with BU72 agonist, Gi coupled'
    },
    'KOR': {
        'pdb_id': '4DJH',  # Human KOR with JDTic (already have)
        'chain': 'A',
        'resolution': 2.9,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'P41145',
        'notes': 'KOR with JDTic antagonist, salvinorin target'
    },
    'DOR': {
        'pdb_id': '4N6H',  # Human DOR with naltrindole
        'chain': 'A',
        'resolution': 1.8,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'P41143',
        'notes': 'DOR with naltrindole, high resolution'
    },
    
    # ============== CANNABINOID RECEPTORS ==============
    'CB1': {
        'pdb_id': '5TGZ',  # Human CB1 with AM6538 (already have)
        'chain': 'A',
        'resolution': 2.8,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'P21554',
        'notes': 'CB1 inactive, THC target'
    },
    
    # ============== GLUTAMATE RECEPTORS ==============
    'NMDA_GluN2A': {
        'pdb_id': '7EU7',  # GluN1/GluN2A NMDA receptor
        'chain': 'A',
        'resolution': 3.2,
        'state': 'active',
        'species': 'human/rat',
        'uniprot': 'Q12879',
        'notes': 'NMDA receptor GluN2A subunit, ketamine target'
    },
    'NMDA_GluN2B': {
        'pdb_id': '7SAD',  # GluN1/GluN2B with ifenprodil
        'chain': 'B',
        'resolution': 2.8,
        'state': 'inhibited',
        'species': 'rat',
        'uniprot': 'Q00960',
        'notes': 'NMDA GluN2B subunit with allosteric modulator'
    },
    'AMPA_GluA1': {
        'pdb_id': '5WEO',  # GluA1 AMPA receptor
        'chain': 'A',
        'resolution': 3.4,
        'state': 'desensitized',
        'species': 'rat',
        'uniprot': 'P19490',
        'notes': 'AMPA GluA1 subunit'
    },
    
    # ============== GABA RECEPTORS ==============
    'GABA_A_alpha1': {
        'pdb_id': '6HUP',  # Human GABA-A alpha1/beta2/gamma2
        'chain': 'A',
        'resolution': 3.1,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P14867',
        'notes': 'GABA-A alpha1 subunit, benzodiazepine site'
    },
    'GABA_A_alpha2': {
        'pdb_id': '6X3T',  # GABA-A with diazepam
        'chain': 'A',
        'resolution': 3.2,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P47869',
        'notes': 'GABA-A alpha2 subunit'
    },
    
    # ============== MUSCARINIC RECEPTORS ==============
    'M1_muscarinic': {
        'pdb_id': '6OIJ',  # Human M1 with iperoxo
        'chain': 'A',
        'resolution': 3.0,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P11229',
        'notes': 'M1 muscarinic, deliriant target'
    },
    
    # ============== ADRENERGIC RECEPTORS ==============
    'alpha2A': {
        'pdb_id': '6KUX',  # Human alpha2A with dexmedetomidine
        'chain': 'A',
        'resolution': 2.9,
        'state': 'active',
        'species': 'human',
        'uniprot': 'P08913',
        'notes': 'Alpha2A adrenergic with dexmedetomidine'
    },
    
    # ============== NICOTINIC RECEPTORS ==============
    'nAChR_alpha4beta2': {
        'pdb_id': '6CNJ',  # Human alpha4beta2 nicotinic
        'chain': 'A',
        'resolution': 3.9,
        'state': 'desensitized',
        'species': 'human',
        'uniprot': 'P43681',
        'notes': 'Nicotinic alpha4beta2, nicotine target'
    },
    
    # ============== SIGMA RECEPTORS ==============
    'Sigma1': {
        'pdb_id': '5HK1',  # Human sigma-1 with PD144418
        'chain': 'A',
        'resolution': 2.5,
        'state': 'inactive',
        'species': 'human',
        'uniprot': 'Q99720',
        'notes': 'Sigma-1 receptor, DMT and other psychedelic target'
    },
}

# Alternative structures for some receptors (backup options)
ALTERNATIVE_STRUCTURES = {
    '5-HT1A': ['7E2X', '7E32'],  # Other 5-HT1A structures
    '5-HT2A': ['6WHA', '6WH4'],  # 6WH4 is with 25-CN-NBOH
    'D2': ['7JVR', '6VMS'],  # Other D2 structures
    'MOR': ['6DDF', '6DDE'],  # Active structures with different ligands
    'KOR': ['6B73', '6VI4'],  # Active KOR structures
    'CB1': ['5XRA', '6N4B'],  # CB1 with different ligands
}


def get_all_receptor_ids():
    """Return list of all PDB IDs needed."""
    return [v['pdb_id'] for v in RECEPTOR_PDB_DATABASE.values()]


def get_receptor_info(receptor_name):
    """Get PDB info for a specific receptor."""
    return RECEPTOR_PDB_DATABASE.get(receptor_name)


def print_database_summary():
    """Print a summary of the receptor database."""
    print("=" * 80)
    print("NEURORECEPTOR PDB DATABASE")
    print("=" * 80)
    
    # Group by receptor family
    families = {
        'Serotonin': ['5-HT1A', '5-HT2A', '5-HT2B', '5-HT2C'],
        'Dopamine': ['D1', 'D2', 'D3'],
        'Opioid': ['MOR', 'KOR', 'DOR'],
        'Cannabinoid': ['CB1'],
        'Glutamate': ['NMDA_GluN2A', 'NMDA_GluN2B', 'AMPA_GluA1'],
        'GABA': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'Muscarinic': ['M1_muscarinic'],
        'Adrenergic': ['alpha2A'],
        'Nicotinic': ['nAChR_alpha4beta2'],
        'Sigma': ['Sigma1'],
    }
    
    for family, receptors in families.items():
        print(f"\n{family} Receptors:")
        print("-" * 40)
        for rec in receptors:
            info = RECEPTOR_PDB_DATABASE.get(rec, {})
            print(f"  {rec:20s} -> {info.get('pdb_id', 'N/A'):6s} "
                  f"({info.get('resolution', '?'):.1f}Å, {info.get('species', '?')})")
    
    print(f"\n{'=' * 80}")
    print(f"Total receptors: {len(RECEPTOR_PDB_DATABASE)}")
    print(f"Unique PDB structures: {len(set(get_all_receptor_ids()))}")


if __name__ == '__main__':
    print_database_summary()
