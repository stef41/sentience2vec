#!/usr/bin/env python3
"""
Comprehensive Psychoactive Ligand SMILES Database

Contains SMILES strings for all major psychoactive substances relevant to
consciousness research, organized by pharmacological class.
"""

# Comprehensive ligand database with SMILES
# Format: compound_name -> {smiles, class, targets, molecular_weight, notes}

LIGAND_DATABASE = {
    # ========================================================================
    # CLASSICAL PSYCHEDELICS - Tryptamines
    # ========================================================================
    'LSD': {
        'smiles': 'CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(C2=C1)c34)C',
        'class': 'lysergamide',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C', '5-HT1A', 'D2'],
        'mw': 323.43,
        'notes': 'Prototypical psychedelic'
    },
    'psilocin': {
        'smiles': 'CN(C)CCc1c[nH]c2cccc(O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C', '5-HT1A'],
        'mw': 204.27,
        'notes': 'Active metabolite of psilocybin'
    },
    'psilocybin': {
        'smiles': 'CN(C)CCc1c[nH]c2cccc(OP(=O)(O)O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C'],
        'mw': 284.25,
        'notes': 'Prodrug of psilocin'
    },
    'DMT': {
        'smiles': 'CN(C)CCc1c[nH]c2ccccc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C', 'Sigma1'],
        'mw': 188.27,
        'notes': 'N,N-Dimethyltryptamine, endogenous'
    },
    '5-MeO-DMT': {
        'smiles': 'CN(C)CCc1c[nH]c2ccc(OC)cc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT1A', 'Sigma1'],
        'mw': 218.30,
        'notes': '5-Methoxy-DMT, potent psychedelic'
    },
    'bufotenin': {
        'smiles': 'CN(C)CCc1c[nH]c2ccc(O)cc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 204.27,
        'notes': '5-HO-DMT, found in toad venom'
    },
    'DPT': {
        'smiles': 'CCCN(CCC)CCc1c[nH]c2ccccc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT1A'],
        'mw': 244.38,
        'notes': 'Dipropyltryptamine'
    },
    'DiPT': {
        'smiles': 'CC(C)N(C(C)C)CCc1c[nH]c2ccccc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A'],
        'mw': 244.38,
        'notes': 'Diisopropyltryptamine, auditory effects'
    },
    '4-AcO-DMT': {
        'smiles': 'CN(C)CCc1c[nH]c2cccc(OC(C)=O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 246.31,
        'notes': 'Psilacetin, prodrug of psilocin'
    },
    '4-HO-MET': {
        'smiles': 'CCN(C)CCc1c[nH]c2cccc(O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 218.30,
        'notes': '4-Hydroxy-N-methyl-N-ethyltryptamine'
    },
    '4-HO-MiPT': {
        'smiles': 'CC(C)N(C)CCc1c[nH]c2cccc(O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 232.32,
        'notes': '4-Hydroxy-N-methyl-N-isopropyltryptamine'
    },
    '4-AcO-MET': {
        'smiles': 'CCN(C)CCc1c[nH]c2cccc(OC(C)=O)c12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 260.33,
        'notes': '4-Acetoxy-N-methyl-N-ethyltryptamine'
    },
    '5-MeO-MiPT': {
        'smiles': 'CC(C)N(C)CCc1c[nH]c2ccc(OC)cc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT1A'],
        'mw': 246.35,
        'notes': '5-Methoxy-N-methyl-N-isopropyltryptamine'
    },
    '5-MeO-DiPT': {
        'smiles': 'CC(C)N(C(C)C)CCc1c[nH]c2ccc(OC)cc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', '5-HT1A'],
        'mw': 274.40,
        'notes': 'Foxy methoxy'
    },
    'AMT': {
        'smiles': 'NCCc1c[nH]c2ccccc12',
        'class': 'tryptamine',
        'targets': ['5-HT2A', 'SERT', 'MAO'],
        'mw': 174.24,
        'notes': 'Alpha-methyltryptamine, MAOI activity'
    },
    
    # ========================================================================
    # LYSERGAMIDES
    # ========================================================================
    '1P-LSD': {
        'smiles': 'CCCN(C=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(C2=C1)c34)C',
        'class': 'lysergamide',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C'],
        'mw': 379.50,
        'notes': '1-Propionyl-LSD, prodrug'
    },
    'AL-LAD': {
        'smiles': 'CC=CCN(CC=C)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(C2=C1)c34)C',
        'class': 'lysergamide',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 349.47,
        'notes': '6-Allyl-6-nor-LSD'
    },
    'ETH-LAD': {
        'smiles': 'CCN(CC)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(C2=C1)c34)CC',
        'class': 'lysergamide',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 351.49,
        'notes': '6-Ethyl-6-nor-LSD'
    },
    'LSA': {
        'smiles': 'NC(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(C2=C1)c34)C',
        'class': 'lysergamide',
        'targets': ['5-HT2A', 'D2'],
        'mw': 267.33,
        'notes': 'D-lysergic acid amide, from morning glory'
    },
    'ALD-52': {
        'smiles': 'CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c(C(C)=O)[nH]c4cccc(C2=C1)c34)C',
        'class': 'lysergamide',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 365.47,
        'notes': '1-Acetyl-LSD'
    },
    
    # ========================================================================
    # PHENETHYLAMINES - Mescaline derivatives
    # ========================================================================
    'mescaline': {
        'smiles': 'COc1cc(CCN)cc(OC)c1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C'],
        'mw': 211.26,
        'notes': 'From peyote cactus'
    },
    'escaline': {
        'smiles': 'CCOc1cc(CCN)cc(OC)c1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 225.29,
        'notes': '4-Ethoxy analog of mescaline'
    },
    'proscaline': {
        'smiles': 'CCCOc1cc(CCN)cc(OC)c1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 239.31,
        'notes': '4-Propoxy analog of mescaline'
    },
    
    # ========================================================================
    # PHENETHYLAMINES - 2C-x series
    # ========================================================================
    '2C-B': {
        'smiles': 'COc1cc(CCN)c(Br)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2B', '5-HT2C'],
        'mw': 260.11,
        'notes': '4-Bromo-2,5-dimethoxyphenethylamine'
    },
    '2C-E': {
        'smiles': 'CCc1cc(OC)c(CCN)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 209.29,
        'notes': '4-Ethyl-2,5-dimethoxyphenethylamine'
    },
    '2C-I': {
        'smiles': 'COc1cc(CCN)c(I)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 307.13,
        'notes': '4-Iodo-2,5-dimethoxyphenethylamine'
    },
    '2C-C': {
        'smiles': 'COc1cc(CCN)c(Cl)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 215.68,
        'notes': '4-Chloro-2,5-dimethoxyphenethylamine'
    },
    '2C-P': {
        'smiles': 'CCCc1cc(OC)c(CCN)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 223.31,
        'notes': '4-Propyl-2,5-dimethoxyphenethylamine, potent'
    },
    '2C-T-2': {
        'smiles': 'CCSc1cc(OC)c(CCN)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 241.35,
        'notes': '4-Ethylthio-2,5-dimethoxyphenethylamine'
    },
    '2C-T-7': {
        'smiles': 'CCCCSc1cc(OC)c(CCN)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 269.40,
        'notes': '4-Propylthio-2,5-dimethoxyphenethylamine'
    },
    '2C-D': {
        'smiles': 'Cc1cc(OC)c(CCN)cc1OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 195.26,
        'notes': '4-Methyl-2,5-dimethoxyphenethylamine'
    },
    
    # ========================================================================
    # PHENETHYLAMINES - DOx series (amphetamines)
    # ========================================================================
    'DOC': {
        'smiles': 'CC(N)Cc1cc(OC)c(Cl)cc1OC',
        'class': 'amphetamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 229.70,
        'notes': '4-Chloro-2,5-dimethoxyamphetamine'
    },
    'DOM': {
        'smiles': 'CC(N)Cc1cc(OC)c(C)cc1OC',
        'class': 'amphetamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 209.29,
        'notes': '4-Methyl-2,5-dimethoxyamphetamine, STP'
    },
    'DOI': {
        'smiles': 'CC(N)Cc1cc(OC)c(I)cc1OC',
        'class': 'amphetamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 321.16,
        'notes': '4-Iodo-2,5-dimethoxyamphetamine, research tool'
    },
    'DOB': {
        'smiles': 'CC(N)Cc1cc(OC)c(Br)cc1OC',
        'class': 'amphetamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 274.15,
        'notes': '4-Bromo-2,5-dimethoxyamphetamine'
    },
    
    # ========================================================================
    # PHENETHYLAMINES - NBOMe series
    # ========================================================================
    '25I-NBOMe': {
        'smiles': 'COc1cc(I)c(OC)cc1CCNCc2ccccc2OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 427.28,
        'notes': 'N-benzyl analog, high 5-HT2A affinity'
    },
    '25C-NBOMe': {
        'smiles': 'COc1cc(Cl)c(OC)cc1CCNCc2ccccc2OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 335.83,
        'notes': 'N-benzyl-2C-C analog'
    },
    '25B-NBOMe': {
        'smiles': 'COc1cc(Br)c(OC)cc1CCNCc2ccccc2OC',
        'class': 'phenethylamine',
        'targets': ['5-HT2A', '5-HT2C'],
        'mw': 380.27,
        'notes': 'N-benzyl-2C-B analog'
    },
    
    # ========================================================================
    # DISSOCIATIVES - Arylcyclohexylamines
    # ========================================================================
    'ketamine': {
        'smiles': 'CNC1(c2ccccc2Cl)CCCCC1=O',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B', 'Sigma1', 'MOR', 'D2'],
        'mw': 237.73,
        'notes': 'NMDA antagonist, anesthetic, antidepressant'
    },
    'PCP': {
        'smiles': 'c1ccc(C2(N3CCCCC3)CCCCC2)cc1',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B', 'Sigma1', 'D2'],
        'mw': 243.39,
        'notes': 'Phencyclidine, potent dissociative'
    },
    'MXE': {
        'smiles': 'CCNC1(c2ccccc2OC)CCCCC1=O',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B', 'Sigma1'],
        'mw': 247.33,
        'notes': 'Methoxetamine'
    },
    '3-MeO-PCP': {
        'smiles': 'COc1cccc(C2(N3CCCCC3)CCCCC2)c1',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'Sigma1', 'MOR'],
        'mw': 273.42,
        'notes': '3-Methoxy-PCP, sigma agonist'
    },
    '3-HO-PCP': {
        'smiles': 'Oc1cccc(C2(N3CCCCC3)CCCCC2)c1',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'MOR', 'Sigma1'],
        'mw': 259.39,
        'notes': '3-Hydroxy-PCP, opioid activity'
    },
    '3-MeO-PCE': {
        'smiles': 'CCNC1(c2cccc(OC)c2)CCCCC1',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'Sigma1'],
        'mw': 247.38,
        'notes': '3-Methoxy-PCE'
    },
    'DCK': {
        'smiles': 'CNC1(c2ccc(Cl)cc2Cl)CCCCC1=O',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B'],
        'mw': 272.17,
        'notes': 'Deschloroketamine (2-Cl-ketamine)'
    },
    '2-FDCK': {
        'smiles': 'CNC1(c2ccccc2F)CCCCC1=O',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B'],
        'mw': 221.30,
        'notes': '2-Fluorodeschloroketamine'
    },
    'O-PCE': {
        'smiles': 'CCNC1(c2ccccc2)CCCCC1=O',
        'class': 'arylcyclohexylamine',
        'targets': ['NMDA_GluN2A', 'Sigma1'],
        'mw': 217.31,
        'notes': '2-Oxo-PCE'
    },
    
    # ========================================================================
    # DISSOCIATIVES - Other
    # ========================================================================
    'DXM': {
        'smiles': 'COc1ccc2c(c1)[C@@H]1[C@@H](CC[C@H]3[C@@H]1C(C)(C)N(C)CC3)O2',
        'class': 'morphinan',
        'targets': ['NMDA_GluN2A', 'Sigma1', 'SERT'],
        'mw': 271.40,
        'notes': 'Dextromethorphan, OTC dissociative'
    },
    'nitrous_oxide': {
        'smiles': '[N-]=[N+]=O',
        'class': 'gas',
        'targets': ['NMDA_GluN2A', 'GABA_A_alpha1'],
        'mw': 44.01,
        'notes': 'N2O, laughing gas'
    },
    'memantine': {
        'smiles': 'CC12CC3CC(C)(C1)CC(N)(C3)C2',
        'class': 'adamantane',
        'targets': ['NMDA_GluN2A', 'Sigma1'],
        'mw': 179.30,
        'notes': 'Low-affinity NMDA antagonist'
    },
    
    # ========================================================================
    # EMPATHOGENS / ENTACTOGENS
    # ========================================================================
    'MDMA': {
        'smiles': 'CC(N)Cc1ccc2c(c1)OCO2',
        'class': 'amphetamine',
        'targets': ['SERT', '5-HT2A', '5-HT2B', 'D2'],
        'mw': 193.24,
        'notes': '3,4-Methylenedioxymethamphetamine, ecstasy'
    },
    'MDA': {
        'smiles': 'CC(N)Cc1ccc2c(c1)OCO2',
        'class': 'amphetamine',
        'targets': ['SERT', '5-HT2A', '5-HT2B'],
        'mw': 179.22,
        'notes': '3,4-Methylenedioxyamphetamine'
    },
    '6-APB': {
        'smiles': 'CC(N)Cc1ccc2occc2c1',
        'class': 'benzofuran',
        'targets': ['SERT', '5-HT2A', '5-HT2B'],
        'mw': 175.23,
        'notes': '6-(2-Aminopropyl)benzofuran'
    },
    '5-MAPB': {
        'smiles': 'CNC(C)Cc1ccc2occc2c1',
        'class': 'benzofuran',
        'targets': ['SERT', '5-HT2B'],
        'mw': 189.25,
        'notes': '5-(2-Methylaminopropyl)benzofuran'
    },
    'methylone': {
        'smiles': 'CNC(C)C(=O)c1ccc2c(c1)OCO2',
        'class': 'cathinone',
        'targets': ['SERT', 'DAT', 'NET'],
        'mw': 207.23,
        'notes': 'bk-MDMA, synthetic cathinone'
    },
    
    # ========================================================================
    # CANNABINOIDS
    # ========================================================================
    'THC': {
        'smiles': 'CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1',
        'class': 'cannabinoid',
        'targets': ['CB1'],
        'mw': 314.46,
        'notes': 'Delta-9-THC, primary psychoactive in cannabis'
    },
    'CBD': {
        'smiles': 'CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1',
        'class': 'cannabinoid',
        'targets': ['CB1', '5-HT1A'],
        'mw': 314.46,
        'notes': 'Cannabidiol, non-psychoactive'
    },
    'delta8_THC': {
        'smiles': 'CCCCCc1cc(O)c2C3CC(C)=CCC3C(C)(C)Oc2c1',
        'class': 'cannabinoid',
        'targets': ['CB1'],
        'mw': 314.46,
        'notes': 'Delta-8-THC isomer'
    },
    
    # ========================================================================
    # KAPPA OPIOID AGONISTS
    # ========================================================================
    'salvinorin_A': {
        'smiles': 'COC(=O)C1CC2C(OC(C)=O)C(=O)O[C@@H]3C[C@H](C)[C@@H]4CC(=O)[C@H](C1)C234',
        'class': 'terpenoid',
        'targets': ['KOR'],
        'mw': 432.46,
        'notes': 'From Salvia divinorum, selective KOR agonist'
    },
    
    # ========================================================================
    # DELIRIANTS
    # ========================================================================
    'scopolamine': {
        'smiles': 'CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C1OC21',
        'class': 'tropane',
        'targets': ['M1_muscarinic'],
        'mw': 303.35,
        'notes': 'Muscarinic antagonist, deliriant'
    },
    'atropine': {
        'smiles': 'CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C(C2)C',
        'class': 'tropane',
        'targets': ['M1_muscarinic'],
        'mw': 289.37,
        'notes': 'Muscarinic antagonist'
    },
    'diphenhydramine': {
        'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2',
        'class': 'antihistamine',
        'targets': ['M1_muscarinic', 'H1'],
        'mw': 255.36,
        'notes': 'DPH, deliriant at high doses'
    },
    
    # ========================================================================
    # STIMULANTS
    # ========================================================================
    'amphetamine': {
        'smiles': 'CC(N)Cc1ccccc1',
        'class': 'amphetamine',
        'targets': ['DAT', 'NET', 'D1', 'D2'],
        'mw': 135.21,
        'notes': 'Releases dopamine and norepinephrine'
    },
    'methamphetamine': {
        'smiles': 'CNC(C)Cc1ccccc1',
        'class': 'amphetamine',
        'targets': ['DAT', 'NET', 'D1', 'D2'],
        'mw': 149.23,
        'notes': 'More potent than amphetamine'
    },
    'cocaine': {
        'smiles': 'COC(=O)C1C(OC(=O)c2ccccc2)CC2CCC1N2C',
        'class': 'tropane',
        'targets': ['DAT', 'NET', 'SERT'],
        'mw': 303.35,
        'notes': 'DAT/NET/SERT reuptake inhibitor'
    },
    'caffeine': {
        'smiles': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
        'class': 'xanthine',
        'targets': ['A1', 'A2A'],
        'mw': 194.19,
        'notes': 'Adenosine receptor antagonist'
    },
    'nicotine': {
        'smiles': 'CN1CCCC1c2cccnc2',
        'class': 'pyridine',
        'targets': ['nAChR_alpha4beta2'],
        'mw': 162.23,
        'notes': 'Nicotinic agonist'
    },
    'modafinil': {
        'smiles': 'NC(=O)CS(=O)C(c1ccccc1)c2ccccc2',
        'class': 'benzhydryl',
        'targets': ['DAT', 'D2'],
        'mw': 273.35,
        'notes': 'Eugeroic, weak DAT inhibitor'
    },
    
    # ========================================================================
    # OPIOIDS
    # ========================================================================
    'morphine': {
        'smiles': 'CN1CCC23C4Oc5c(O)ccc(C12CC(O)C=CC34)c5',
        'class': 'opioid',
        'targets': ['MOR', 'KOR', 'DOR'],
        'mw': 285.34,
        'notes': 'Prototypical MOR agonist'
    },
    'heroin': {
        'smiles': 'CC(=O)Oc1ccc2C3CC=CC4C(Oc1c24)C(OC(C)=O)CC3N(C)C',
        'class': 'opioid',
        'targets': ['MOR', 'KOR', 'DOR'],
        'mw': 369.41,
        'notes': 'Diacetylmorphine, prodrug'
    },
    'oxycodone': {
        'smiles': 'CN1CCC23C4Oc5c(OC)ccc(C12C(=O)C(O)CC34)c5',
        'class': 'opioid',
        'targets': ['MOR', 'KOR'],
        'mw': 315.36,
        'notes': 'Semi-synthetic opioid'
    },
    'fentanyl': {
        'smiles': 'CCC(=O)N(c1ccccc1)C2CCN(CCc3ccccc3)CC2',
        'class': 'opioid',
        'targets': ['MOR'],
        'mw': 336.47,
        'notes': 'Highly potent synthetic opioid'
    },
    'kratom_mitragynine': {
        'smiles': 'CCC1CN2CCc3c([nH]c4cccc(OC)c34)C2CC1C(=O)OC',
        'class': 'indole',
        'targets': ['MOR', 'KOR', 'DOR', 'alpha2A'],
        'mw': 398.50,
        'notes': 'Main alkaloid from kratom'
    },
    'buprenorphine': {
        'smiles': 'COC1=CC2=C(C=C1)[C@]34CC[C@@H]5[C@H](CC[C@@]6([C@@H]5[C@@H]3O)CCN(C)CC7=CC=CC=C7)[C@@]4(C(C)(C)C)OC6=O',
        'class': 'opioid',
        'targets': ['MOR', 'KOR', 'DOR'],
        'mw': 467.64,
        'notes': 'Partial MOR agonist, KOR antagonist'
    },
    
    # ========================================================================
    # GABAERGICS
    # ========================================================================
    'GHB': {
        'smiles': 'OCCCC(=O)O',
        'class': 'fatty_acid',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2', 'GHB_receptor'],
        'mw': 104.10,
        'notes': 'Gamma-hydroxybutyrate'
    },
    'phenibut': {
        'smiles': 'NCCCc1ccc(Cl)cc1',
        'class': 'gaba_analog',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 179.22,
        'notes': 'GABA-B agonist, anxiolytic'
    },
    'alprazolam': {
        'smiles': 'Cc1nnc2CN=C(c3ccccc3Cl)c4cc(ccc4n12)Cl',
        'class': 'benzodiazepine',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 308.76,
        'notes': 'Xanax, high-potency benzo'
    },
    'diazepam': {
        'smiles': 'CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc31',
        'class': 'benzodiazepine',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 284.74,
        'notes': 'Valium, prototypical benzo'
    },
    'clonazepam': {
        'smiles': 'O=[N+]([O-])c1ccc2NC(=O)CN=C(c3ccccc3Cl)c2c1',
        'class': 'benzodiazepine',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 315.71,
        'notes': 'Klonopin, anticonvulsant'
    },
    'lorazepam': {
        'smiles': 'OC1N=C(c2ccccc2Cl)c3cc(Cl)ccc3NC1=O',
        'class': 'benzodiazepine',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 321.16,
        'notes': 'Ativan'
    },
    'etizolam': {
        'smiles': 'CCc1cc2c(s1)NC(=O)CN=C(c3ccccc3Cl)c3cccnc32',
        'class': 'thienodiazepine',
        'targets': ['GABA_A_alpha1'],
        'mw': 342.85,
        'notes': 'Thienodiazepine, not a true benzo'
    },
    'zolpidem': {
        'smiles': 'Cc1ccc(nc1)Cc2c(C)nc3ccc(C)cn23',
        'class': 'imidazopyridine',
        'targets': ['GABA_A_alpha1'],
        'mw': 307.39,
        'notes': 'Ambien, GABA-A α1 selective'
    },
    
    # ========================================================================
    # IBOGA ALKALOIDS
    # ========================================================================
    'ibogaine': {
        'smiles': 'CCC1CN2CCc3c([nH]c4cccc(OC)c34)C2CC1C(=O)OC',
        'class': 'indole',
        'targets': ['NMDA_GluN2A', 'Sigma1', 'KOR', 'SERT'],
        'mw': 310.43,
        'notes': 'From Tabernanthe iboga, anti-addiction'
    },
    
    # ========================================================================
    # AMANITA / GABAergic psychedelics
    # ========================================================================
    'muscimol': {
        'smiles': 'OC1=CC(=NO1)CN',
        'class': 'isoxazole',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 114.10,
        'notes': 'From Amanita muscaria, GABA-A agonist'
    },
    'ibotenic_acid': {
        'smiles': 'NC(Cc1cc(=O)[nH]o1)C(=O)O',
        'class': 'isoxazole',
        'targets': ['NMDA_GluN2A', 'AMPA_GluA1'],
        'mw': 158.11,
        'notes': 'Excitotoxin, glutamate agonist'
    },
    
    # ========================================================================
    # NEUROTRANSMITTERS (for reference/control)
    # ========================================================================
    'serotonin': {
        'smiles': 'NCCc1c[nH]c2ccc(O)cc12',
        'class': 'neurotransmitter',
        'targets': ['5-HT1A', '5-HT2A', '5-HT2B', '5-HT2C'],
        'mw': 176.22,
        'notes': '5-HT, endogenous'
    },
    'dopamine': {
        'smiles': 'NCCc1ccc(O)c(O)c1',
        'class': 'neurotransmitter',
        'targets': ['D1', 'D2', 'D3'],
        'mw': 153.18,
        'notes': 'Endogenous catecholamine'
    },
    'norepinephrine': {
        'smiles': 'NCC(O)c1ccc(O)c(O)c1',
        'class': 'neurotransmitter',
        'targets': ['alpha2A'],
        'mw': 169.18,
        'notes': 'Endogenous catecholamine'
    },
    'GABA': {
        'smiles': 'NCCCC(=O)O',
        'class': 'neurotransmitter',
        'targets': ['GABA_A_alpha1', 'GABA_A_alpha2'],
        'mw': 103.12,
        'notes': 'Gamma-aminobutyric acid'
    },
    'glutamate': {
        'smiles': 'NC(CCC(=O)O)C(=O)O',
        'class': 'neurotransmitter',
        'targets': ['NMDA_GluN2A', 'NMDA_GluN2B', 'AMPA_GluA1'],
        'mw': 147.13,
        'notes': 'Primary excitatory neurotransmitter'
    },
    'acetylcholine': {
        'smiles': 'CC(=O)OCC[N+](C)(C)C',
        'class': 'neurotransmitter',
        'targets': ['M1_muscarinic', 'nAChR_alpha4beta2'],
        'mw': 146.21,
        'notes': 'Endogenous cholinergic ligand'
    },
}


def get_ligand_smiles(name):
    """Get SMILES for a ligand by name."""
    ligand = LIGAND_DATABASE.get(name)
    return ligand['smiles'] if ligand else None


def get_ligands_for_receptor(receptor_name):
    """Get all ligands that target a specific receptor."""
    return {
        name: data for name, data in LIGAND_DATABASE.items()
        if receptor_name in data.get('targets', [])
    }


def get_all_ligand_names():
    """Get list of all ligand names."""
    return list(LIGAND_DATABASE.keys())


def print_database_summary():
    """Print summary of the ligand database."""
    print("=" * 80)
    print("PSYCHOACTIVE LIGAND DATABASE")
    print("=" * 80)
    
    # Count by class
    classes = {}
    for name, data in LIGAND_DATABASE.items():
        cls = data.get('class', 'unknown')
        classes[cls] = classes.get(cls, 0) + 1
    
    print("\nLigands by class:")
    for cls, count in sorted(classes.items(), key=lambda x: -x[1]):
        print(f"  {cls:25s}: {count}")
    
    print(f"\nTotal ligands: {len(LIGAND_DATABASE)}")
    
    # Show receptor coverage
    print("\n" + "-" * 40)
    print("Ligands per receptor target:")
    
    receptor_counts = {}
    for name, data in LIGAND_DATABASE.items():
        for target in data.get('targets', []):
            receptor_counts[target] = receptor_counts.get(target, 0) + 1
    
    for rec, count in sorted(receptor_counts.items(), key=lambda x: -x[1]):
        print(f"  {rec:20s}: {count}")


if __name__ == '__main__':
    print_database_summary()
