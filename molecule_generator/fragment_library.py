#!/usr/bin/env python3
"""
Advanced molecular fragment library for drug generation.

Contains validated fragments known to affect receptor binding,
organized by receptor affinity and drug class.
"""

from dataclasses import dataclass
from typing import Dict, List, Set


@dataclass
class Fragment:
    """A molecular fragment with metadata."""
    smiles: str
    name: str
    description: str
    receptor_effects: Dict[str, str]  # receptor -> "enhance" or "reduce"
    drug_classes: List[str]
    attachment_points: int  # Number of positions where this can connect


# Core scaffolds - the backbone structures
CORE_SCAFFOLDS = {
    # Indole family (tryptamines, LSD)
    "indole": Fragment(
        smiles="c1ccc2[nH]ccc2c1",
        name="Indole",
        description="Core structure of tryptamines and ergolines",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2B": "enhance", 
            "5-HT2C": "enhance",
            "5-HT1A": "enhance"
        },
        drug_classes=["psychedelic", "tryptamine"],
        attachment_points=4
    ),
    
    "tryptamine": Fragment(
        smiles="NCCc1c[nH]c2ccccc12",
        name="Tryptamine",
        description="Indole with ethylamine side chain",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2B": "enhance",
            "5-HT1A": "enhance"
        },
        drug_classes=["psychedelic", "tryptamine"],
        attachment_points=3
    ),
    
    # Phenethylamine family (amphetamines, mescaline, 2C-x)
    "phenethylamine": Fragment(
        smiles="NCCc1ccccc1",
        name="Phenethylamine",
        description="Base structure of many stimulants and psychedelics",
        receptor_effects={
            "D2": "enhance",
            "NET": "enhance",
            "DAT": "enhance"
        },
        drug_classes=["stimulant", "phenethylamine"],
        attachment_points=5
    ),
    
    "amphetamine": Fragment(
        smiles="CC(N)Cc1ccccc1",
        name="Amphetamine",
        description="Alpha-methylated phenethylamine",
        receptor_effects={
            "DAT": "enhance",
            "NET": "enhance",
            "D2": "enhance",
            "TAAR1": "enhance"
        },
        drug_classes=["stimulant"],
        attachment_points=4
    ),
    
    # Piperidine/morpholine (many drugs)
    "piperidine": Fragment(
        smiles="C1CCNCC1",
        name="Piperidine",
        description="Six-membered nitrogen heterocycle",
        receptor_effects={
            "MOR": "enhance",
            "D2": "enhance"
        },
        drug_classes=["opioid", "stimulant"],
        attachment_points=2
    ),
    
    "piperazine": Fragment(
        smiles="C1CNCCN1",
        name="Piperazine",
        description="Two nitrogen heterocycle",
        receptor_effects={
            "5-HT2A": "neutral",
            "D2": "enhance"
        },
        drug_classes=["stimulant"],
        attachment_points=2
    ),
    
    "morpholine": Fragment(
        smiles="C1COCCN1",
        name="Morpholine",
        description="Oxygen-containing heterocycle",
        receptor_effects={
            "MOR": "enhance"
        },
        drug_classes=["opioid"],
        attachment_points=2
    ),
    
    # Arylcyclohexylamine (ketamine, PCP)
    "arylcyclohexylamine": Fragment(
        smiles="NC1CCCCC1c1ccccc1",
        name="Arylcyclohexylamine",
        description="Dissociative scaffold",
        receptor_effects={
            "NMDA": "enhance",
            "Sigma1": "enhance",
            "D2": "enhance"
        },
        drug_classes=["dissociative"],
        attachment_points=3
    ),
    
    # Benzofuran (6-APB, etc)
    "benzofuran": Fragment(
        smiles="c1ccc2occc2c1",
        name="Benzofuran",
        description="Fused furan-benzene ring",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2B": "enhance"
        },
        drug_classes=["entactogen"],
        attachment_points=4
    ),
    
    # Cathinone
    "cathinone": Fragment(
        smiles="CC(N)C(=O)c1ccccc1",
        name="Cathinone",
        description="Beta-keto amphetamine",
        receptor_effects={
            "DAT": "enhance",
            "NET": "enhance",
            "SERT": "enhance"
        },
        drug_classes=["stimulant"],
        attachment_points=4
    ),
}


# Substituent fragments that modify activity
SUBSTITUENTS = {
    # Methyl groups
    "n_methyl": Fragment(
        smiles="CN",
        name="N-Methyl",
        description="Methylation of amine",
        receptor_effects={
            "5-HT2A": "neutral",
            "DAT": "enhance"
        },
        drug_classes=[],
        attachment_points=1
    ),
    
    "n_dimethyl": Fragment(
        smiles="CN(C)",
        name="N,N-Dimethyl",
        description="Dimethylation of amine (like DMT)",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT1A": "enhance"
        },
        drug_classes=["tryptamine"],
        attachment_points=1
    ),
    
    # Methoxy groups
    "4_methoxy": Fragment(
        smiles="COc1ccc(*)cc1",
        name="4-Methoxy",
        description="Para-methoxyphenyl",
        receptor_effects={
            "5-HT2A": "enhance"
        },
        drug_classes=["psychedelic"],
        attachment_points=1
    ),
    
    "3_4_methylenedioxy": Fragment(
        smiles="c1cc2OCOc2cc1",
        name="3,4-Methylenedioxy",
        description="MDMA-like substitution pattern",
        receptor_effects={
            "SERT": "enhance",
            "5-HT2B": "enhance"
        },
        drug_classes=["entactogen"],
        attachment_points=2
    ),
    
    "2_5_dimethoxy": Fragment(
        smiles="COc1ccc(OC)cc1",
        name="2,5-Dimethoxy",
        description="2C-x substitution pattern",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2C": "enhance"
        },
        drug_classes=["psychedelic"],
        attachment_points=2
    ),
    
    "3_4_5_trimethoxy": Fragment(
        smiles="COc1cc(*)cc(OC)c1OC",
        name="3,4,5-Trimethoxy",
        description="Mescaline substitution pattern",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2C": "enhance"
        },
        drug_classes=["psychedelic"],
        attachment_points=1
    ),
    
    # Hydroxy groups
    "4_hydroxy": Fragment(
        smiles="Oc1ccc(*)cc1",
        name="4-Hydroxy",
        description="Para-hydroxyphenyl (psilocin pattern)",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT1A": "enhance"
        },
        drug_classes=["tryptamine"],
        attachment_points=1
    ),
    
    "5_hydroxy": Fragment(
        smiles="*c1ccc(O)cc1",
        name="5-Hydroxy",
        description="5-HT (serotonin) pattern",
        receptor_effects={
            "SERT": "enhance"
        },
        drug_classes=[],
        attachment_points=1
    ),
    
    # Halogens
    "4_fluoro": Fragment(
        smiles="Fc1ccc(*)cc1",
        name="4-Fluoro",
        description="Para-fluorophenyl",
        receptor_effects={
            "D2": "enhance",
            "DAT": "enhance"
        },
        drug_classes=["stimulant"],
        attachment_points=1
    ),
    
    "4_chloro": Fragment(
        smiles="Clc1ccc(*)cc1",
        name="4-Chloro",
        description="Para-chlorophenyl",
        receptor_effects={
            "D2": "enhance"
        },
        drug_classes=[],
        attachment_points=1
    ),
    
    "4_bromo": Fragment(
        smiles="Brc1ccc(*)cc1",
        name="4-Bromo",
        description="Para-bromophenyl (2C-B pattern)",
        receptor_effects={
            "5-HT2A": "enhance",
            "5-HT2C": "enhance"
        },
        drug_classes=["psychedelic"],
        attachment_points=1
    ),
    
    "4_iodo": Fragment(
        smiles="Ic1ccc(*)cc1",
        name="4-Iodo",
        description="Para-iodophenyl (DOI pattern)",
        receptor_effects={
            "5-HT2A": "enhance"
        },
        drug_classes=["psychedelic"],
        attachment_points=1
    ),
    
    # Acetyl/carbonyl
    "acetyl": Fragment(
        smiles="CC(=O)*",
        name="Acetyl",
        description="Acetyl group",
        receptor_effects={},
        drug_classes=[],
        attachment_points=1
    ),
    
    # Phosphate (prodrug)
    "phosphate": Fragment(
        smiles="OP(=O)(O)O*",
        name="Phosphate",
        description="Phosphate ester (psilocybin pattern)",
        receptor_effects={},
        drug_classes=["prodrug"],
        attachment_points=1
    ),
}


# Linker fragments
LINKERS = {
    "ethyl": Fragment(
        smiles="*CC*",
        name="Ethyl linker",
        description="Two-carbon linker",
        receptor_effects={},
        drug_classes=[],
        attachment_points=2
    ),
    
    "propyl": Fragment(
        smiles="*CCC*",
        name="Propyl linker",
        description="Three-carbon linker",
        receptor_effects={},
        drug_classes=[],
        attachment_points=2
    ),
    
    "alpha_methyl": Fragment(
        smiles="*C(C)*",
        name="Alpha-methyl",
        description="Branched carbon (amphetamine pattern)",
        receptor_effects={
            "DAT": "enhance"
        },
        drug_classes=["stimulant"],
        attachment_points=2
    ),
}


def get_fragments_for_receptor(receptor: str, effect: str = "enhance") -> List[Fragment]:
    """Get all fragments that have the specified effect on a receptor."""
    fragments = []
    
    for frag in list(CORE_SCAFFOLDS.values()) + list(SUBSTITUENTS.values()):
        if receptor in frag.receptor_effects:
            if frag.receptor_effects[receptor] == effect:
                fragments.append(frag)
    
    return fragments


def get_fragments_for_class(drug_class: str) -> List[Fragment]:
    """Get all fragments associated with a drug class."""
    fragments = []
    
    for frag in list(CORE_SCAFFOLDS.values()) + list(SUBSTITUENTS.values()):
        if drug_class in frag.drug_classes:
            fragments.append(frag)
    
    return fragments


def suggest_modifications(current_profile: Dict[str, float],
                         target_profile: Dict[str, float]) -> List[str]:
    """
    Suggest fragment modifications to move from current to target profile.
    """
    suggestions = []
    
    for receptor, target_val in target_profile.items():
        current_val = current_profile.get(receptor, 0)
        
        if target_val < current_val:  # Need stronger binding (more negative)
            enhancers = get_fragments_for_receptor(receptor, "enhance")
            for frag in enhancers[:3]:
                suggestions.append(f"Add {frag.name} to enhance {receptor} binding")
        elif target_val > current_val:  # Need weaker binding
            reducers = get_fragments_for_receptor(receptor, "reduce")
            for frag in reducers[:2]:
                suggestions.append(f"Add {frag.name} to reduce {receptor} binding")
    
    return suggestions


# Example drug blueprints showing how fragments combine
DRUG_BLUEPRINTS = {
    "DMT": {
        "core": "tryptamine",
        "substituents": ["n_dimethyl"],
        "profile": {"5-HT2A": -85000, "5-HT1A": -60000}
    },
    "psilocin": {
        "core": "tryptamine",
        "substituents": ["4_hydroxy", "n_dimethyl"],
        "profile": {"5-HT2A": -80000, "5-HT1A": -55000}
    },
    "mescaline": {
        "core": "phenethylamine",
        "substituents": ["3_4_5_trimethoxy"],
        "profile": {"5-HT2A": -75000, "5-HT2C": -65000}
    },
    "MDMA": {
        "core": "amphetamine",
        "substituents": ["3_4_methylenedioxy"],
        "profile": {"SERT": -90000, "5-HT2B": -70000}
    },
    "2C-B": {
        "core": "phenethylamine",
        "substituents": ["2_5_dimethoxy", "4_bromo"],
        "profile": {"5-HT2A": -82000, "5-HT2C": -75000}
    },
    "ketamine": {
        "core": "arylcyclohexylamine",
        "substituents": ["4_chloro"],
        "profile": {"NMDA": -85000, "Sigma1": -70000}
    },
}


if __name__ == "__main__":
    print("Fragment Library Summary")
    print("=" * 50)
    print(f"\nCore scaffolds: {len(CORE_SCAFFOLDS)}")
    for name, frag in CORE_SCAFFOLDS.items():
        print(f"  - {frag.name}: {frag.smiles}")
    
    print(f"\nSubstituents: {len(SUBSTITUENTS)}")
    for name, frag in SUBSTITUENTS.items():
        print(f"  - {frag.name}: {frag.smiles}")
    
    print(f"\nLinkers: {len(LINKERS)}")
    
    print("\n" + "=" * 50)
    print("Fragments that enhance 5-HT2A binding:")
    for frag in get_fragments_for_receptor("5-HT2A", "enhance"):
        print(f"  - {frag.name}")
    
    print("\nFragments that enhance D2 binding:")
    for frag in get_fragments_for_receptor("D2", "enhance"):
        print(f"  - {frag.name}")
