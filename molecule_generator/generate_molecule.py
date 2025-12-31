#!/usr/bin/env python3
"""
Molecule Generator from Target Binding Profile

Given a desired receptor binding profile, generate novel valid molecules
that are predicted to have similar properties.

Supports two modes:
1. Fast mode: Predict binding from structural similarity (default)
2. Full docking mode: Run Chai-1 + OpenMM pipeline (--dock flag)
"""

import json
import os
import sys
import random
import argparse
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from collections import defaultdict
from datetime import datetime

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.rdchem import BondType
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Install with: pip install rdkit")

# Try to import SELFIES for guaranteed-valid mutations
try:
    import selfies as sf
    SELFIES_AVAILABLE = True
except ImportError:
    SELFIES_AVAILABLE = False
    print("Warning: SELFIES not available. Install with: pip install selfies")

import numpy as np

# Path to data
DATA_DIR = Path(__file__).parent.parent / "website"
PROFILE_MATCHER_DIR = Path(__file__).parent.parent / "profile_matcher"
CACHE_DIR = Path(__file__).parent / "cache"


@dataclass
class Molecule:
    """Represents a generated or existing molecule."""
    smiles: str
    name: str = ""
    predicted_profile: Dict[str, float] = field(default_factory=dict)
    docked_profile: Dict[str, float] = field(default_factory=dict)  # From full docking
    source_molecules: List[str] = field(default_factory=list)
    score: float = 0.0
    qed: float = 0.0
    lipinski_violations: int = 0
    docking_completed: bool = False


# ============================================================================
# CACHING SYSTEM
# ============================================================================

class ProfileCache:
    """
    Cache for binding profile predictions and docking results.
    Persists to disk in JSON format.
    """
    
    def __init__(self, cache_dir: Path = None):
        self.cache_dir = cache_dir or CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.prediction_cache: Dict[str, Dict] = {}
        self.docking_cache: Dict[str, Dict] = {}
        self._load_caches()
    
    def _load_caches(self):
        """Load caches from disk."""
        pred_file = self.cache_dir / "prediction_cache.json"
        dock_file = self.cache_dir / "docking_cache.json"
        
        if pred_file.exists():
            try:
                with open(pred_file, 'r') as f:
                    self.prediction_cache = json.load(f)
            except Exception:
                self.prediction_cache = {}
        
        if dock_file.exists():
            try:
                with open(dock_file, 'r') as f:
                    self.docking_cache = json.load(f)
            except Exception:
                self.docking_cache = {}
    
    def _save_caches(self):
        """Save caches to disk."""
        pred_file = self.cache_dir / "prediction_cache.json"
        dock_file = self.cache_dir / "docking_cache.json"
        
        with open(pred_file, 'w') as f:
            json.dump(self.prediction_cache, f)
        
        with open(dock_file, 'w') as f:
            json.dump(self.docking_cache, f)
    
    def _get_key(self, smiles: str) -> str:
        """Get cache key for a SMILES string."""
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        return hashlib.md5(smiles.encode()).hexdigest()
    
    def get_prediction(self, smiles: str) -> Optional[Dict[str, float]]:
        """Get cached prediction for SMILES."""
        key = self._get_key(smiles)
        return self.prediction_cache.get(key)
    
    def put_prediction(self, smiles: str, profile: Dict[str, float]):
        """Cache a prediction."""
        key = self._get_key(smiles)
        self.prediction_cache[key] = profile
        self._save_caches()
    
    def get_docking(self, smiles: str, receptor: str) -> Optional[float]:
        """Get cached docking result."""
        key = self._get_key(smiles)
        if key in self.docking_cache:
            return self.docking_cache[key].get(receptor)
        return None
    
    def put_docking(self, smiles: str, receptor: str, energy: float):
        """Cache a docking result."""
        key = self._get_key(smiles)
        if key not in self.docking_cache:
            self.docking_cache[key] = {}
        self.docking_cache[key][receptor] = energy
        self._save_caches()
    
    def get_full_docking_profile(self, smiles: str) -> Optional[Dict[str, float]]:
        """Get full cached docking profile if available."""
        key = self._get_key(smiles)
        return self.docking_cache.get(key)
    
    def get_stats(self) -> Dict:
        """Get cache statistics."""
        return {
            'predictions_cached': len(self.prediction_cache),
            'docking_entries': len(self.docking_cache),
            'total_docking_pairs': sum(len(v) for v in self.docking_cache.values())
        }


# Global cache instance
_cache = None

def get_cache() -> ProfileCache:
    """Get or create the global cache instance."""
    global _cache
    if _cache is None:
        _cache = ProfileCache()
    return _cache


def load_drug_database() -> Dict:
    """Load drug data from the website data.js file."""
    data_file = DATA_DIR / "data.js"
    
    if not data_file.exists():
        # Try profile matcher directory
        data_file = PROFILE_MATCHER_DIR / "data.js"
    
    if not data_file.exists():
        print(f"Error: Could not find data.js")
        return {"drugs": {}, "receptors": {}}
    
    with open(data_file, 'r') as f:
        content = f.read()
    
    import re
    
    drugs = {}
    receptors = {}
    
    # Find the drugs section - it starts with "drugs: {" and we need to find matching brace
    drugs_start = content.find('drugs: {')
    if drugs_start == -1:
        drugs_start = content.find('drugs:{')
    
    if drugs_start != -1:
        # Find the opening brace of drugs
        brace_start = content.find('{', drugs_start)
        
        # Now parse drug entries - they look like:
        # "DrugName": { "receptors": {...}, "class": "..." }
        
        # Use a regex to find each drug entry
        # Match: "name": { ... "class": "..." }
        drug_pattern = re.compile(
            r'"([^"]+)":\s*\{\s*"receptors":\s*\{([^}]*)\}[^}]*"class":\s*"([^"]*)"',
            re.DOTALL
        )
        
        for match in drug_pattern.finditer(content[drugs_start:]):
            drug_name = match.group(1)
            receptors_str = match.group(2)
            drug_class = match.group(3)
            
            # Parse receptors
            receptors_dict = {}
            receptor_pattern = re.compile(r'"([^"]+)":\s*(-?\d+\.?\d*)')
            for r_match in receptor_pattern.finditer(receptors_str):
                receptors_dict[r_match.group(1)] = float(r_match.group(2))
            
            drugs[drug_name] = {
                "smiles": "",  # Will be loaded from smiles database
                "class": drug_class,
                "receptors": receptors_dict
            }
    
    return {"drugs": drugs, "receptors": receptors}


def load_smiles_database() -> Dict[str, str]:
    """Load SMILES from the ligand database."""
    smiles_file = Path(__file__).parent.parent / "md_simulations" / "ligand_smiles_database.py"
    
    smiles_db = {}
    
    if smiles_file.exists():
        # Import the module directly
        import importlib.util
        spec = importlib.util.spec_from_file_location("ligand_smiles_database", smiles_file)
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
            
            if hasattr(module, 'LIGAND_DATABASE'):
                for name, data in module.LIGAND_DATABASE.items():
                    if 'smiles' in data:
                        smiles_db[name] = data['smiles']
            
            if hasattr(module, 'LIGAND_SMILES'):
                smiles_db.update(module.LIGAND_SMILES)
                
        except Exception as e:
            print(f"Warning: Could not load SMILES database: {e}")
            # Fallback to regex parsing
            with open(smiles_file, 'r') as f:
                content = f.read()
            
            import re
            # Parse LIGAND_DATABASE entries
            pattern = re.compile(r"'([^']+)':\s*\{[^}]*'smiles':\s*'([^']+)'", re.DOTALL)
            for match in pattern.finditer(content):
                smiles_db[match.group(1)] = match.group(2)
    
    return smiles_db


# Common pharmacophore fragments found in psychoactive compounds
PHARMACOPHORE_FRAGMENTS = {
    # Indole core (tryptamines, LSD)
    "indole": "c1ccc2[nH]ccc2c1",
    "tryptamine": "NCCc1c[nH]c2ccccc12",
    
    # Phenethylamine core (amphetamines, mescaline)
    "phenethylamine": "NCCc1ccccc1",
    "amphetamine": "CC(N)Cc1ccccc1",
    
    # Phenyl ring substitutions
    "4-hydroxy-phenyl": "Oc1ccc(CC)cc1",
    "3,4,5-trimethoxy-phenyl": "COc1cc(CC)cc(OC)c1OC",
    "2,5-dimethoxy-phenyl": "COc1ccc(OC)c(CC)c1",
    
    # Ergoline core (LSD)
    "ergoline_partial": "CN1CC(C=C2)C3=CC=CC4=C3C2C1CN4",
    
    # Piperidine (many drugs)
    "piperidine": "C1CCNCC1",
    "morpholine": "C1COCCN1",
    
    # Ketone/amide groups
    "acetyl": "CC(=O)",
    "methylamino": "CNC",
    "dimethylamino": "CN(C)",
    
    # Hydroxyl substitutions
    "4-hydroxy": "Oc1ccc(*)cc1",
    "5-hydroxy": "*c1ccc(O)cc1",
    
    # Halogen substitutions
    "4-fluoro": "Fc1ccc(*)cc1",
    "4-chloro": "Clc1ccc(*)cc1",
    "4-bromo": "Brc1ccc(*)cc1",
}


# Fragment modifications that tend to affect receptor binding
RECEPTOR_MODIFYING_FRAGMENTS = {
    "5-HT2A": {
        "enhance": ["indole", "tryptamine", "2,5-dimethoxy-phenyl", "4-bromo"],
        "reduce": ["piperidine", "morpholine"],
    },
    "D2": {
        "enhance": ["phenethylamine", "piperidine", "4-fluoro"],
        "reduce": ["indole", "3,4,5-trimethoxy-phenyl"],
    },
    "5-HT2B": {
        "enhance": ["indole", "tryptamine"],
        "reduce": ["amphetamine"],
    },
    "MOR": {
        "enhance": ["piperidine", "morpholine", "4-hydroxy"],
        "reduce": ["indole"],
    },
    "CB1": {
        "enhance": ["long_alkyl_chain", "4-hydroxy"],
        "reduce": ["indole", "tryptamine"],
    },
}


def calculate_profile_similarity(profile1: Dict[str, float], 
                                  profile2: Dict[str, float]) -> float:
    """Calculate cosine similarity between two binding profiles."""
    # Get common receptors
    common_receptors = set(profile1.keys()) & set(profile2.keys())
    
    if not common_receptors:
        return 0.0
    
    vec1 = np.array([profile1[r] for r in common_receptors])
    vec2 = np.array([profile2[r] for r in common_receptors])
    
    # Normalize
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    return np.dot(vec1, vec2) / (norm1 * norm2)


def find_similar_drugs(target_profile: Dict[str, float], 
                       drug_db: Dict,
                       smiles_db: Dict[str, str] = None,
                       top_k: int = 10) -> List[Tuple[str, float, Dict]]:
    """Find drugs with similar binding profiles to the target."""
    similarities = []
    
    for drug_name, drug_data in drug_db.get("drugs", {}).items():
        if "receptors" in drug_data and drug_data["receptors"]:
            # Check if we have SMILES for this drug
            has_smiles = False
            if smiles_db:
                # Try exact match and normalized match
                if drug_name in smiles_db:
                    has_smiles = True
                else:
                    norm_name = drug_name.lower().replace('-', '').replace(' ', '')
                    for sname in smiles_db:
                        if sname.lower().replace('-', '').replace(' ', '') == norm_name:
                            has_smiles = True
                            break
            
            if smiles_db and not has_smiles:
                continue  # Skip drugs without SMILES
                
            sim = calculate_profile_similarity(target_profile, drug_data["receptors"])
            similarities.append((drug_name, sim, drug_data))
    
    # Sort by similarity (descending)
    similarities.sort(key=lambda x: -x[1])
    
    return similarities[:top_k]


def calculate_drug_likeness(mol) -> Tuple[float, int]:
    """Calculate QED score and Lipinski violations."""
    if not RDKIT_AVAILABLE or mol is None:
        return 0.0, 5
    
    try:
        # Calculate QED (Quantitative Estimate of Drug-likeness)
        qed = Descriptors.qed(mol)
        
        # Calculate Lipinski violations
        violations = 0
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        
        return qed, violations
    except:
        return 0.0, 5


def get_molecular_fingerprint(mol):
    """Get Morgan fingerprint for a molecule."""
    if not RDKIT_AVAILABLE or mol is None:
        return None
    
    try:
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return fpgen.GetFingerprint(mol)
    except:
        return None


def calculate_tanimoto_similarity(mol1, mol2) -> float:
    """Calculate Tanimoto similarity between two molecules."""
    if not RDKIT_AVAILABLE:
        return 0.0
    
    fp1 = get_molecular_fingerprint(mol1)
    fp2 = get_molecular_fingerprint(mol2)
    
    if fp1 is None or fp2 is None:
        return 0.0
    
    return DataStructs.TanimotoSimilarity(fp1, fp2)


# ============================================================================
# SOTA MUTATION ENGINE
# ============================================================================

# Matched Molecular Pair (MMP) bioisosteric replacements
# These are validated medicinal chemistry transformations
BIOISOSTERIC_REPLACEMENTS = [
    # Halogen swaps
    ("[F]", "[Cl]"),
    ("[Cl]", "[Br]"),
    ("[Br]", "[I]"),
    ("[F]", "[Br]"),
    # Hydrogen/halogen
    ("[CH3]", "[CF3]"),
    ("[CH3]", "[Cl]"),
    # Carboxylic acid bioisosteres
    ("C(=O)O", "c1nnn[nH]1"),  # COOH -> tetrazole
    ("C(=O)O", "S(=O)(=O)N"),  # COOH -> sulfonamide
    # Amide bioisosteres
    ("C(=O)N", "C(=N)N"),  # amide -> amidine
    ("C(=O)N", "c1ncnc1"),  # amide -> pyrimidine
    # Phenyl ring modifications
    ("c1ccccc1", "c1ccncc1"),  # benzene -> pyridine
    ("c1ccccc1", "c1ccc2[nH]ccc2c1"),  # benzene -> indole
    ("c1ccccc1", "c1ccoc1"),  # benzene -> furan
    ("c1ccccc1", "c1ccsc1"),  # benzene -> thiophene
    # Ether/thioether
    ("CO", "CS"),
    ("COC", "CSC"),
    # Amine modifications (common in psychoactives)
    ("CN(C)C", "CN1CCCC1"),  # dimethylamine -> pyrrolidine
    ("CNCC", "CN1CCCCC1"),  # ethylmethylamine -> piperidine
    ("NCC", "N1CCCC1"),  # ethylamine -> pyrrolidine
    ("CN", "C[N+]([O-])"),  # amine -> N-oxide
    # Hydroxyl bioisosteres
    ("O", "S"),  # OH -> SH
    ("O", "N"),  # OH -> NH2 (simplified)
    # Methoxy modifications (important for phenethylamines)
    ("COc", "CSc"),  # OMe -> SMe
    ("COc", "Fc"),   # OMe -> F
    # Ring expansions/contractions
    ("C1CC1", "C1CCC1"),  # cyclopropyl -> cyclobutyl
    ("C1CCC1", "C1CCCC1"),  # cyclobutyl -> cyclopentyl
    ("C1CCCC1", "C1CCCCC1"),  # cyclopentyl -> cyclohexyl
]

# SELFIES vocabulary for random mutations
SELFIES_ATOMS = ['[C]', '[N]', '[O]', '[S]', '[F]', '[Cl]', '[Br]', '[I]', '[P]']
SELFIES_BONDS = ['[=C]', '[=N]', '[=O]', '[#C]', '[#N]']
SELFIES_RINGS = ['[Ring1]', '[Ring2]', '[Branch1]', '[Branch2]']


def smiles_to_selfies(smiles: str) -> Optional[str]:
    """Convert SMILES to SELFIES representation."""
    if not SELFIES_AVAILABLE:
        return None
    try:
        return sf.encoder(smiles)
    except:
        return None


def selfies_to_smiles(selfies_str: str) -> Optional[str]:
    """Convert SELFIES back to SMILES (always valid)."""
    if not SELFIES_AVAILABLE:
        return None
    try:
        smiles = sf.decoder(selfies_str)
        # Verify with RDKit
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
        return smiles
    except:
        return None


def mutate_selfies(selfies_str: str) -> Optional[str]:
    """Apply random mutation at SELFIES token level - always produces valid molecules."""
    if not SELFIES_AVAILABLE:
        return None
    
    try:
        # Parse SELFIES into tokens
        tokens = list(sf.split_selfies(selfies_str))
        if len(tokens) < 3:
            return selfies_str
        
        mutation_type = random.choice(['substitute', 'insert', 'delete', 'swap'])
        
        if mutation_type == 'substitute':
            # Replace a random token
            idx = random.randint(0, len(tokens) - 1)
            # Choose replacement based on token type
            if tokens[idx] in SELFIES_ATOMS:
                tokens[idx] = random.choice(SELFIES_ATOMS)
            elif tokens[idx] in SELFIES_BONDS:
                tokens[idx] = random.choice(SELFIES_BONDS + SELFIES_ATOMS)
            else:
                tokens[idx] = random.choice(SELFIES_ATOMS + ['[=C]', '[=N]'])
        
        elif mutation_type == 'insert':
            # Insert a new token
            idx = random.randint(0, len(tokens))
            new_token = random.choice(SELFIES_ATOMS + ['[=C]', '[Branch1]', '[Ring1]'])
            tokens.insert(idx, new_token)
        
        elif mutation_type == 'delete' and len(tokens) > 5:
            # Delete a random token (keep minimum size)
            idx = random.randint(0, len(tokens) - 1)
            tokens.pop(idx)
        
        elif mutation_type == 'swap' and len(tokens) > 2:
            # Swap two adjacent tokens
            idx = random.randint(0, len(tokens) - 2)
            tokens[idx], tokens[idx + 1] = tokens[idx + 1], tokens[idx]
        
        return ''.join(tokens)
    except:
        return selfies_str


def apply_bioisosteric_replacement(smiles: str) -> Optional[str]:
    """Apply a medicinal chemistry bioisosteric replacement."""
    if not RDKIT_AVAILABLE:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Shuffle replacements to try random ones
    replacements = BIOISOSTERIC_REPLACEMENTS.copy()
    random.shuffle(replacements)
    
    for pattern, replacement in replacements:
        try:
            # Try substructure replacement
            patt = Chem.MolFromSmarts(pattern)
            repl = Chem.MolFromSmiles(replacement)
            
            if patt is None or repl is None:
                continue
            
            if mol.HasSubstructMatch(patt):
                # Perform replacement
                new_mol = AllChem.ReplaceSubstructs(mol, patt, repl, replaceAll=False)[0]
                try:
                    Chem.SanitizeMol(new_mol)
                    new_smiles = Chem.MolToSmiles(new_mol, canonical=True)
                    if new_smiles != smiles:  # Ensure we actually changed something
                        return new_smiles
                except:
                    continue
        except:
            continue
    
    return None


def mutate_functional_group(smiles: str) -> Optional[str]:
    """Add or modify functional groups relevant to neuroactive compounds."""
    if not RDKIT_AVAILABLE:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Functional groups important for receptor binding
    neuroactive_groups = [
        ('C', 'methyl'),
        ('CC', 'ethyl'),
        ('C(C)C', 'isopropyl'),
        ('OC', 'methoxy'),
        ('O', 'hydroxy'),
        ('F', 'fluoro'),
        ('Cl', 'chloro'),
        ('Br', 'bromo'),
        ('I', 'iodo'),
        ('N(C)C', 'dimethylamino'),
        ('NC', 'methylamino'),
        ('N', 'amino'),
        ('C#N', 'cyano'),
        ('C(=O)N', 'amide'),
        ('S(=O)(=O)N', 'sulfonamide'),
        ('c1ccccc1', 'phenyl'),
    ]
    
    try:
        rwmol = Chem.RWMol(mol)
        
        # Find aromatic carbons or carbons with H for substitution
        candidates = []
        for atom in rwmol.GetAtoms():
            if atom.GetSymbol() == 'C':
                if atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
                    candidates.append((atom.GetIdx(), 'aromatic'))
                elif atom.GetTotalNumHs() > 0:
                    candidates.append((atom.GetIdx(), 'aliphatic'))
        
        if not candidates:
            return None
        
        target_idx, target_type = random.choice(candidates)
        group_smiles, group_name = random.choice(neuroactive_groups)
        
        # Create the functional group
        group = Chem.MolFromSmiles(group_smiles)
        if group is None:
            return None
        
        # Simple addition: add the group and connect
        combined = Chem.CombineMols(rwmol, group)
        ed_combined = Chem.RWMol(combined)
        
        # Find connection point in the group (first carbon or heteroatom)
        group_connect = rwmol.GetNumAtoms()  # First atom of added group
        
        # Add bond
        ed_combined.AddBond(target_idx, group_connect, BondType.SINGLE)
        
        try:
            Chem.SanitizeMol(ed_combined)
            return Chem.MolToSmiles(ed_combined, canonical=True)
        except:
            return None
            
    except:
        return None


def mutate_molecule(smiles: str, mutation_type: str = "random") -> Optional[str]:
    """
    Apply SOTA mutation to a molecule.
    
    Uses a combination of:
    1. SELFIES mutations (guaranteed valid)
    2. Bioisosteric replacements (medicinal chemistry knowledge)
    3. Functional group modifications (neuroactive-relevant)
    """
    if not RDKIT_AVAILABLE:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    if mutation_type == "random":
        # Weight towards more sophisticated mutations
        mutation_type = random.choices(
            ['selfies', 'bioisostere', 'functional', 'legacy'],
            weights=[0.35, 0.30, 0.25, 0.10]
        )[0]
    
    result = None
    
    if mutation_type == 'selfies' and SELFIES_AVAILABLE:
        # SELFIES-based mutation (always valid)
        selfies_str = smiles_to_selfies(smiles)
        if selfies_str:
            mutated_selfies = mutate_selfies(selfies_str)
            if mutated_selfies:
                result = selfies_to_smiles(mutated_selfies)
    
    elif mutation_type == 'bioisostere':
        # Bioisosteric replacement
        result = apply_bioisosteric_replacement(smiles)
    
    elif mutation_type == 'functional':
        # Functional group modification
        result = mutate_functional_group(smiles)
    
    elif mutation_type == 'legacy':
        # Fall back to simple mutations
        result = _legacy_mutate(smiles)
    
    # Validate result
    if result:
        mol = Chem.MolFromSmiles(result)
        if mol and result != smiles:
            return Chem.MolToSmiles(mol, canonical=True)
    
    # If primary mutation failed, try another method
    for fallback in ['selfies', 'functional', 'legacy']:
        if fallback == mutation_type:
            continue
        
        if fallback == 'selfies' and SELFIES_AVAILABLE:
            selfies_str = smiles_to_selfies(smiles)
            if selfies_str:
                mutated_selfies = mutate_selfies(selfies_str)
                if mutated_selfies:
                    result = selfies_to_smiles(mutated_selfies)
                    if result and result != smiles:
                        return result
        
        elif fallback == 'functional':
            result = mutate_functional_group(smiles)
            if result and result != smiles:
                return result
        
        elif fallback == 'legacy':
            result = _legacy_mutate(smiles)
            if result and result != smiles:
                return result
    
    return smiles


def _legacy_mutate(smiles: str) -> Optional[str]:
    """Legacy simple mutation (fallback)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    try:
        rwmol = Chem.RWMol(mol)
        
        mutation = random.choice([
            "add_methyl", "add_hydroxy", "add_fluoro", "add_methoxy"
        ])
        
        carbons = [a for a in rwmol.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
        if not carbons:
            return smiles
        
        target = random.choice(carbons)
        target_idx = target.GetIdx()
        
        if mutation == "add_methyl":
            new_idx = rwmol.AddAtom(Chem.Atom(6))
            rwmol.AddBond(target_idx, new_idx, BondType.SINGLE)
        elif mutation == "add_hydroxy":
            new_idx = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(target_idx, new_idx, BondType.SINGLE)
        elif mutation == "add_fluoro":
            new_idx = rwmol.AddAtom(Chem.Atom(9))
            rwmol.AddBond(target_idx, new_idx, BondType.SINGLE)
        elif mutation == "add_methoxy":
            o_idx = rwmol.AddAtom(Chem.Atom(8))
            c_idx = rwmol.AddAtom(Chem.Atom(6))
            rwmol.AddBond(target_idx, o_idx, BondType.SINGLE)
            rwmol.AddBond(o_idx, c_idx, BondType.SINGLE)
        
        Chem.SanitizeMol(rwmol)
        return Chem.MolToSmiles(rwmol)
    except:
        return smiles


def combine_fragments(smiles1: str, smiles2: str) -> Optional[str]:
    """
    Combine two molecular fragments using SOTA approaches.
    
    Strategies:
    1. BRICS decomposition and recombination
    2. Ring-preserving fragment merging
    3. Linker-based combination
    """
    if not RDKIT_AVAILABLE:
        return None
    
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return None
    
    strategy = random.choice(['brics', 'scaffold', 'linker', 'simple'])
    
    # Strategy 1: BRICS decomposition and recombination
    if strategy == 'brics':
        try:
            from rdkit.Chem import BRICS
            
            # Get BRICS fragments from both molecules
            frags1 = list(BRICS.BRICSDecompose(mol1, minFragmentSize=3))
            frags2 = list(BRICS.BRICSDecompose(mol2, minFragmentSize=3))
            
            if frags1 and frags2:
                # Try to build new molecules from fragments
                frag1 = random.choice(frags1)
                frag2 = random.choice(frags2)
                
                # Clean up fragment notation
                frag1_clean = frag1.replace('[*]', '').replace('()', '')
                frag2_clean = frag2.replace('[*]', '').replace('()', '')
                
                # Try direct combination
                for combined_smiles in [f"{frag1_clean}{frag2_clean}", f"{frag1_clean}C{frag2_clean}"]:
                    mol = Chem.MolFromSmiles(combined_smiles)
                    if mol:
                        try:
                            Chem.SanitizeMol(mol)
                            return Chem.MolToSmiles(mol, canonical=True)
                        except:
                            continue
        except:
            pass
    
    # Strategy 2: Scaffold-preserving combination
    if strategy == 'scaffold':
        try:
            from rdkit.Chem.Scaffolds import MurckoScaffold
            
            # Get core scaffolds
            core1 = MurckoScaffold.GetScaffoldForMol(mol1)
            core2 = MurckoScaffold.GetScaffoldForMol(mol2)
            
            # Use smaller scaffold and decorate with larger molecule's substituents
            if core1.GetNumAtoms() < core2.GetNumAtoms():
                core, donor = core1, mol2
            else:
                core, donor = core2, mol1
            
            # Try to connect core with a fragment from donor
            combined = Chem.CombineMols(core, donor)
            rwmol = Chem.RWMol(combined)
            
            # Find connection points
            core_carbons = [a for a in core.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
            donor_offset = core.GetNumAtoms()
            donor_carbons = [a for a in donor.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
            
            if core_carbons and donor_carbons:
                idx1 = random.choice(core_carbons).GetIdx()
                idx2 = random.choice(donor_carbons).GetIdx() + donor_offset
                rwmol.AddBond(idx1, idx2, BondType.SINGLE)
                
                try:
                    Chem.SanitizeMol(rwmol)
                    return Chem.MolToSmiles(rwmol, canonical=True)
                except:
                    pass
        except:
            pass
    
    # Strategy 3: Linker-based combination
    if strategy == 'linker':
        try:
            # Common linkers in drug design
            linkers = ['C', 'CC', 'CCC', 'C(=O)', 'C(=O)N', 'NC(=O)', 'O', 'S', 'NC', 'CN']
            linker = random.choice(linkers)
            
            # Find attachment points
            atoms1 = [a for a in mol1.GetAtoms() if a.GetSymbol() in ['C', 'N'] and a.GetTotalNumHs() > 0]
            atoms2 = [a for a in mol2.GetAtoms() if a.GetSymbol() in ['C', 'N'] and a.GetTotalNumHs() > 0]
            
            if atoms1 and atoms2:
                # Get attachment point indices
                attach1 = random.choice(atoms1).GetIdx()
                attach2 = random.choice(atoms2).GetIdx()
                
                # Create linker molecule
                linker_mol = Chem.MolFromSmiles(linker)
                if linker_mol:
                    # Combine mol1 + linker + mol2
                    combined = Chem.CombineMols(mol1, linker_mol)
                    combined = Chem.CombineMols(combined, mol2)
                    rwmol = Chem.RWMol(combined)
                    
                    # Connect mol1 to linker
                    linker_start = mol1.GetNumAtoms()
                    rwmol.AddBond(attach1, linker_start, BondType.SINGLE)
                    
                    # Connect linker to mol2
                    linker_end = linker_start + linker_mol.GetNumAtoms() - 1
                    mol2_start = mol1.GetNumAtoms() + linker_mol.GetNumAtoms()
                    rwmol.AddBond(linker_end, mol2_start + attach2, BondType.SINGLE)
                    
                    try:
                        Chem.SanitizeMol(rwmol)
                        return Chem.MolToSmiles(rwmol, canonical=True)
                    except:
                        pass
        except:
            pass
    
    # Fallback: Simple combination
    try:
        combined = Chem.CombineMols(mol1, mol2)
        rwmol = Chem.RWMol(combined)
        
        atoms1 = [a for a in mol1.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
        atoms2_offset = mol1.GetNumAtoms()
        atoms2 = [a for a in mol2.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
        
        if atoms1 and atoms2:
            idx1 = random.choice(atoms1).GetIdx()
            idx2 = random.choice(atoms2).GetIdx() + atoms2_offset
            rwmol.AddBond(idx1, idx2, BondType.SINGLE)
            
            try:
                Chem.SanitizeMol(rwmol)
                return Chem.MolToSmiles(rwmol, canonical=True)
            except:
                pass
    except:
        pass
    
    return None


def predict_binding_profile(smiles: str, 
                           drug_db: Dict,
                           smiles_db: Dict[str, str]) -> Dict[str, float]:
    """
    Predict binding profile for a molecule based on structural similarity
    to known drugs.
    """
    if not RDKIT_AVAILABLE:
        return {}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    
    # Calculate similarity to all drugs with known profiles
    weighted_profile = defaultdict(float)
    total_weight = 0
    
    for drug_name, drug_data in drug_db.get("drugs", {}).items():
        if "receptors" not in drug_data or not drug_data["receptors"]:
            continue
        
        # Get SMILES for comparison
        drug_smiles = drug_data.get("smiles", "") or smiles_db.get(drug_name, "")
        if not drug_smiles:
            continue
        
        drug_mol = Chem.MolFromSmiles(drug_smiles)
        if drug_mol is None:
            continue
        
        # Calculate similarity
        similarity = calculate_tanimoto_similarity(mol, drug_mol)
        
        if similarity > 0.3:  # Only consider reasonably similar molecules
            weight = similarity ** 2  # Square for stronger emphasis on similar molecules
            
            for receptor, affinity in drug_data["receptors"].items():
                weighted_profile[receptor] += affinity * weight
            
            total_weight += weight
    
    # Normalize
    if total_weight > 0:
        for receptor in weighted_profile:
            weighted_profile[receptor] /= total_weight
    
    return dict(weighted_profile)


def predict_binding_profile_cached(smiles: str, 
                                   drug_db: Dict,
                                   smiles_db: Dict[str, str]) -> Dict[str, float]:
    """
    Predict binding profile with caching.
    """
    cache = get_cache()
    
    # Check cache first
    cached = cache.get_prediction(smiles)
    if cached is not None:
        return cached
    
    # Compute and cache
    profile = predict_binding_profile(smiles, drug_db, smiles_db)
    if profile:
        cache.put_prediction(smiles, profile)
    
    return profile


def run_full_docking(smiles: str, receptors: List[str] = None,
                     gpu_id: int = 0, use_cache: bool = True) -> Dict[str, float]:
    """
    Run full docking pipeline (Chai-1 + OpenMM) for a molecule.
    
    Args:
        smiles: SMILES string
        receptors: List of receptors to dock to (None = all available)
        gpu_id: GPU to use
        use_cache: Whether to use cached results
    
    Returns:
        Dict of receptor -> binding energy (kJ/mol)
    """
    try:
        from docking_pipeline import DockingPipeline
    except ImportError:
        print("    Warning: docking_pipeline not available, using prediction")
        return {}
    
    cache = get_cache()
    profile = {}
    
    with DockingPipeline(gpu_id=gpu_id) as pipeline:
        if receptors is None:
            receptors = pipeline.get_available_receptors()
        
        for receptor in receptors:
            # Check cache
            if use_cache:
                cached_energy = cache.get_docking(smiles, receptor)
                if cached_energy is not None:
                    profile[receptor] = cached_energy
                    continue
            
            # Run docking
            result = pipeline.dock_molecule(smiles, receptor, use_cache=use_cache)
            
            if result.success or result.binding_energy_kJ != 0:
                profile[receptor] = result.binding_energy_kJ
                cache.put_docking(smiles, receptor, result.binding_energy_kJ)
    
    return profile


def dock_molecules(molecules: List[Molecule], 
                   receptors: List[str] = None,
                   gpu_id: int = 0,
                   use_cache: bool = True) -> List[Molecule]:
    """
    Run full docking on a list of generated molecules.
    
    Updates the molecules in-place with docked_profile.
    """
    print(f"\nRunning full docking pipeline on {len(molecules)} molecules...")
    
    for i, mol in enumerate(molecules):
        print(f"\n  [{i+1}/{len(molecules)}] Docking {mol.name}: {mol.smiles[:50]}...")
        
        docked_profile = run_full_docking(
            mol.smiles, 
            receptors=receptors,
            gpu_id=gpu_id,
            use_cache=use_cache
        )
        
        if docked_profile:
            mol.docked_profile = docked_profile
            mol.docking_completed = True
            print(f"    ✓ Docked to {len(docked_profile)} receptors")
        else:
            print(f"    ✗ Docking failed, using predicted profile")
    
    return molecules


def generate_molecules(target_profile: Dict[str, float],
                       drug_db: Dict,
                       smiles_db: Dict[str, str],
                       n_molecules: int = 10,
                       n_iterations: int = 100) -> List[Molecule]:
    """
    Generate molecules targeting a specific binding profile.
    
    Uses an evolutionary approach:
    1. Start with similar known drugs
    2. Apply mutations
    3. Combine fragments
    4. Score and select best candidates
    """
    if not RDKIT_AVAILABLE:
        print("RDKit required for molecule generation")
        return []
    
    # Find seed molecules (similar to target profile)
    similar_drugs = find_similar_drugs(target_profile, drug_db, smiles_db, top_k=50)
    
    print(f"\nFound {len(similar_drugs)} similar drugs as seeds")
    
    # Initialize population with seed molecules
    population = []
    for drug_name, similarity, drug_data in similar_drugs[:30]:  # Use more seeds
        smiles = drug_data.get("smiles", "") or smiles_db.get(drug_name, "")
        if smiles and Chem.MolFromSmiles(smiles):
            population.append({
                "smiles": smiles,
                "source": drug_name,
                "score": similarity
            })
    
    if not population:
        print("No valid seed molecules found")
        return []
    
    print(f"Starting with {len(population)} seed molecules")
    
    # Track all unique molecules found
    all_molecules = {}  # smiles -> mol_data
    
    # Evolutionary loop
    max_pop_size = max(200, n_molecules * 2)  # Dynamic population size
    mutations_per_mol = max(5, n_molecules // 20)  # More mutations for larger requests
    
    for iteration in range(n_iterations):
        new_candidates = []
        
        # Mutate existing molecules (more mutations from top molecules)
        for mol_data in population[:20]:
            for _ in range(mutations_per_mol):
                mutated = mutate_molecule(mol_data["smiles"])
                if mutated and mutated != mol_data["smiles"]:
                    new_candidates.append({
                        "smiles": mutated,
                        "source": f"mutation of {mol_data['source']}",
                        "score": 0
                    })
        
        # Combine pairs of molecules (more combinations)
        if len(population) >= 2:
            for _ in range(mutations_per_mol * 2):
                mol1, mol2 = random.sample(population[:15], 2)
                combined = combine_fragments(mol1["smiles"], mol2["smiles"])
                if combined:
                    new_candidates.append({
                        "smiles": combined,
                        "source": f"hybrid of {mol1['source']} + {mol2['source']}",
                        "score": 0
                    })
        
        # Score new candidates
        for candidate in new_candidates:
            mol = Chem.MolFromSmiles(candidate["smiles"])
            if mol is None:
                continue
            
            # Predict binding profile (with caching)
            predicted = predict_binding_profile_cached(candidate["smiles"], drug_db, smiles_db)
            
            if predicted:
                # Calculate score (similarity to target + drug-likeness)
                profile_sim = calculate_profile_similarity(target_profile, predicted)
                qed, violations = calculate_drug_likeness(mol)
                
                # Penalize Lipinski violations
                drug_likeness_penalty = violations * 0.1
                
                candidate["score"] = profile_sim * (1 - drug_likeness_penalty) * (0.5 + 0.5 * qed)
                candidate["predicted_profile"] = predicted
                candidate["qed"] = qed
                candidate["violations"] = violations
                
                population.append(candidate)
                
                # Track all valid molecules with minimum score
                canonical = Chem.CanonSmiles(candidate["smiles"])
                if canonical and candidate["score"] > 0.5:
                    if canonical not in all_molecules or all_molecules[canonical]["score"] < candidate["score"]:
                        all_molecules[canonical] = candidate
        
        # Sort by score and keep top molecules
        population = [p for p in population if p.get("score", 0) > 0]
        population.sort(key=lambda x: -x.get("score", 0))
        population = population[:max_pop_size]
        
        if (iteration + 1) % 20 == 0:
            print(f"Iteration {iteration + 1}: Best score = {population[0]['score']:.3f}, unique molecules = {len(all_molecules)}")
    
    # Convert to Molecule objects and return top results
    # Use all_molecules instead of just population for more diversity
    results = []
    seen_smiles = set()
    
    # Sort all molecules by score
    sorted_molecules = sorted(all_molecules.items(), key=lambda x: -x[1].get("score", 0))
    
    for canonical, mol_data in sorted_molecules:
        if canonical in seen_smiles:
            continue
        seen_smiles.add(canonical)
        
        mol = Molecule(
            smiles=canonical,
            name=f"Generated_{len(results)+1}",
            predicted_profile=mol_data.get("predicted_profile", {}),
            source_molecules=[mol_data.get("source", "unknown")],
            score=mol_data.get("score", 0),
            qed=mol_data.get("qed", 0),
            lipinski_violations=mol_data.get("violations", 0)
        )
        results.append(mol)
        
        if len(results) >= n_molecules:
            break
    
    # If we don't have enough from all_molecules, also check population
    if len(results) < n_molecules:
        for mol_data in population:
            smiles = mol_data["smiles"]
            canonical = Chem.CanonSmiles(smiles) if RDKIT_AVAILABLE else smiles
            if canonical in seen_smiles:
                continue
            seen_smiles.add(canonical)
            
            mol = Molecule(
                smiles=canonical,
                name=f"Generated_{len(results)+1}",
                predicted_profile=mol_data.get("predicted_profile", {}),
                source_molecules=[mol_data.get("source", "unknown")],
                score=mol_data.get("score", 0),
                qed=mol_data.get("qed", 0),
                lipinski_violations=mol_data.get("violations", 0)
            )
            results.append(mol)
            
            if len(results) >= n_molecules:
                break
    
    return results


def parse_profile_string(profile_str: str) -> Dict[str, float]:
    """Parse a profile string like '5-HT2A:-90000,D2:-50000' into a dict."""
    profile = {}
    
    for part in profile_str.split(","):
        if ":" in part:
            receptor, value = part.split(":", 1)
            try:
                profile[receptor.strip()] = float(value.strip())
            except ValueError:
                pass
    
    return profile


def save_molecule_image(mol, filename: str):
    """Save a 2D structure image of a molecule."""
    if not RDKIT_AVAILABLE:
        return
    
    try:
        img = Draw.MolToImage(mol, size=(400, 400))
        img.save(filename)
    except Exception as e:
        print(f"Could not save image: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate molecules for a target binding profile"
    )
    parser.add_argument(
        "--profile", "-p",
        help="Target binding profile (e.g., '5-HT2A:-90000,D2:-50000')"
    )
    parser.add_argument(
        "--like", "-l",
        help="Generate molecules similar to this drug"
    )
    parser.add_argument(
        "--modify", "-m",
        help="Modification instructions (e.g., 'increase D2, decrease 5-HT2B')"
    )
    parser.add_argument(
        "--n-molecules", "-n", type=int, default=10,
        help="Number of molecules to generate"
    )
    parser.add_argument(
        "--iterations", "-i", type=int, default=100,
        help="Number of optimization iterations"
    )
    parser.add_argument(
        "--output", "-o", default="generated_molecules",
        help="Output directory for results"
    )
    parser.add_argument(
        "--interactive", action="store_true",
        help="Interactive mode"
    )
    # Docking options
    parser.add_argument(
        "--dock", "-d", action="store_true",
        help="Run full docking pipeline (Chai-1 + OpenMM) on generated molecules"
    )
    parser.add_argument(
        "--dock-receptors", 
        help="Comma-separated list of receptors to dock to (default: all available)"
    )
    parser.add_argument(
        "--gpu", "-g", type=int, default=0,
        help="GPU ID for docking (default: 0)"
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Disable caching"
    )
    parser.add_argument(
        "--cache-stats", action="store_true",
        help="Show cache statistics and exit"
    )
    
    args = parser.parse_args()
    
    # Show cache stats if requested
    if args.cache_stats:
        cache = get_cache()
        stats = cache.get_stats()
        print("Cache Statistics:")
        print(f"  Predictions cached: {stats['predictions_cached']}")
        print(f"  Docking entries: {stats['docking_entries']}")
        print(f"  Total docking pairs: {stats['total_docking_pairs']}")
        return
    
    # Load data
    print("Loading drug database...")
    drug_db = load_drug_database()
    smiles_db = load_smiles_database()
    
    print(f"Loaded {len(drug_db.get('drugs', {}))} drugs")
    print(f"Loaded {len(smiles_db)} SMILES strings")
    
    # Determine target profile
    target_profile = {}
    
    if args.like:
        # Get profile from reference drug
        drug_name = args.like
        drugs = drug_db.get("drugs", {})
        
        # Try to find the drug
        ref_drug = None
        for name, data in drugs.items():
            if name.lower() == drug_name.lower():
                ref_drug = data
                break
        
        if ref_drug and "receptors" in ref_drug:
            target_profile = ref_drug["receptors"].copy()
            print(f"\nUsing profile from {drug_name}:")
            for r, v in sorted(target_profile.items(), key=lambda x: x[1])[:10]:
                print(f"  {r}: {v:.0f}")
        else:
            print(f"Could not find drug: {drug_name}")
            return
        
        # Apply modifications if specified
        if args.modify:
            print(f"\nApplying modifications: {args.modify}")
            # Parse modification string
            for mod in args.modify.split(","):
                mod = mod.strip().lower()
                if "increase" in mod:
                    receptor = mod.replace("increase", "").strip()
                    if receptor in target_profile:
                        target_profile[receptor] *= 1.5  # More negative = stronger
                    else:
                        # Add new receptor target
                        target_profile[receptor] = -80000
                elif "decrease" in mod:
                    receptor = mod.replace("decrease", "").strip()
                    if receptor in target_profile:
                        target_profile[receptor] *= 0.5
    
    elif args.profile:
        target_profile = parse_profile_string(args.profile)
        print(f"\nTarget profile:")
        for r, v in sorted(target_profile.items(), key=lambda x: x[1]):
            print(f"  {r}: {v:.0f}")
    
    elif args.interactive:
        print("\n=== Interactive Molecule Generator ===")
        print("\nAvailable receptors:")
        
        # Get all receptors from database
        all_receptors = set()
        for drug_data in drug_db.get("drugs", {}).values():
            if "receptors" in drug_data:
                all_receptors.update(drug_data["receptors"].keys())
        
        print(", ".join(sorted(all_receptors)[:20]))
        
        print("\nEnter target binding values (more negative = stronger binding)")
        print("Example: 5-HT2A:-90000")
        print("Type 'done' when finished, 'quit' to exit\n")
        
        while True:
            line = input("> ").strip()
            if line.lower() == "quit":
                return
            if line.lower() == "done":
                break
            
            if ":" in line:
                parts = line.split(":")
                receptor = parts[0].strip()
                try:
                    value = float(parts[1].strip())
                    target_profile[receptor] = value
                    print(f"  Added: {receptor} = {value}")
                except ValueError:
                    print("  Invalid value, try again")
    
    else:
        # Default: generate molecules for a classic psychedelic profile
        print("\nNo profile specified. Using default psychedelic profile...")
        target_profile = {
            "5-HT2A": -85000,
            "5-HT2B": -70000,
            "5-HT2C": -65000,
            "D2": -45000,
            "5-HT1A": -55000
        }
    
    if not target_profile:
        print("No target profile specified")
        return
    
    # Generate molecules
    print(f"\nGenerating {args.n_molecules} molecules...")
    print(f"Running {args.iterations} optimization iterations...")
    
    molecules = generate_molecules(
        target_profile,
        drug_db,
        smiles_db,
        n_molecules=args.n_molecules,
        n_iterations=args.iterations
    )
    
    if not molecules:
        print("\nNo molecules generated")
        return
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    # Save results
    print(f"\n{'='*60}")
    print(f"Generated {len(molecules)} molecules")
    print(f"{'='*60}\n")
    
    results = []
    for i, mol in enumerate(molecules):
        print(f"\n{i+1}. {mol.name}")
        print(f"   SMILES: {mol.smiles}")
        print(f"   Score: {mol.score:.3f}")
        print(f"   QED: {mol.qed:.3f}")
        print(f"   Lipinski violations: {mol.lipinski_violations}")
        print(f"   Source: {mol.source_molecules[0] if mol.source_molecules else 'unknown'}")
        
        if mol.predicted_profile:
            print("   Predicted binding profile (top 5):")
            sorted_receptors = sorted(mol.predicted_profile.items(), key=lambda x: x[1])
            for receptor, value in sorted_receptors[:5]:
                print(f"      {receptor}: {value:.0f}")
        
        # Save image
        if RDKIT_AVAILABLE:
            rdkit_mol = Chem.MolFromSmiles(mol.smiles)
            if rdkit_mol:
                img_path = output_dir / f"molecule_{i+1}.png"
                save_molecule_image(rdkit_mol, str(img_path))
        
        results.append({
            "name": mol.name,
            "smiles": mol.smiles,
            "score": mol.score,
            "qed": mol.qed,
            "lipinski_violations": mol.lipinski_violations,
            "source": mol.source_molecules,
            "predicted_profile": mol.predicted_profile
        })
    
    # Save JSON results (initial, before docking)
    results_file = output_dir / "generated_molecules.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n\nResults saved to: {output_dir}/")
    print(f"  - generated_molecules.json")
    print(f"  - molecule_*.png (2D structures)")
    
    # Run full docking pipeline if requested
    if args.dock:
        print(f"\n{'='*60}")
        print("Running Full Docking Pipeline (Chai-1 + OpenMM)")
        print(f"{'='*60}\n")
        
        # Parse receptor list if provided
        dock_receptors = None
        if args.dock_receptors:
            dock_receptors = [r.strip() for r in args.dock_receptors.split(",")]
            print(f"Docking to specified receptors: {dock_receptors}")
        else:
            # Dock to receptors in target profile
            dock_receptors = list(target_profile.keys())
            print(f"Docking to profile receptors: {dock_receptors}")
        
        # Run docking
        use_cache = not args.no_cache
        docked_molecules = dock_molecules(
            molecules,
            dock_receptors,
            gpu_id=args.gpu,
            use_cache=use_cache
        )
        
        if docked_molecules:
            # Update results with docking data from molecules
            smiles_to_docked = {m.smiles: m.docked_profile for m in docked_molecules if m.docked_profile}
            
            for result in results:
                smiles = result["smiles"]
                if smiles in smiles_to_docked:
                    result["docked_profile"] = smiles_to_docked[smiles]
            
            # Save updated results
            docked_file = output_dir / "docked_molecules.json"
            with open(docked_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            # Print docking summary
            print(f"\n{'='*60}")
            print("Docking Results Summary")
            print(f"{'='*60}\n")
            
            for result in results:
                if "docked_profile" in result and result["docked_profile"]:
                    print(f"\n{result['name']}")
                    print(f"  SMILES: {result['smiles'][:50]}...")
                    sorted_receptors = sorted(
                        result["docked_profile"].items(), 
                        key=lambda x: x[1] if x[1] is not None else 0
                    )
                    for receptor, energy in sorted_receptors[:5]:
                        if energy is not None:
                            print(f"    {receptor}: {energy:.1f} kJ/mol")
            
            print(f"\n\nDocked results saved to: {docked_file}")
        else:
            print("\nNo docking results obtained")


if __name__ == "__main__":
    main()
