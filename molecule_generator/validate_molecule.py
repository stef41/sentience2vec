#!/usr/bin/env python3
"""
Validate generated molecules and score their drug-likeness.

Checks:
- Chemical validity (RDKit sanitization)
- Lipinski's Rule of Five
- Veber's rules (rotatable bonds, TPSA)
- QED (Quantitative Estimate of Drug-likeness)
- Synthetic accessibility
- PAINS filters (Pan-Assay Interference Compounds)
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem import FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available. Install with: pip install rdkit")


@dataclass
class ValidationResult:
    """Results of molecule validation."""
    smiles: str
    canonical_smiles: str
    is_valid: bool
    
    # Properties
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0  # hydrogen bond donors
    hba: int = 0  # hydrogen bond acceptors
    tpsa: float = 0.0  # topological polar surface area
    rotatable_bonds: int = 0
    ring_count: int = 0
    aromatic_rings: int = 0
    heavy_atoms: int = 0
    
    # Scores
    qed: float = 0.0
    sa_score: float = 0.0  # synthetic accessibility (1=easy, 10=hard)
    
    # Rule violations
    lipinski_violations: int = 0
    veber_violations: int = 0
    pains_alerts: List[str] = None
    
    # Overall
    drug_likeness_score: float = 0.0
    warnings: List[str] = None
    
    def __post_init__(self):
        if self.pains_alerts is None:
            self.pains_alerts = []
        if self.warnings is None:
            self.warnings = []


def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """Check if SMILES is valid and return canonical form."""
    if not RDKIT_AVAILABLE:
        return True, smiles
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, None
        
        # Try to sanitize
        Chem.SanitizeMol(mol)
        
        # Get canonical SMILES
        canonical = Chem.MolToSmiles(mol, canonical=True)
        return True, canonical
    except Exception as e:
        return False, None


def calculate_properties(mol) -> Dict:
    """Calculate molecular properties."""
    props = {}
    
    try:
        props["molecular_weight"] = Descriptors.MolWt(mol)
        props["logp"] = Descriptors.MolLogP(mol)
        props["hbd"] = rdMolDescriptors.CalcNumHBD(mol)
        props["hba"] = rdMolDescriptors.CalcNumHBA(mol)
        props["tpsa"] = Descriptors.TPSA(mol)
        props["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        props["ring_count"] = rdMolDescriptors.CalcNumRings(mol)
        props["aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        props["heavy_atoms"] = mol.GetNumHeavyAtoms()
        props["qed"] = Descriptors.qed(mol)
    except Exception as e:
        print(f"Warning: Could not calculate some properties: {e}")
    
    return props


def check_lipinski(props: Dict) -> Tuple[int, List[str]]:
    """
    Check Lipinski's Rule of Five.
    
    Rules:
    - MW ≤ 500 Da
    - LogP ≤ 5
    - HBD ≤ 5
    - HBA ≤ 10
    """
    violations = []
    
    if props.get("molecular_weight", 0) > 500:
        violations.append(f"MW > 500 ({props['molecular_weight']:.1f})")
    
    if props.get("logp", 0) > 5:
        violations.append(f"LogP > 5 ({props['logp']:.2f})")
    
    if props.get("hbd", 0) > 5:
        violations.append(f"HBD > 5 ({props['hbd']})")
    
    if props.get("hba", 0) > 10:
        violations.append(f"HBA > 10 ({props['hba']})")
    
    return len(violations), violations


def check_veber(props: Dict) -> Tuple[int, List[str]]:
    """
    Check Veber's rules for oral bioavailability.
    
    Rules:
    - Rotatable bonds ≤ 10
    - TPSA ≤ 140 Å²
    """
    violations = []
    
    if props.get("rotatable_bonds", 0) > 10:
        violations.append(f"Rotatable bonds > 10 ({props['rotatable_bonds']})")
    
    if props.get("tpsa", 0) > 140:
        violations.append(f"TPSA > 140 ({props['tpsa']:.1f})")
    
    return len(violations), violations


def check_pains(mol) -> List[str]:
    """Check for PAINS (Pan-Assay Interference Compounds) alerts."""
    if not RDKIT_AVAILABLE:
        return []
    
    alerts = []
    
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        
        entries = catalog.GetMatches(mol)
        for entry in entries:
            alerts.append(entry.GetDescription())
    except Exception as e:
        pass  # PAINS filter not critical
    
    return alerts


def calculate_sa_score(mol) -> float:
    """
    Calculate synthetic accessibility score.
    
    Returns a score from 1 (easy to synthesize) to 10 (hard).
    This is a simplified version - full SA score requires additional data.
    """
    if not RDKIT_AVAILABLE:
        return 5.0
    
    try:
        # Simplified SA score based on:
        # - Ring complexity
        # - Stereocenters
        # - Unusual elements
        # - Fragment complexity
        
        score = 1.0  # Start at easy
        
        # Penalize many rings
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        if num_rings > 4:
            score += (num_rings - 4) * 0.5
        
        # Penalize fused rings
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            # Check for complex ring systems
            for ring in ring_info.AtomRings():
                if len(ring) > 6:
                    score += 0.5
        
        # Penalize stereocenters
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        score += len(chiral_centers) * 0.3
        
        # Penalize unusual atoms
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'H']:
                score += 0.5
        
        # Cap at 10
        return min(10.0, max(1.0, score))
    
    except Exception as e:
        return 5.0


def calculate_drug_likeness_score(result: ValidationResult) -> float:
    """
    Calculate overall drug-likeness score (0-1).
    
    Combines:
    - QED score
    - Lipinski compliance
    - Veber compliance
    - SA score
    - PAINS alerts
    """
    if not result.is_valid:
        return 0.0
    
    score = result.qed  # Start with QED (0-1)
    
    # Penalize Lipinski violations
    score *= (1 - 0.1 * result.lipinski_violations)
    
    # Penalize Veber violations
    score *= (1 - 0.1 * result.veber_violations)
    
    # Penalize high SA score
    if result.sa_score > 5:
        score *= (1 - (result.sa_score - 5) * 0.05)
    
    # Penalize PAINS alerts
    score *= (1 - 0.15 * len(result.pains_alerts))
    
    return max(0.0, min(1.0, score))


def validate_molecule(smiles: str) -> ValidationResult:
    """Perform full validation on a molecule."""
    result = ValidationResult(
        smiles=smiles,
        canonical_smiles="",
        is_valid=False
    )
    
    # Check validity
    is_valid, canonical = validate_smiles(smiles)
    if not is_valid:
        result.warnings.append("Invalid SMILES")
        return result
    
    result.is_valid = True
    result.canonical_smiles = canonical
    
    if not RDKIT_AVAILABLE:
        return result
    
    mol = Chem.MolFromSmiles(smiles)
    
    # Calculate properties
    props = calculate_properties(mol)
    result.molecular_weight = props.get("molecular_weight", 0)
    result.logp = props.get("logp", 0)
    result.hbd = props.get("hbd", 0)
    result.hba = props.get("hba", 0)
    result.tpsa = props.get("tpsa", 0)
    result.rotatable_bonds = props.get("rotatable_bonds", 0)
    result.ring_count = props.get("ring_count", 0)
    result.aromatic_rings = props.get("aromatic_rings", 0)
    result.heavy_atoms = props.get("heavy_atoms", 0)
    result.qed = props.get("qed", 0)
    
    # Check rules
    lipinski_count, lipinski_warnings = check_lipinski(props)
    result.lipinski_violations = lipinski_count
    result.warnings.extend(lipinski_warnings)
    
    veber_count, veber_warnings = check_veber(props)
    result.veber_violations = veber_count
    result.warnings.extend(veber_warnings)
    
    # Check PAINS
    result.pains_alerts = check_pains(mol)
    if result.pains_alerts:
        result.warnings.extend([f"PAINS alert: {a}" for a in result.pains_alerts])
    
    # Calculate SA score
    result.sa_score = calculate_sa_score(mol)
    if result.sa_score > 6:
        result.warnings.append(f"High synthetic difficulty (SA={result.sa_score:.1f})")
    
    # Calculate overall score
    result.drug_likeness_score = calculate_drug_likeness_score(result)
    
    return result


def validate_molecules_file(filepath: str) -> List[ValidationResult]:
    """Validate all molecules in a JSON file."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    results = []
    for mol_data in data:
        smiles = mol_data.get("smiles", "")
        if smiles:
            result = validate_molecule(smiles)
            results.append(result)
    
    return results


def print_validation_report(result: ValidationResult):
    """Print a detailed validation report for a molecule."""
    print(f"\n{'='*60}")
    print(f"SMILES: {result.smiles}")
    print(f"Canonical: {result.canonical_smiles}")
    print(f"Valid: {'✓' if result.is_valid else '✗'}")
    
    if not result.is_valid:
        print("Validation failed!")
        return
    
    print(f"\n{'─'*30} Properties {'─'*30}")
    print(f"  Molecular Weight: {result.molecular_weight:.1f} Da")
    print(f"  LogP: {result.logp:.2f}")
    print(f"  H-Bond Donors: {result.hbd}")
    print(f"  H-Bond Acceptors: {result.hba}")
    print(f"  TPSA: {result.tpsa:.1f} Å²")
    print(f"  Rotatable Bonds: {result.rotatable_bonds}")
    print(f"  Ring Count: {result.ring_count}")
    print(f"  Aromatic Rings: {result.aromatic_rings}")
    print(f"  Heavy Atoms: {result.heavy_atoms}")
    
    print(f"\n{'─'*30} Scores {'─'*30}")
    print(f"  QED: {result.qed:.3f}")
    print(f"  Synthetic Accessibility: {result.sa_score:.1f}/10")
    print(f"  Drug-likeness Score: {result.drug_likeness_score:.3f}")
    
    print(f"\n{'─'*30} Rule Compliance {'─'*30}")
    print(f"  Lipinski violations: {result.lipinski_violations}")
    print(f"  Veber violations: {result.veber_violations}")
    print(f"  PAINS alerts: {len(result.pains_alerts)}")
    
    if result.warnings:
        print(f"\n{'─'*30} Warnings {'─'*30}")
        for warning in result.warnings:
            print(f"  ⚠ {warning}")
    
    # Overall assessment
    print(f"\n{'─'*30} Assessment {'─'*30}")
    if result.drug_likeness_score >= 0.7:
        print("  ✓ GOOD: Molecule shows good drug-like properties")
    elif result.drug_likeness_score >= 0.4:
        print("  ~ MODERATE: Some drug-likeness concerns")
    else:
        print("  ✗ POOR: Significant drug-likeness issues")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate molecules")
    parser.add_argument("input", help="SMILES string or path to JSON file")
    parser.add_argument("--json", action="store_true", help="Output as JSON")
    
    args = parser.parse_args()
    
    if args.input.endswith(".json"):
        # Validate file
        results = validate_molecules_file(args.input)
        
        if args.json:
            output = []
            for r in results:
                output.append({
                    "smiles": r.smiles,
                    "canonical_smiles": r.canonical_smiles,
                    "is_valid": r.is_valid,
                    "molecular_weight": r.molecular_weight,
                    "logp": r.logp,
                    "qed": r.qed,
                    "drug_likeness_score": r.drug_likeness_score,
                    "lipinski_violations": r.lipinski_violations,
                    "warnings": r.warnings
                })
            print(json.dumps(output, indent=2))
        else:
            for result in results:
                print_validation_report(result)
    else:
        # Validate single SMILES
        result = validate_molecule(args.input)
        
        if args.json:
            print(json.dumps({
                "smiles": result.smiles,
                "canonical_smiles": result.canonical_smiles,
                "is_valid": result.is_valid,
                "molecular_weight": result.molecular_weight,
                "logp": result.logp,
                "qed": result.qed,
                "drug_likeness_score": result.drug_likeness_score,
                "lipinski_violations": result.lipinski_violations,
                "warnings": result.warnings
            }, indent=2))
        else:
            print_validation_report(result)


if __name__ == "__main__":
    main()
