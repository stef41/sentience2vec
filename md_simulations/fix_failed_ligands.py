#!/usr/bin/env python3
"""Fix failed ligands with correct PubChem names or SMILES"""

import os
import json
import urllib.request
from pathlib import Path

# Map failed compounds to correct PubChem names or SMILES
FIXES = {
    # Actual compounds with alternate names
    "DET": "N,N-Diethyltryptamine",
    "DPT": "N,N-Dipropyltryptamine", 
    "EPT": "N-Ethyl-N-propyltryptamine",
    "MPT": "N-Methyl-N-propyltryptamine",
    "MDA": "3,4-Methylenedioxyamphetamine",
    "PCP": "Phencyclidine",
    "PMMA": "para-Methoxymethamphetamine",
    "NEP": "N-Ethylpentedrone",
    "A-PHP": "alpha-Pyrrolidinohexiophenone",
    "HHC": "Hexahydrocannabinol",
    "HXE": "Hydroxynorketamine",  # Actually 3-HO-2'-Oxo-PCE
    "F-Phenibut": "Fluorophenibut",
    "Alpha-GPC": "L-Alpha glycerylphosphorylcholine",
    "3-Cl-PCP": "3-Chlorophencyclidine",
    "4F-EPH": "4-Fluoroethylphenidate",
    "ΑMT": "alpha-Methyltryptamine",  # Greek alpha
    "Βk-2C-B": "beta-Keto 2C-B",
    "Δ-8-Tetrahydrocannabinol": "Delta-8-THC",
    "Δ-10-Tetrahydrocannabinol": "Delta-10-THC",
    "Δ-11-Tetrahydrocannabinol": "Delta-11-THC",
    "Tryptamine (compound)": "Tryptamine",
    "Phenethylamine (compound)": "Phenethylamine",
    "Harmala alkaloid": "Harmine",  # Main harmala alkaloid
    
    # Plant extracts - use main active compound
    "Kratom": "Mitragynine",
    "Salvia divinorum": "Salvinorin A",  # Already have this
    "Ayahuasca": "N,N-Dimethyltryptamine",  # Already have DMT
    "Changa": "N,N-Dimethyltryptamine",  # DMT mixture
    "Datura": "Scopolamine",
    "Amanita muscaria": "Muscimol",
    "Poppers": "Isoamyl nitrite",
    "Nitrous": "Nitrous oxide",
    
    # Generic categories - skip
    "25x-NBOMe": None,  # Generic category
    "2C-T-x": None,
    "2C-x": None,
    "MDxx": None,
    "Lysergamides": None,
    "LSA adducts": None,
    "Oneirogen": None,  # Category
    
    # Botanical names - skip (need active compound)
    "Anadenanthera peregrina": "Bufotenin",  # Main active
    "Argyreia nervosa": "Ergine",  # LSA
    "Banisteriopsis caapi": "Harmine",
    "Calea ternifolia": None,  # Complex mix
    "Echinopsis pachanoi": "Mescaline",
    "Entada rheedii": None,  # Unknown actives
    "Ipomoea tricolor": "Ergine",  # LSA
    "Lophophora williamsii": "Mescaline",
    "Peganum harmala": "Harmine",
    "Tabernanthe iboga (botany)": "Ibogaine",
    "Nicotiana (botany)": "Nicotine",
    
    # Novel cannabinoids
    "Tetrahydrocannabihexol": "THCH",  # Try abbreviation
    "Tetrahydrocannabutol": "THCB",
}

# Direct SMILES for compounds not in PubChem
SMILES_FIXES = {
    "N-Methylbisfluoromodafinil": "CC(C(=O)NC)C(C1=CC=C(C=C1)F)SC2=CC=C(C=C2)F",  # Approximate
    "HXE": "CCC(NC)C1=CC=CC=C1O",  # 3-HO-2'-oxo-PCE approximation
}

def download_pubchem_sdf(name, outdir):
    """Download 3D SDF from PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/SDF?record_type=3d"
    sdf_path = outdir / f"{name.replace(' ', '_')}.sdf"
    try:
        urllib.request.urlretrieve(url, sdf_path)
        if sdf_path.stat().st_size > 100:
            return sdf_path
        sdf_path.unlink()
    except:
        pass
    return None

def get_smiles(name):
    """Get SMILES from PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/property/CanonicalSMILES/JSON"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
            props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
            return props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
    except:
        return None

def sdf_to_pdb(sdf_path, pdb_path):
    """Convert SDF to PDB using RDKit"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_path))
        if mol:
            Chem.MolToPDBFile(mol, str(pdb_path))
            return True
    except:
        pass
    return False

def generate_3d_from_smiles(smiles, name, outdir):
    """Generate 3D structure from SMILES"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
            
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            if AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42) != 0:
                return False
        
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        
        sdf_path = outdir / f"{name.replace(' ', '_')}.sdf"
        pdb_path = outdir / f"{name.replace(' ', '_')}.pdb"
        
        Chem.MolToMolFile(mol, str(sdf_path))
        Chem.MolToPDBFile(mol, str(pdb_path))
        
        return True
    except Exception as e:
        print(f"    RDKit error: {e}")
        return False

def main():
    import urllib.parse
    
    base_dir = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/ligands_full")
    
    fixed = 0
    skipped = 0
    failed = 0
    
    for original_name, fix in FIXES.items():
        if fix is None:
            print(f"  ⏭ {original_name} (generic category)")
            skipped += 1
            continue
            
        outdir = base_dir / original_name.replace(" ", "_")
        outdir.mkdir(exist_ok=True)
        
        # Check if already has structure
        has_sdf = any(outdir.glob("*.sdf"))
        if has_sdf:
            print(f"  ✓ {original_name} (already done)")
            continue
        
        print(f"  → {original_name} => {fix}")
        
        # Try PubChem 3D
        sdf = download_pubchem_sdf(fix, outdir)
        if sdf:
            pdb_path = outdir / f"{fix.replace(' ', '_')}.pdb"
            if sdf_to_pdb(sdf, pdb_path):
                print(f"    ✓ Downloaded from PubChem")
                fixed += 1
                continue
        
        # Try SMILES
        smiles = SMILES_FIXES.get(original_name) or get_smiles(fix)
        if smiles and generate_3d_from_smiles(smiles, fix, outdir):
            print(f"    ✓ Generated from SMILES")
            fixed += 1
            continue
            
        print(f"    ✗ Failed")
        failed += 1
    
    print(f"\n{'='*60}")
    print(f"Fixed: {fixed}, Skipped: {skipped}, Failed: {failed}")

if __name__ == "__main__":
    main()
