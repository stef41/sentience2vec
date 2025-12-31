#!/usr/bin/env python3
"""
FULL COMPREHENSIVE Preparation Pipeline

Prepares ALL ~320 psychoactive compounds and ALL ~77 neuroreceptor structures
for molecular dynamics simulations.
"""

import json
import os
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# RDKit for ligand 3D generation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

BASE_DIR = Path(__file__).parent
LIGANDS_DIR = BASE_DIR / "ligands_full"
RECEPTORS_DIR = BASE_DIR / "prepared_receptors_full"


def get_pubchem_smiles(name, retries=2):
    """Fetch SMILES from PubChem."""
    # Name variations to try
    variations = [name]
    
    # Add common variations
    name_lower = name.lower()
    if 'thc' in name_lower:
        variations.append('tetrahydrocannabinol')
    if name == 'LSD':
        variations.append('lysergic acid diethylamide')
    if name == 'DMT':
        variations.append('N,N-dimethyltryptamine')
    if name == 'MDMA':
        variations.append('3,4-methylenedioxymethamphetamine')
    
    for var in variations:
        for attempt in range(retries):
            try:
                encoded = urllib.parse.quote(var)
                url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES/JSON'
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=10) as response:
                    data = json.loads(response.read())
                    props = data['PropertyTable']['Properties'][0]
                    # Try different SMILES fields
                    return props.get('CanonicalSMILES') or props.get('ConnectivitySMILES') or props.get('IsomericSMILES')
            except:
                time.sleep(0.5)
    return None


def download_pubchem_sdf(name, output_path, retries=2):
    """Download 3D SDF from PubChem."""
    variations = [name]
    if name == 'LSD':
        variations.append('lysergic acid diethylamide')
    if name == 'DMT':
        variations.append('N,N-dimethyltryptamine')
    if name == 'THC':
        variations.append('tetrahydrocannabinol')
    
    for var in variations:
        for attempt in range(retries):
            try:
                # First get CID
                encoded = urllib.parse.quote(var)
                url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/JSON'
                with urllib.request.urlopen(url, timeout=10) as r:
                    data = json.loads(r.read())
                    cid = data['IdentifierList']['CID'][0]
                
                # Download 3D SDF
                url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d'
                urllib.request.urlretrieve(url, output_path)
                return True
            except:
                time.sleep(0.5)
    return False


def generate_3d_from_smiles(smiles, output_pdb, output_sdf):
    """Generate 3D structure from SMILES using RDKit."""
    if not HAS_RDKIT:
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        mol = Chem.AddHs(mol)
        
        # Try embedding
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        
        if result != 0:
            return False
        
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        
        Chem.MolToPDBFile(mol, str(output_pdb))
        
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(mol)
        writer.close()
        
        return True
    except:
        return False


def prepare_ligand(name):
    """Prepare a single ligand - try PubChem first, then RDKit from SMILES."""
    ligand_dir = LIGANDS_DIR / name.replace('/', '_').replace(' ', '_')
    ligand_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_file = ligand_dir / f"{name.replace('/', '_')}.pdb"
    sdf_file = ligand_dir / f"{name.replace('/', '_')}.sdf"
    
    # Skip if already done
    if pdb_file.exists() or sdf_file.exists():
        return 'cached'
    
    # Try downloading 3D from PubChem
    if download_pubchem_sdf(name, sdf_file):
        # Convert to PDB if we have RDKit
        if HAS_RDKIT:
            try:
                supplier = Chem.SDMolSupplier(str(sdf_file))
                mol = next(iter(supplier))
                if mol:
                    mol = Chem.AddHs(mol, addCoords=True)
                    Chem.MolToPDBFile(mol, str(pdb_file))
            except:
                pass
        return 'pubchem'
    
    # Fall back to generating from SMILES
    smiles = get_pubchem_smiles(name)
    if smiles and generate_3d_from_smiles(smiles, pdb_file, sdf_file):
        return 'rdkit'
    
    return 'failed'


def download_pdb(pdb_id, output_path):
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except:
        return False


def main():
    print("=" * 70)
    print("FULL COMPREHENSIVE PREPARATION")
    print("=" * 70)
    
    # Load all compounds from PsychonautWiki
    pw_file = BASE_DIR.parent / 'data/raw/psychonautwiki/substances.json'
    with open(pw_file) as f:
        all_substances = json.load(f)
    
    # Filter to actual compounds
    CATEGORY_KEYWORDS = [
        'substituted', 'list/', 'classical', 'deliriant', 'depressant',
        'stimulant', 'dissociative', 'opioid', 'psychedelic', 'entheogen',
        'nootropic', 'eugeroic', 'neurotransmitter', 'monoamine', 'cannabinoid',
        'barbiturate', 'benzodiazepine', 'racetam', 'adamantane', 'xanthine',
        'inhalant', 'arylcyclohexylamine', 'diarylethylamine', 'sigmaergic',
        'alcohols', 'gabapentinoid', 'thienodiazepine', 'antihist'
    ]
    
    compounds = []
    for s in all_substances:
        name = s.get('name', '')
        if not name or len(name) < 2:
            continue
        if any(kw in name.lower() for kw in CATEGORY_KEYWORDS):
            continue
        compounds.append(name)
    
    compounds = list(set(compounds))  # Remove duplicates
    print(f"\n📊 Total compounds to prepare: {len(compounds)}")
    
    # Load receptor database
    from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE
    receptors_with_pdb = {k: v for k, v in RECEPTOR_PDB_DATABASE.items() if v.get('pdb_id')}
    print(f"📊 Total receptors with PDB: {len(receptors_with_pdb)}")
    
    # Create directories
    LIGANDS_DIR.mkdir(exist_ok=True)
    RECEPTORS_DIR.mkdir(exist_ok=True)
    
    # Prepare ligands
    print("\n" + "=" * 70)
    print("LIGAND PREPARATION")
    print("=" * 70)
    
    results = {'success': 0, 'cached': 0, 'failed': 0}
    failed_ligands = []
    
    for i, name in enumerate(sorted(compounds), 1):
        print(f"  [{i}/{len(compounds)}] {name}...", end=' ', flush=True)
        status = prepare_ligand(name)
        
        if status in ['pubchem', 'rdkit']:
            print(f"✓ ({status})")
            results['success'] += 1
        elif status == 'cached':
            print("(cached)")
            results['cached'] += 1
        else:
            print("✗")
            results['failed'] += 1
            failed_ligands.append(name)
        
        time.sleep(0.1)  # Rate limiting
    
    print(f"\nLigand Results: {results['success']} new, {results['cached']} cached, {results['failed']} failed")
    
    # Save results
    with open(BASE_DIR / 'full_preparation_results.json', 'w') as f:
        json.dump({
            'ligands': {
                'total': len(compounds),
                'prepared': results['success'] + results['cached'],
                'failed': failed_ligands
            },
            'receptors': {
                'total': len(receptors_with_pdb),
                'list': list(receptors_with_pdb.keys())
            }
        }, f, indent=2)
    
    print(f"\n✓ Results saved to full_preparation_results.json")


if __name__ == '__main__':
    main()
