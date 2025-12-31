#!/usr/bin/env python3
"""
Fetch InChI (International Chemical Identifier) for all substances
==================================================================
Maps drug names to their InChI codes using PubChem API.

InChI provides:
- Unique chemical identification
- Cross-database linking (PubChem, ChEMBL, DrugBank)
- Molecular structure verification
- Isomer differentiation
"""

import json
import time
import re
import requests
from pathlib import Path
from collections import defaultdict
from typing import Optional
from datetime import datetime

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "raw"
OUTPUT_DIR = BASE_DIR / "data" / "chemical_identifiers"

# PubChem API endpoints
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Common name variations/aliases
SUBSTANCE_ALIASES = {
    # Psychedelics
    'lsd': ['LSD', 'lysergic acid diethylamide', 'LSD-25', 'd-lysergic acid diethylamide'],
    'psilocybin': ['psilocybin', '4-phosphoryloxy-N,N-dimethyltryptamine'],
    'psilocin': ['psilocin', '4-hydroxy-N,N-dimethyltryptamine', '4-HO-DMT'],
    'dmt': ['DMT', 'N,N-Dimethyltryptamine', 'dimethyltryptamine'],
    '5-meo-dmt': ['5-MeO-DMT', '5-methoxy-N,N-dimethyltryptamine'],
    'mescaline': ['mescaline', '3,4,5-trimethoxyphenethylamine'],
    '2c-b': ['2C-B', '4-bromo-2,5-dimethoxyphenethylamine'],
    '2c-e': ['2C-E', '4-ethyl-2,5-dimethoxyphenethylamine'],
    '2c-i': ['2C-I', '4-iodo-2,5-dimethoxyphenethylamine'],
    
    # Tryptamines
    '4-aco-dmt': ['4-AcO-DMT', '4-acetoxy-N,N-dimethyltryptamine', 'psilacetin'],
    '4-ho-met': ['4-HO-MET', '4-hydroxy-N-methyl-N-ethyltryptamine'],
    'dpt': ['DPT', 'N,N-Dipropyltryptamine'],
    'amt': ['AMT', 'alpha-methyltryptamine'],
    
    # Dissociatives
    'ketamine': ['ketamine', '2-(2-chlorophenyl)-2-(methylamino)cyclohexanone'],
    'pcp': ['PCP', 'phencyclidine', '1-(1-phenylcyclohexyl)piperidine'],
    'dxm': ['DXM', 'dextromethorphan'],
    'mxe': ['MXE', 'methoxetamine', '2-(ethylamino)-2-(3-methoxyphenyl)cyclohexanone'],
    '3-meo-pcp': ['3-MeO-PCP', '3-methoxyphencyclidine'],
    'nitrous': ['nitrous oxide', 'N2O', 'dinitrogen monoxide'],
    
    # Empathogens
    'mdma': ['MDMA', '3,4-methylenedioxymethamphetamine', 'ecstasy'],
    'mda': ['MDA', '3,4-methylenedioxyamphetamine'],
    '6-apb': ['6-APB', '6-(2-aminopropyl)benzofuran'],
    
    # Cannabis
    'cannabis': ['THC', 'delta-9-tetrahydrocannabinol', 'tetrahydrocannabinol'],
    'thc': ['THC', 'delta-9-tetrahydrocannabinol', 'dronabinol'],
    'cbd': ['CBD', 'cannabidiol'],
    
    # Stimulants
    'amphetamine': ['amphetamine', 'alpha-methylphenethylamine'],
    'methamphetamine': ['methamphetamine', 'N-methyl-alpha-methylphenethylamine'],
    'cocaine': ['cocaine', 'benzoylmethylecgonine'],
    'caffeine': ['caffeine', '1,3,7-trimethylxanthine'],
    'nicotine': ['nicotine', '3-(1-methylpyrrolidin-2-yl)pyridine'],
    'modafinil': ['modafinil', '2-[(diphenylmethyl)sulfinyl]acetamide'],
    
    # Depressants
    'alcohol': ['ethanol', 'ethyl alcohol'],
    'ghb': ['GHB', 'gamma-hydroxybutyric acid', '4-hydroxybutanoic acid'],
    'gbl': ['GBL', 'gamma-butyrolactone'],
    'phenibut': ['phenibut', '4-amino-3-phenylbutanoic acid'],
    
    # Opioids
    'morphine': ['morphine'],
    'heroin': ['heroin', 'diacetylmorphine', 'diamorphine'],
    'codeine': ['codeine'],
    'oxycodone': ['oxycodone'],
    'fentanyl': ['fentanyl'],
    'kratom': ['mitragynine', '7-hydroxymitragynine'],
    
    # Benzodiazepines
    'alprazolam': ['alprazolam', 'xanax'],
    'diazepam': ['diazepam', 'valium'],
    'clonazepam': ['clonazepam', 'klonopin'],
    'lorazepam': ['lorazepam', 'ativan'],
    
    # Deliriants
    'dph': ['diphenhydramine', 'benadryl'],
    'scopolamine': ['scopolamine', 'hyoscine'],
    
    # Atypical
    'salvia': ['salvinorin A'],
    'ibogaine': ['ibogaine'],
    'muscimol': ['muscimol'],
}


def query_pubchem_by_name(name: str) -> Optional[dict]:
    """Query PubChem for compound information by name."""
    try:
        # First get CID
        url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
        resp = requests.get(url, timeout=10)
        
        if resp.status_code != 200:
            return None
        
        data = resp.json()
        cids = data.get('IdentifierList', {}).get('CID', [])
        
        if not cids:
            return None
        
        cid = cids[0]  # Take first match
        
        # Now get full properties including InChI
        props_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
        props_resp = requests.get(props_url, timeout=10)
        
        if props_resp.status_code != 200:
            return None
        
        props_data = props_resp.json()
        properties = props_data.get('PropertyTable', {}).get('Properties', [{}])[0]
        
        return {
            'cid': cid,
            'inchi': properties.get('InChI'),
            'inchi_key': properties.get('InChIKey'),
            'molecular_formula': properties.get('MolecularFormula'),
            'molecular_weight': properties.get('MolecularWeight'),
            'iupac_name': properties.get('IUPACName'),
            'canonical_smiles': properties.get('CanonicalSMILES'),
            'isomeric_smiles': properties.get('IsomericSMILES'),
        }
        
    except Exception as e:
        return None


def query_pubchem_synonyms(cid: int) -> list[str]:
    """Get all synonyms for a compound."""
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
        resp = requests.get(url, timeout=10)
        
        if resp.status_code != 200:
            return []
        
        data = resp.json()
        synonyms = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
        return synonyms[:50]  # Limit to 50 synonyms
        
    except Exception:
        return []


def normalize_name(name: str) -> str:
    """Normalize substance name for lookup."""
    name = name.lower().strip()
    # Remove common suffixes
    name = re.sub(r'\s*\(.*\)$', '', name)
    name = re.sub(r'\s*-\s*(oral|smoked|iv|insufflated)$', '', name, flags=re.I)
    return name


def collect_all_substances() -> dict[str, list[str]]:
    """Collect all substance names from our databases."""
    substances = defaultdict(set)
    
    # From PsychonautWiki
    pw_path = DATA_DIR / "psychonautwiki" / "substances.json"
    if pw_path.exists():
        with open(pw_path) as f:
            pw = json.load(f)
        for s in pw:
            name = s.get('name', '')
            if name:
                norm = normalize_name(name)
                substances[norm].add(name)
    
    # From TripSit
    ts_path = DATA_DIR / "tripsit" / "all_drugs_full.json"
    if ts_path.exists():
        with open(ts_path) as f:
            ts = json.load(f)
        for name, data in ts.items():
            if name:
                norm = normalize_name(name)
                substances[norm].add(name)
            if isinstance(data, dict):
                pretty = data.get('pretty_name', '')
                if pretty:
                    substances[norm].add(pretty)
                # Also add aliases from TripSit
                for alias in data.get('aliases', []):
                    substances[norm].add(alias)
    
    # From Ferrari
    ferrari_path = DATA_DIR / "ferrari_processed" / "summary.json"
    if ferrari_path.exists():
        with open(ferrari_path) as f:
            ferrari = json.load(f)
        for sub in ferrari.get('substance_counts', {}).keys():
            if sub:
                norm = normalize_name(sub)
                substances[norm].add(sub)
    
    # Convert sets to lists
    return {k: list(v) for k, v in substances.items() if k}


def fetch_all_inchi():
    """Fetch InChI for all substances."""
    print("=" * 70)
    print("FETCHING InChI IDENTIFIERS FROM PUBCHEM")
    print("=" * 70)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Collect substances
    print("\n📊 Collecting substance names...")
    substances = collect_all_substances()
    print(f"   Found {len(substances)} unique substances")
    
    # Load existing data if any
    output_file = OUTPUT_DIR / "substance_inchi.json"
    existing_data = {}
    if output_file.exists():
        with open(output_file) as f:
            existing_data = json.load(f)
        print(f"   Loaded {len(existing_data)} cached entries")
    
    # Fetch InChI for each substance
    results = dict(existing_data)
    success = 0
    failed = []
    skipped = 0
    
    substance_list = sorted(substances.keys())
    total = len(substance_list)
    
    print(f"\n📡 Querying PubChem for {total} substances...")
    
    for idx, norm_name in enumerate(substance_list):
        # Skip if already have
        if norm_name in results and results[norm_name].get('inchi'):
            skipped += 1
            continue
        
        if idx % 50 == 0:
            print(f"   Processing {idx}/{total}...")
        
        aliases = substances[norm_name]
        
        # Try normalized name first
        search_names = [norm_name] + list(aliases)
        
        # Also try known aliases
        if norm_name in SUBSTANCE_ALIASES:
            search_names.extend(SUBSTANCE_ALIASES[norm_name])
        
        found = False
        for search_name in search_names:
            result = query_pubchem_by_name(search_name)
            if result and result.get('inchi'):
                # Get synonyms
                synonyms = query_pubchem_synonyms(result['cid'])
                
                results[norm_name] = {
                    'normalized_name': norm_name,
                    'original_names': list(aliases),
                    'search_name_used': search_name,
                    'pubchem_cid': result['cid'],
                    'inchi': result['inchi'],
                    'inchi_key': result['inchi_key'],
                    'molecular_formula': result['molecular_formula'],
                    'molecular_weight': result['molecular_weight'],
                    'iupac_name': result['iupac_name'],
                    'canonical_smiles': result['canonical_smiles'],
                    'isomeric_smiles': result['isomeric_smiles'],
                    'synonyms': synonyms,
                }
                success += 1
                found = True
                break
            
            time.sleep(0.2)  # Rate limiting
        
        if not found:
            failed.append(norm_name)
            results[norm_name] = {
                'normalized_name': norm_name,
                'original_names': list(aliases),
                'inchi': None,
                'error': 'Not found in PubChem'
            }
        
        # Save periodically
        if idx % 100 == 0:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
    
    # Final save
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Statistics
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    with_inchi = sum(1 for r in results.values() if r.get('inchi'))
    
    print(f"\n✅ Successfully identified: {with_inchi}")
    print(f"⏭️  Skipped (cached): {skipped}")
    print(f"❌ Not found: {len(failed)}")
    
    # Save failed list
    failed_file = OUTPUT_DIR / "failed_inchi_lookups.json"
    with open(failed_file, 'w') as f:
        json.dump(failed, f, indent=2)
    
    # Create summary by class
    print("\n📊 Sample results:")
    samples = ['lsd', 'psilocybin', 'mdma', 'ketamine', 'cannabis', 'dmt']
    for sample in samples:
        if sample in results and results[sample].get('inchi'):
            r = results[sample]
            print(f"\n   {sample.upper()}:")
            print(f"      Formula: {r.get('molecular_formula')}")
            print(f"      InChI Key: {r.get('inchi_key')}")
            print(f"      IUPAC: {r.get('iupac_name', 'N/A')[:60]}...")
    
    # Create a simplified mapping file
    simple_mapping = {}
    for name, data in results.items():
        if data.get('inchi'):
            simple_mapping[name] = {
                'inchi': data['inchi'],
                'inchi_key': data['inchi_key'],
                'formula': data.get('molecular_formula'),
                'weight': data.get('molecular_weight'),
                'smiles': data.get('canonical_smiles'),
                'cid': data.get('pubchem_cid'),
            }
    
    simple_file = OUTPUT_DIR / "inchi_mapping.json"
    with open(simple_file, 'w') as f:
        json.dump(simple_mapping, f, indent=2)
    
    print(f"\n📁 Output files:")
    print(f"   {output_file}")
    print(f"   {simple_file}")
    print(f"   {failed_file}")
    
    # Summary stats
    summary = {
        'generated_at': datetime.now().isoformat(),
        'total_substances': len(substances),
        'with_inchi': with_inchi,
        'without_inchi': len(failed),
        'coverage_percent': round(100 * with_inchi / len(substances), 1),
        'failed_substances': failed[:50],
    }
    
    summary_file = OUTPUT_DIR / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✅ InChI coverage: {summary['coverage_percent']}%")


if __name__ == "__main__":
    fetch_all_inchi()
