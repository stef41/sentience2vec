#!/usr/bin/env python3
"""
Consciousness Study - Receptor-Ligand Pose Prediction Pipeline
===============================================================

Workflow:
1. Chai-1: AI-based pose prediction (structure prediction)
2. PoseBusters: Validate predicted poses (quality check)
3. OpenFreeEnergy: Relative binding free energy calculations

"""

import os
import sys
import json
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

# Base directories
BASE_DIR = Path("/data/users/zacharie/whc/consciousness_study/md_simulations")
RECEPTORS_DIR = BASE_DIR / "receptors_full"
LIGANDS_DIR = BASE_DIR / "ligands_full"
OUTPUT_DIR = BASE_DIR / "docking_results"
OUTPUT_DIR.mkdir(exist_ok=True)


# =============================================================================
# STEP 1: CHAI-1 POSE PREDICTION
# =============================================================================

def run_chai1_prediction(receptor_pdb: Path, ligand_sdf: Path, output_dir: Path) -> Optional[Path]:
    """
    Use Chai-1 for protein-ligand pose prediction.
    
    Chai-1 is a multi-modal foundation model that can predict:
    - Protein structures
    - Protein-ligand complexes
    - Protein-nucleic acid complexes
    """
    from chai_lab.chai1 import run_inference
    import shutil
    
    # Read receptor sequence from PDB
    receptor_sequence = extract_sequence_from_pdb(receptor_pdb)
    if not receptor_sequence:
        print(f"    Could not extract sequence from {receptor_pdb}")
        return None
    
    # Read ligand SMILES from SDF
    ligand_smiles = extract_smiles_from_sdf(ligand_sdf)
    if not ligand_smiles:
        print(f"    Could not extract SMILES from {ligand_sdf}")
        return None
    
    # Create parent dir and put FASTA outside the output_dir
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    
    # Prepare FASTA input for Chai-1 (must be outside output_dir!)
    fasta_content = f">protein|name=receptor\n{receptor_sequence}\n>ligand|name=ligand\n{ligand_smiles}\n"
    
    fasta_path = output_dir.parent / f"{output_dir.name}_input.fasta"
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    # Clear output directory (Chai-1 requires completely empty directory)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Run Chai-1 inference
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=output_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device="cuda:0",  # Use GPU
            use_esm_embeddings=True,
        )
        
        # Get best pose (highest confidence)
        if candidates and len(candidates.cif_paths) > 0:
            best_pose = candidates.cif_paths[0]
            print(f"    ✓ Chai-1 prediction complete: {best_pose}")
            return best_pose
        
    except Exception as e:
        print(f"    ✗ Chai-1 error: {str(e)[:60]}")
    
    return None


def extract_sequence_from_pdb(pdb_path: Path) -> Optional[str]:
    """Extract amino acid sequence from PDB file."""
    
    AA_MAP = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C'
    }
    
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                chain = line[21]
                if chain == 'A' and res_name in AA_MAP:
                    residues[res_num] = AA_MAP[res_name]
    
    if not residues:
        return None
    
    sequence = ''.join(residues[k] for k in sorted(residues.keys()))
    return sequence


def extract_smiles_from_sdf(sdf_path: Path) -> Optional[str]:
    """Extract SMILES from SDF file using RDKit."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(sdf_path))
        if mol:
            return Chem.MolToSmiles(mol)
    except:
        pass
    return None


# =============================================================================
# STEP 2: POSEBUSTERS VALIDATION
# =============================================================================

def run_posebusters_validation(predicted_pose: Path, receptor_pdb: Path, 
                                ligand_sdf: Path) -> Dict:
    """
    Use PoseBusters to validate the predicted pose.
    
    PoseBusters checks for:
    - Molecular validity
    - Intramolecular validity (bond lengths, angles, planarity)
    - Intermolecular validity (clashes, interactions)
    """
    from posebusters import PoseBusters
    
    results = {
        'valid': False,
        'checks': {},
        'score': 0.0
    }
    
    try:
        # Initialize PoseBusters
        buster = PoseBusters(config="redock")
        
        # Run validation
        df = buster.bust(
            mol_pred=str(predicted_pose),
            mol_true=str(ligand_sdf),
            mol_cond=str(receptor_pdb)
        )
        
        # Extract results
        if not df.empty:
            row = df.iloc[0]
            
            # Get all check results
            check_columns = [c for c in df.columns if c not in ['mol_pred', 'mol_true', 'mol_cond']]
            for col in check_columns:
                results['checks'][col] = bool(row[col]) if col in row else None
            
            # Calculate overall validity
            bool_checks = [v for v in results['checks'].values() if isinstance(v, bool)]
            results['valid'] = all(bool_checks) if bool_checks else False
            results['score'] = sum(bool_checks) / len(bool_checks) if bool_checks else 0.0
            
            print(f"    ✓ PoseBusters: {results['score']*100:.1f}% checks passed")
        
    except Exception as e:
        print(f"    ✗ PoseBusters error: {str(e)[:60]}")
    
    return results


# =============================================================================
# STEP 3: FREE ENERGY CALCULATIONS (OpenMM-based FEP)
# =============================================================================

def setup_fep_calculation(receptor_pdb: Path, ligand_sdf: Path, 
                          reference_ligand_sdf: Path, output_dir: Path) -> Dict:
    """
    Set up Free Energy Perturbation calculation using OpenMM.
    
    This provides:
    - Atom mapping between ligands using RDKit MCS
    - Alchemical transformation setup
    - FEP protocol configuration
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    import numpy as np
    
    results = {
        'transformation_set': False,
        'atom_mapping': None,
        'protocol': 'OpenMM_FEP',
        'lambda_windows': 11
    }
    
    try:
        # Load molecules
        ligand_mol = Chem.MolFromMolFile(str(ligand_sdf))
        reference_mol = Chem.MolFromMolFile(str(reference_ligand_sdf))
        
        if ligand_mol is None or reference_mol is None:
            raise ValueError("Could not load molecules")
        
        # Find Maximum Common Substructure (MCS) for atom mapping
        mcs = rdFMCS.FindMCS(
            [ligand_mol, reference_mol],
            bondCompare=rdFMCS.BondCompare.CompareAny,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=60
        )
        
        if mcs.numAtoms > 0:
            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            
            # Get atom mappings
            match_lig = ligand_mol.GetSubstructMatch(mcs_mol)
            match_ref = reference_mol.GetSubstructMatch(mcs_mol)
            
            atom_mapping = {match_lig[i]: match_ref[i] for i in range(len(match_lig))}
            
            results['atom_mapping'] = {
                'n_mapped_atoms': len(atom_mapping),
                'mcs_smarts': mcs.smartsString,
                'ligand_atoms': ligand_mol.GetNumAtoms(),
                'reference_atoms': reference_mol.GetNumAtoms(),
                'mapping': atom_mapping
            }
            
            # Define lambda schedule for FEP
            results['lambda_schedule'] = {
                'electrostatics': np.linspace(0, 1, 11).tolist(),
                'sterics': np.linspace(0, 1, 11).tolist()
            }
            
            # Save FEP setup
            output_dir.mkdir(parents=True, exist_ok=True)
            setup_file = output_dir / "fep_setup.json"
            
            fep_config = {
                'receptor': str(receptor_pdb),
                'ligand_A': str(ligand_sdf),
                'ligand_B': str(reference_ligand_sdf),
                'atom_mapping': results['atom_mapping'],
                'lambda_windows': results['lambda_windows'],
                'lambda_schedule': results['lambda_schedule'],
                'simulation_settings': {
                    'temperature_K': 300,
                    'pressure_atm': 1.0,
                    'timestep_fs': 2.0,
                    'equilibration_steps': 50000,
                    'production_steps': 500000,
                    'n_replicas': 3
                }
            }
            
            with open(setup_file, 'w') as f:
                json.dump(fep_config, f, indent=2, default=str)
            
            results['transformation_set'] = True
            print(f"    ✓ FEP Setup: {results['atom_mapping']['n_mapped_atoms']} mapped atoms")
            print(f"       MCS: {mcs.smartsString[:50]}...")
        else:
            print(f"    ✗ No common substructure found")
        
    except Exception as e:
        print(f"    ✗ FEP setup error: {str(e)[:60]}")
    
    return results


# =============================================================================
# FULL PIPELINE
# =============================================================================

def run_full_pipeline(receptor_name: str, ligand_name: str, 
                      reference_ligand: str = None) -> Dict:
    """
    Run the complete Chai-1 → PoseBusters → FEP pipeline.
    """
    print(f"\n{'='*60}")
    print(f"PIPELINE: {receptor_name} + {ligand_name}")
    print(f"{'='*60}")
    
    results = {
        'receptor': receptor_name,
        'ligand': ligand_name,
        'chai1': None,
        'posebusters': None,
        'fep': None,
        'success': False
    }
    
    # Find files
    receptor_pdb = RECEPTORS_DIR / receptor_name / f"{receptor_name}_prepared.pdb"
    ligand_sdf = LIGANDS_DIR / ligand_name / f"{ligand_name}.sdf"
    
    if not receptor_pdb.exists():
        print(f"  ✗ Receptor not found: {receptor_pdb}")
        return results
    
    if not ligand_sdf.exists():
        print(f"  ✗ Ligand not found: {ligand_sdf}")
        return results
    
    output_dir = OUTPUT_DIR / f"{receptor_name}_{ligand_name}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Chai-1 pose prediction
    print("\n[1/3] Chai-1 Pose Prediction...")
    chai_output = output_dir / "chai1"
    predicted_pose = run_chai1_prediction(receptor_pdb, ligand_sdf, chai_output)
    results['chai1'] = {'predicted_pose': str(predicted_pose) if predicted_pose else None}
    
    # Step 2: PoseBusters validation
    if predicted_pose and predicted_pose.exists():
        print("\n[2/3] PoseBusters Validation...")
        pb_results = run_posebusters_validation(predicted_pose, receptor_pdb, ligand_sdf)
        results['posebusters'] = pb_results
    else:
        print("\n[2/3] PoseBusters Validation... SKIPPED (no pose)")
    
    # Step 3: FEP setup (if reference ligand provided)
    if reference_ligand:
        ref_sdf = LIGANDS_DIR / reference_ligand / f"{reference_ligand}.sdf"
        if ref_sdf.exists():
            print("\n[3/3] FEP Free Energy Setup...")
            fep_output = output_dir / "fep"
            fep_results = setup_fep_calculation(receptor_pdb, ligand_sdf, ref_sdf, fep_output)
            results['fep'] = fep_results
        else:
            print(f"\n[3/3] FEP... SKIPPED (reference {reference_ligand} not found)")
    else:
        print("\n[3/3] FEP... SKIPPED (no reference ligand)")
    
    # Overall success
    results['success'] = (
        results['chai1'] is not None and 
        results['chai1'].get('predicted_pose') is not None
    )
    
    # Save results
    with open(output_dir / "pipeline_results.json", 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE: {'SUCCESS' if results['success'] else 'FAILED'}")
    print(f"Results saved to: {output_dir}")
    
    return results


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def run_batch_pipeline(pairs: List[Tuple[str, str, Optional[str]]]):
    """
    Run pipeline for multiple receptor-ligand pairs.
    
    Args:
        pairs: List of (receptor_name, ligand_name, reference_ligand) tuples
    """
    all_results = []
    
    for receptor, ligand, reference in pairs:
        try:
            result = run_full_pipeline(receptor, ligand, reference)
            all_results.append(result)
        except Exception as e:
            print(f"Error processing {receptor}+{ligand}: {e}")
            all_results.append({
                'receptor': receptor,
                'ligand': ligand,
                'error': str(e),
                'success': False
            })
    
    # Save batch results
    batch_file = OUTPUT_DIR / "batch_results.json"
    with open(batch_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Summary
    success_count = sum(1 for r in all_results if r.get('success'))
    print(f"\n{'='*60}")
    print(f"BATCH COMPLETE: {success_count}/{len(all_results)} successful")
    print(f"Results saved to: {batch_file}")
    
    return all_results


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Example: Key psychedelic-receptor pairs
    PRIORITY_PAIRS = [
        # Serotonin receptors + psychedelics
        ("5-HT2A", "LSD", "Serotonin"),
        ("5-HT2A", "DMT", "Serotonin"),
        ("5-HT2A", "Psilocin", "Serotonin"),
        ("5-HT2A", "Mescaline", "Serotonin"),
        ("5-HT1A", "LSD", "Serotonin"),
        ("5-HT1A", "DMT", "Serotonin"),
        
        # Dissociatives + NMDA/Sigma
        ("NMDA_GluN2B", "Ketamine", "Glutamate"),
        ("Sigma1", "Ketamine", None),
        ("Sigma1", "DMT", None),
        
        # Opioids
        ("MOR", "Morphine", None),
        ("KOR", "Salvinorin_A", None),
        
        # Cannabinoids
        ("CB1", "Tetrahydrocannabinol", None),
        
        # Dopamine + stimulants
        ("D2", "Amphetamine", "Dopamine"),
        ("DAT", "Cocaine", "Dopamine"),
        ("SERT", "MDMA", "Serotonin"),
    ]
    
    print("="*60)
    print("CHAI-1 → POSEBUSTERS → OPENFE PIPELINE")
    print("="*60)
    print(f"\nReceptors: {RECEPTORS_DIR}")
    print(f"Ligands: {LIGANDS_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"\nPriority pairs: {len(PRIORITY_PAIRS)}")
    
    # Run single example first
    if len(sys.argv) > 1 and sys.argv[1] == "--batch":
        run_batch_pipeline(PRIORITY_PAIRS)
    else:
        # Test with single pair
        run_full_pipeline("5-HT2A", "DMT", "Serotonin")
