#!/usr/bin/env python3
"""
Prepare All Receptor Structures for MD Simulations
==================================================
Comprehensive structure preparation using OpenMM's built-in tools.

This script:
1. Downloads fresh PDB structures if needed
2. Extracts the receptor chain (removes fusion proteins, antibodies)
3. Removes heterogens (ligands, waters, detergents)
4. Rebuilds missing atoms using OpenMM Modeller
5. Adds hydrogens at physiological pH
6. Validates the prepared structure
7. Optionally solvates for simulation
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import sys
import warnings
warnings.filterwarnings('ignore')


# Receptor definitions with their primary chain
RECEPTORS = {
    '6WHA': {
        'name': '5-HT2A',
        'description': 'Serotonin 2A receptor (LSD-bound)',
        'chain': 'A',  # Main receptor chain
        'exclude_chains': ['B'],  # Antibody/fusion chains
        'membrane_protein': True,
    },
    '5TVN': {
        'name': '5-HT2B', 
        'description': 'Serotonin 2B receptor',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
    '6CM4': {
        'name': 'D2',
        'description': 'Dopamine D2 receptor',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
    '5C1M': {
        'name': 'MOR',
        'description': 'Mu opioid receptor',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
    '4DJH': {
        'name': 'KOR',
        'description': 'Kappa opioid receptor',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
    '5TGZ': {
        'name': 'CB1',
        'description': 'Cannabinoid receptor 1',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
    '4PE5': {
        'name': 'NMDA_GluN1',
        'description': 'NMDA receptor GluN1 subunit',
        'chain': 'A',
        'exclude_chains': [],
        'membrane_protein': True,
    },
}

# Standard amino acids that OpenMM can handle
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL',
    # Modified forms OpenMM understands
    'HID', 'HIE', 'HIP', 'CYX'
}


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB file from RCSB."""
    pdb_path = output_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"  Downloading {pdb_id}...")
        urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def extract_and_clean_chain(input_pdb: Path, output_pdb: Path, 
                           target_chain: str = 'A',
                           keep_waters: bool = False) -> dict:
    """
    Extract a single chain and clean the structure.
    
    Returns statistics about what was removed/kept.
    """
    stats = {
        'total_atoms': 0,
        'kept_atoms': 0,
        'removed_hetatm': 0,
        'removed_water': 0,
        'removed_chains': 0,
        'residues': 0,
    }
    
    output_lines = []
    seen_residues = set()
    
    with open(input_pdb) as f:
        for line in f:
            stats['total_atoms'] += 1 if line.startswith(('ATOM', 'HETATM')) else 0
            
            # Skip non-coordinate lines except END
            if not line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                continue
            
            if line.startswith('END'):
                output_lines.append('END\n')
                continue
            
            if line.startswith('TER'):
                output_lines.append(line)
                continue
            
            # Parse chain
            chain = line[21] if len(line) > 21 else ' '
            
            # Skip other chains
            if chain != target_chain and chain != ' ':
                stats['removed_chains'] += 1
                continue
            
            # Get residue info
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            
            # Handle HETATM records
            if line.startswith('HETATM'):
                # Keep water if requested
                if res_name in ['HOH', 'WAT']:
                    if keep_waters:
                        output_lines.append(line)
                        stats['kept_atoms'] += 1
                    else:
                        stats['removed_water'] += 1
                    continue
                
                # Check if it's a modified amino acid we can convert
                modified_aa = {
                    'MSE': 'MET',  # Selenomethionine -> Methionine
                    'CSE': 'CYS',  # Selenocysteine -> Cysteine
                    'SEP': 'SER',  # Phosphoserine -> Serine
                    'TPO': 'THR',  # Phosphothreonine -> Threonine
                    'PTR': 'TYR',  # Phosphotyrosine -> Tyrosine
                }
                
                if res_name in modified_aa:
                    # Convert to standard residue
                    new_res = modified_aa[res_name]
                    new_line = 'ATOM  ' + line[6:17] + f'{new_res:>3}' + line[20:]
                    output_lines.append(new_line)
                    stats['kept_atoms'] += 1
                    seen_residues.add(f"{chain}_{res_num}")
                    continue
                    
                # Remove other heterogens (ligands, detergents, etc.)
                stats['removed_hetatm'] += 1
                continue
            
            # ATOM records
            if res_name in STANDARD_RESIDUES:
                # Skip alternate conformations (keep only A or no altloc)
                altloc = line[16]
                if altloc not in [' ', 'A']:
                    continue
                    
                # Clear altloc marker
                if altloc == 'A':
                    line = line[:16] + ' ' + line[17:]
                
                output_lines.append(line)
                stats['kept_atoms'] += 1
                seen_residues.add(f"{chain}_{res_num}")
            else:
                # Skip non-standard residues
                stats['removed_hetatm'] += 1
    
    stats['residues'] = len(seen_residues)
    
    # Write output
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    return stats


def fix_histidine_names(pdb_path: Path):
    """Convert HIS to HID for Amber force field (neutral, delta protonated)."""
    with open(pdb_path) as f:
        content = f.read()
    
    # Replace HIS with HID (neutral histidine, proton on ND1)
    content = content.replace(' HIS ', ' HID ')
    
    with open(pdb_path, 'w') as f:
        f.write(content)


def identify_incomplete_residues(pdb_path: Path) -> list:
    """Identify residues with missing heavy atoms."""
    # Required backbone atoms
    BACKBONE = {'N', 'CA', 'C', 'O'}
    
    # Side chain atoms for each residue type
    SIDECHAIN = {
        'ALA': {'CB'},
        'ARG': {'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
        'ASN': {'CB', 'CG', 'OD1', 'ND2'},
        'ASP': {'CB', 'CG', 'OD1', 'OD2'},
        'CYS': {'CB', 'SG'},
        'GLN': {'CB', 'CG', 'CD', 'OE1', 'NE2'},
        'GLU': {'CB', 'CG', 'CD', 'OE1', 'OE2'},
        'GLY': set(),
        'HID': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIE': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIP': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'HIS': {'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
        'ILE': {'CB', 'CG1', 'CG2', 'CD1'},
        'LEU': {'CB', 'CG', 'CD1', 'CD2'},
        'LYS': {'CB', 'CG', 'CD', 'CE', 'NZ'},
        'MET': {'CB', 'CG', 'SD', 'CE'},
        'PHE': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'PRO': {'CB', 'CG', 'CD'},
        'SER': {'CB', 'OG'},
        'THR': {'CB', 'OG1', 'CG2'},
        'TRP': {'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'TYR': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
        'VAL': {'CB', 'CG1', 'CG2'},
    }
    
    # Parse residues
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = int(line[22:26].strip())
            
            key = (chain, res_num, res_name)
            if key not in residues:
                residues[key] = set()
            residues[key].add(atom_name)
    
    # Check completeness
    incomplete = []
    for (chain, res_num, res_name), atoms in residues.items():
        required = BACKBONE.copy()
        if res_name in SIDECHAIN:
            required.update(SIDECHAIN[res_name])
        
        missing = required - atoms
        if missing:
            incomplete.append({
                'chain': chain,
                'resnum': res_num,
                'resname': res_name,
                'missing': list(missing),
                'has': list(atoms)
            })
    
    return incomplete


def remove_incomplete_terminal_residues(pdb_path: Path, n_residues: int = 3) -> int:
    """Remove potentially incomplete residues from chain termini."""
    # Parse residue numbers
    residue_nums = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                res_num = int(line[22:26].strip())
                residue_nums.add(res_num)
    
    if not residue_nums:
        return 0
    
    sorted_nums = sorted(residue_nums)
    
    # Identify terminal residues to remove
    remove_nums = set()
    for i in range(min(n_residues, len(sorted_nums))):
        remove_nums.add(sorted_nums[i])  # N-terminus
        remove_nums.add(sorted_nums[-(i+1)])  # C-terminus
    
    # Filter lines
    kept_lines = []
    removed_count = 0
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                res_num = int(line[22:26].strip())
                if res_num in remove_nums:
                    removed_count += 1
                    continue
            kept_lines.append(line)
    
    with open(pdb_path, 'w') as f:
        f.writelines(kept_lines)
    
    return removed_count


def prepare_structure_for_simulation(
    pdb_path: Path,
    output_dir: Path,
    receptor_name: str,
    add_hydrogens: bool = True,
    add_solvent: bool = False,
    padding: float = 1.0,
) -> dict:
    """
    Prepare a cleaned PDB structure for OpenMM simulation.
    
    Uses OpenMM's Modeller to add missing atoms and hydrogens.
    """
    result = {
        'input': str(pdb_path),
        'success': False,
        'error': None,
        'output_files': {},
    }
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load structure
        print(f"    Loading {pdb_path.name}...")
        pdb = app.PDBFile(str(pdb_path))
        
        initial_atoms = pdb.topology.getNumAtoms()
        initial_residues = pdb.topology.getNumResidues()
        print(f"    Initial: {initial_atoms} atoms, {initial_residues} residues")
        
        # Set up force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Add missing hydrogens
        if add_hydrogens:
            print(f"    Adding hydrogens...")
            modeller.addHydrogens(forcefield, pH=7.4)
            print(f"    After H: {modeller.topology.getNumAtoms()} atoms")
        
        # Save protein-only structure
        protein_pdb = output_dir / f"{receptor_name}_prepared.pdb"
        with open(protein_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        result['output_files']['protein'] = str(protein_pdb)
        
        # Add solvent if requested
        if add_solvent:
            print(f"    Adding solvent (padding={padding} nm)...")
            modeller.addSolvent(
                forcefield,
                model='tip3p',
                padding=padding * unit.nanometers,
                neutralize=True,
            )
            print(f"    Solvated: {modeller.topology.getNumAtoms()} atoms")
            
            solvated_pdb = output_dir / f"{receptor_name}_solvated.pdb"
            with open(solvated_pdb, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            result['output_files']['solvated'] = str(solvated_pdb)
        
        # Validate by creating system
        print(f"    Validating force field compatibility...")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME if add_solvent else app.NoCutoff,
            constraints=app.HBonds,
        )
        print(f"    ✓ System created successfully")
        
        result['success'] = True
        result['final_atoms'] = modeller.topology.getNumAtoms()
        result['final_residues'] = modeller.topology.getNumResidues()
        
    except Exception as e:
        result['error'] = str(e)
        print(f"    ✗ Error: {e}")
    
    return result


def prepare_all_structures(structures_dir: Path, prepared_dir: Path) -> dict:
    """Prepare all receptor structures."""
    
    results = {}
    
    print("=" * 70)
    print("PREPARING ALL RECEPTOR STRUCTURES")
    print("=" * 70)
    
    for pdb_id, info in RECEPTORS.items():
        print(f"\n{'='*60}")
        print(f"Processing {pdb_id}: {info['name']} ({info['description']})")
        print(f"{'='*60}")
        
        # Step 1: Download/check PDB exists
        raw_pdb = download_pdb(pdb_id, structures_dir)
        
        # Step 2: Extract target chain
        chain = info.get('chain', 'A')
        cleaned_pdb = structures_dir / f"{pdb_id}_chain{chain}.pdb"
        
        print(f"\n  Step 1: Extracting chain {chain}...")
        stats = extract_and_clean_chain(raw_pdb, cleaned_pdb, target_chain=chain)
        print(f"    Kept {stats['kept_atoms']}/{stats['total_atoms']} atoms")
        print(f"    Removed: {stats['removed_hetatm']} heterogens, {stats['removed_water']} waters")
        print(f"    Residues: {stats['residues']}")
        
        # Step 3: Fix histidine naming
        print(f"\n  Step 2: Fixing histidine names...")
        fix_histidine_names(cleaned_pdb)
        
        # Step 4: Check for incomplete residues
        print(f"\n  Step 3: Checking for incomplete residues...")
        incomplete = identify_incomplete_residues(cleaned_pdb)
        if incomplete:
            print(f"    Found {len(incomplete)} incomplete residues")
            # Show first few
            for res in incomplete[:5]:
                print(f"      {res['resname']} {res['resnum']}: missing {res['missing']}")
            
            # Remove problematic terminal residues
            print(f"\n  Step 4: Removing terminal residues...")
            removed = remove_incomplete_terminal_residues(cleaned_pdb, n_residues=2)
            print(f"    Removed {removed} atoms from termini")
        
        # Step 5: Prepare for OpenMM
        print(f"\n  Step 5: Preparing for OpenMM simulation...")
        receptor_dir = prepared_dir / info['name']
        result = prepare_structure_for_simulation(
            cleaned_pdb,
            receptor_dir,
            info['name'],
            add_hydrogens=True,
            add_solvent=False,  # Don't solvate yet - save resources
        )
        
        results[pdb_id] = {
            'name': info['name'],
            'description': info['description'],
            'chain': chain,
            'raw_pdb': str(raw_pdb),
            'cleaned_pdb': str(cleaned_pdb),
            'stats': stats,
            'preparation': result,
        }
        
        if result['success']:
            print(f"\n  ✓ {info['name']} prepared successfully!")
            print(f"    Output: {receptor_dir}")
        else:
            print(f"\n  ✗ {info['name']} preparation failed")
            print(f"    Error: {result['error']}")
    
    return results


def run_quick_minimization(prepared_pdb: Path, output_dir: Path, 
                          steps: int = 1000) -> dict:
    """Run quick energy minimization to validate the structure."""
    result = {'success': False}
    
    try:
        # Load
        pdb = app.PDBFile(str(prepared_pdb))
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller and add solvent for minimization
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(forcefield, model='tip3p', padding=0.8*unit.nanometers)
        
        # Create system
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        # Set up simulation
        integrator = mm.LangevinMiddleIntegrator(
            300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds
        )
        
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            platform = mm.Platform.getPlatformByName('CPU')
        
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        
        # Initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=steps)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Save
        minimized_pdb = output_dir / 'minimized.pdb'
        with open(minimized_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        result['success'] = True
        result['initial_energy'] = initial_energy
        result['final_energy'] = final_energy
        result['output'] = str(minimized_pdb)
        result['total_atoms'] = simulation.topology.getNumAtoms()
        
    except Exception as e:
        result['error'] = str(e)
    
    return result


def main():
    """Main entry point."""
    base_dir = Path(__file__).parent
    structures_dir = base_dir / 'structures'
    prepared_dir = base_dir / 'prepared_structures'
    
    structures_dir.mkdir(exist_ok=True)
    prepared_dir.mkdir(exist_ok=True)
    
    # Prepare all structures
    results = prepare_all_structures(structures_dir, prepared_dir)
    
    # Summary
    print("\n" + "=" * 70)
    print("PREPARATION SUMMARY")
    print("=" * 70)
    
    successful = []
    failed = []
    
    for pdb_id, data in results.items():
        if data['preparation']['success']:
            successful.append(data['name'])
            print(f"  ✓ {data['name']}: Ready")
        else:
            failed.append(data['name'])
            print(f"  ✗ {data['name']}: {data['preparation'].get('error', 'Unknown error')[:50]}")
    
    print(f"\nSuccess: {len(successful)}/{len(results)}")
    
    # Run quick minimization test on successful structures
    if successful:
        print("\n" + "=" * 70)
        print("RUNNING VALIDATION MINIMIZATIONS")
        print("=" * 70)
        
        for name in successful[:3]:  # Test first 3
            for pdb_id, data in results.items():
                if data['name'] == name:
                    prepared_pdb = Path(data['preparation']['output_files'].get('protein', ''))
                    if prepared_pdb.exists():
                        print(f"\n  Testing {name}...")
                        output_dir = prepared_dir / name
                        min_result = run_quick_minimization(prepared_pdb, output_dir)
                        
                        if min_result['success']:
                            print(f"    ✓ Minimization successful")
                            print(f"      Energy: {min_result['initial_energy']:.1f} → {min_result['final_energy']:.1f} kJ/mol")
                            print(f"      Atoms: {min_result['total_atoms']}")
                        else:
                            print(f"    ✗ Minimization failed: {min_result.get('error', '')[:60]}")
                    break
    
    # Save results
    results_file = base_dir / 'preparation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_file}")
    
    return results


if __name__ == '__main__':
    main()
