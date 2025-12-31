#!/usr/bin/env python3
"""
Aggressive PDB Structure Preparation for MD Simulations
========================================================
Uses a multi-pass approach to get experimental structures working:
1. Extract target chain
2. Remove ALL incomplete residues (not just terminal)
3. Rebuild continuous fragments
4. Use OpenMM's force field to add missing hydrogens
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import numpy as np
import warnings
warnings.filterwarnings('ignore')


# Standard amino acids
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'
}

# Required atoms for each residue type (backbone + CB for most)
REQUIRED_ATOMS = {
    'ALA': {'N', 'CA', 'C', 'O', 'CB'},
    'ARG': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
    'ASN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'},
    'ASP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'},
    'CYS': {'N', 'CA', 'C', 'O', 'CB', 'SG'},
    'GLN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'},
    'GLU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
    'GLY': {'N', 'CA', 'C', 'O'},
    'HIS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HID': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HIE': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HIP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'ILE': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'},
    'LEU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'},
    'LYS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'},
    'MET': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'},
    'PHE': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'PRO': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'},
    'SER': {'N', 'CA', 'C', 'O', 'CB', 'OG'},
    'THR': {'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'},
    'TRP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'TYR': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
    'VAL': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'},
}

# Receptor definitions
RECEPTORS = {
    '6WHA': {'name': '5-HT2A', 'chain': 'A', 'desc': 'Serotonin 2A receptor'},
    '5TVN': {'name': '5-HT2B', 'chain': 'A', 'desc': 'Serotonin 2B receptor'},
    '6CM4': {'name': 'D2', 'chain': 'A', 'desc': 'Dopamine D2 receptor'},
    '5C1M': {'name': 'MOR', 'chain': 'A', 'desc': 'Mu opioid receptor'},
    '4DJH': {'name': 'KOR', 'chain': 'A', 'desc': 'Kappa opioid receptor'},
    '5TGZ': {'name': 'CB1', 'chain': 'A', 'desc': 'Cannabinoid receptor 1'},
    '4PE5': {'name': 'NMDA', 'chain': 'A', 'desc': 'NMDA receptor'},
}


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB from RCSB."""
    pdb_path = output_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"  Downloading {pdb_id}...")
        urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def parse_pdb(pdb_path: Path) -> dict:
    """Parse PDB file into structured format."""
    residues = {}  # (chain, resnum) -> {atoms, lines, resname}
    
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            atom_name = line[12:16].strip()
            altloc = line[16]
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = int(line[22:26].strip())
            
            # Skip alternate locations
            if altloc not in [' ', 'A']:
                continue
            
            # Skip non-standard residues
            if res_name not in STANDARD_AA:
                continue
            
            key = (chain, res_num)
            if key not in residues:
                residues[key] = {
                    'resname': res_name,
                    'atoms': set(),
                    'lines': []
                }
            
            residues[key]['atoms'].add(atom_name)
            # Remove altloc marker
            if altloc == 'A':
                line = line[:16] + ' ' + line[17:]
            residues[key]['lines'].append(line)
    
    return residues


def check_residue_complete(res_name: str, atoms: set) -> bool:
    """Check if a residue has all required heavy atoms."""
    # Normalize histidine
    if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
        required = REQUIRED_ATOMS.get('HIS', set())
    else:
        required = REQUIRED_ATOMS.get(res_name, set())
    
    if not required:
        return True
    
    return required.issubset(atoms)


def find_continuous_fragments(residues: dict, target_chain: str) -> list:
    """
    Find continuous stretches of complete residues.
    Returns list of fragments, each being a list of (chain, resnum) keys.
    """
    # Filter to target chain and check completeness
    chain_residues = {}
    for (chain, resnum), data in residues.items():
        if chain != target_chain:
            continue
        if check_residue_complete(data['resname'], data['atoms']):
            chain_residues[resnum] = data
    
    if not chain_residues:
        return []
    
    # Find continuous stretches
    sorted_nums = sorted(chain_residues.keys())
    fragments = []
    current_fragment = [sorted_nums[0]]
    
    for i in range(1, len(sorted_nums)):
        # Allow small gaps (1-2 residues) to be bridged
        if sorted_nums[i] - sorted_nums[i-1] <= 3:
            current_fragment.append(sorted_nums[i])
        else:
            if len(current_fragment) >= 10:  # Only keep fragments with 10+ residues
                fragments.append([(target_chain, n) for n in current_fragment])
            current_fragment = [sorted_nums[i]]
    
    # Don't forget last fragment
    if len(current_fragment) >= 10:
        fragments.append([(target_chain, n) for n in current_fragment])
    
    return fragments


def write_fragment_pdb(residues: dict, fragment: list, output_path: Path, 
                       fix_histidines: bool = True):
    """Write a single fragment to PDB file."""
    lines = []
    atom_num = 1
    
    for key in fragment:
        if key not in residues:
            continue
        
        data = residues[key]
        
        for line in data['lines']:
            # Renumber atoms
            new_line = f"ATOM  {atom_num:5d}" + line[11:]
            
            # Fix histidine naming
            if fix_histidines:
                new_line = new_line.replace(' HIS ', ' HID ')
            
            lines.append(new_line)
            atom_num += 1
    
    lines.append('TER\n')
    lines.append('END\n')
    
    with open(output_path, 'w') as f:
        f.writelines(lines)
    
    return len(fragment)


def prepare_with_openmm(input_pdb: Path, output_dir: Path, name: str) -> dict:
    """Use OpenMM to add hydrogens and validate the structure."""
    result = {
        'success': False,
        'input': str(input_pdb),
        'output_files': {}
    }
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load PDB
        pdb = app.PDBFile(str(input_pdb))
        initial_atoms = pdb.topology.getNumAtoms()
        initial_res = pdb.topology.getNumResidues()
        
        print(f"      Loaded: {initial_atoms} atoms, {initial_res} residues")
        
        # Set up force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Add hydrogens
        print(f"      Adding hydrogens...")
        modeller.addHydrogens(forcefield, pH=7.4)
        
        h_atoms = modeller.topology.getNumAtoms()
        print(f"      After H: {h_atoms} atoms")
        
        # Save protein-only
        protein_pdb = output_dir / f"{name}_prepared.pdb"
        with open(protein_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        result['output_files']['protein'] = str(protein_pdb)
        
        # Validate - create system
        print(f"      Validating...")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds
        )
        
        result['success'] = True
        result['atoms'] = h_atoms
        result['residues'] = modeller.topology.getNumResidues()
        
    except Exception as e:
        result['error'] = str(e)
        print(f"      Error: {e}")
    
    return result


def run_minimization(prepared_pdb: Path, output_dir: Path, 
                    steps: int = 1000) -> dict:
    """Run energy minimization with solvation."""
    result = {'success': False}
    
    try:
        # Load
        pdb = app.PDBFile(str(prepared_pdb))
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller and solvate
        modeller = app.Modeller(pdb.topology, pdb.positions)
        print(f"      Adding solvent...")
        modeller.addSolvent(forcefield, model='tip3p', padding=0.8*unit.nanometers)
        
        total_atoms = modeller.topology.getNumAtoms()
        print(f"      Total: {total_atoms} atoms")
        
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
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        init_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Initial E: {init_E:.1f} kJ/mol")
        
        # Minimize
        print(f"      Minimizing...")
        simulation.minimizeEnergy(maxIterations=steps)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Final E: {final_E:.1f} kJ/mol")
        
        # Save
        minimized_pdb = output_dir / 'minimized.pdb'
        with open(minimized_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        result['success'] = True
        result['initial_energy'] = init_E
        result['final_energy'] = final_E
        result['atoms'] = total_atoms
        result['output'] = str(minimized_pdb)
        
    except Exception as e:
        result['error'] = str(e)
        print(f"      Error: {e}")
    
    return result


def process_receptor(pdb_id: str, info: dict, structures_dir: Path, 
                    output_dir: Path) -> dict:
    """Process a single receptor structure."""
    result = {
        'pdb_id': pdb_id,
        'name': info['name'],
        'success': False,
        'fragments': []
    }
    
    print(f"\n{'='*60}")
    print(f"{pdb_id}: {info['name']} - {info['desc']}")
    print(f"{'='*60}")
    
    # Download
    raw_pdb = download_pdb(pdb_id, structures_dir)
    
    # Parse
    print(f"  Parsing structure...")
    residues = parse_pdb(raw_pdb)
    
    chain = info['chain']
    chain_residues = {k: v for k, v in residues.items() if k[0] == chain}
    print(f"  Chain {chain}: {len(chain_residues)} residues")
    
    # Check completeness
    complete = sum(1 for k, v in chain_residues.items() 
                   if check_residue_complete(v['resname'], v['atoms']))
    incomplete = len(chain_residues) - complete
    print(f"  Complete: {complete}, Incomplete: {incomplete}")
    
    # Find continuous fragments
    print(f"  Finding continuous fragments...")
    fragments = find_continuous_fragments(residues, chain)
    print(f"  Found {len(fragments)} fragments (≥10 residues)")
    
    if not fragments:
        result['error'] = 'No continuous fragments found'
        return result
    
    # Process each fragment
    receptor_dir = output_dir / info['name']
    receptor_dir.mkdir(parents=True, exist_ok=True)
    
    for i, fragment in enumerate(fragments):
        frag_name = f"{info['name']}_frag{i+1}"
        frag_dir = receptor_dir / f"fragment_{i+1}"
        frag_dir.mkdir(exist_ok=True)
        
        start_res = fragment[0][1]
        end_res = fragment[-1][1]
        print(f"\n  Fragment {i+1}: residues {start_res}-{end_res} ({len(fragment)} res)")
        
        # Write fragment
        frag_pdb = frag_dir / f"{frag_name}.pdb"
        write_fragment_pdb(residues, fragment, frag_pdb)
        
        # Prepare with OpenMM
        print(f"    Preparing with OpenMM...")
        prep_result = prepare_with_openmm(frag_pdb, frag_dir, frag_name)
        
        if prep_result['success']:
            print(f"    ✓ Preparation successful!")
            
            # Run minimization
            print(f"    Running minimization test...")
            min_result = run_minimization(
                Path(prep_result['output_files']['protein']),
                frag_dir,
                steps=500
            )
            
            if min_result['success']:
                print(f"    ✓ Minimization successful!")
                result['fragments'].append({
                    'fragment': i+1,
                    'residues': f"{start_res}-{end_res}",
                    'count': len(fragment),
                    'prepared_pdb': prep_result['output_files']['protein'],
                    'minimized_pdb': min_result['output'],
                    'energy': min_result['final_energy'],
                    'atoms': min_result['atoms']
                })
            else:
                print(f"    ✗ Minimization failed: {min_result.get('error', '')[:50]}")
        else:
            print(f"    ✗ Preparation failed: {prep_result.get('error', '')[:50]}")
    
    result['success'] = len(result['fragments']) > 0
    return result


def main():
    """Main entry point."""
    base_dir = Path(__file__).parent
    structures_dir = base_dir / 'structures'
    output_dir = base_dir / 'prepared_structures'
    
    structures_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("AGGRESSIVE PDB STRUCTURE PREPARATION")
    print("=" * 70)
    print("\nStrategy: Extract only COMPLETE residues, find continuous fragments")
    
    results = {}
    successful_receptors = []
    
    for pdb_id, info in RECEPTORS.items():
        result = process_receptor(pdb_id, info, structures_dir, output_dir)
        results[pdb_id] = result
        
        if result['success']:
            successful_receptors.append(result['name'])
    
    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    
    total_fragments = 0
    for pdb_id, result in results.items():
        name = result['name']
        if result['success']:
            frags = len(result['fragments'])
            total_fragments += frags
            print(f"  ✓ {name}: {frags} fragment(s) prepared")
            for frag in result['fragments']:
                print(f"      Fragment {frag['fragment']}: {frag['count']} residues, "
                      f"{frag['atoms']} atoms, E={frag['energy']:.1f} kJ/mol")
        else:
            print(f"  ✗ {name}: {result.get('error', 'Failed')}")
    
    print(f"\nTotal: {len(successful_receptors)}/{len(RECEPTORS)} receptors, "
          f"{total_fragments} usable fragments")
    
    # Save results
    results_file = base_dir / 'aggressive_preparation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_file}")
    
    # List output files
    print(f"\nPrepared structures in: {output_dir}")
    for d in sorted(output_dir.iterdir()):
        if d.is_dir():
            print(f"  {d.name}/")
    
    return results


if __name__ == '__main__':
    main()
