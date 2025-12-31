#!/usr/bin/env python3
"""
PDB Structure Preparation with Terminal Capping
================================================
Adds proper terminal capping groups for OpenMM compatibility.

The issue is that OpenMM's Amber14 force field requires either:
1. Standard N-terminal (with NH3+) and C-terminal (with COO-)
2. Or explicit capping groups (ACE for N-term, NME for C-term)

This script uses OpenMM's Modeller.addHydrogens() with proper terminal handling.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import numpy as np
from io import StringIO
import warnings
warnings.filterwarnings('ignore')

# Standard amino acids
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'
}

# Required heavy atoms for complete residues
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

RECEPTORS = {
    '6WHA': {'name': '5-HT2A', 'chain': 'A', 'desc': 'Serotonin 2A receptor'},
    '5TVN': {'name': '5-HT2B', 'chain': 'A', 'desc': 'Serotonin 2B receptor'},
    '6CM4': {'name': 'D2', 'chain': 'A', 'desc': 'Dopamine D2 receptor'},
    '5C1M': {'name': 'MOR', 'chain': 'A', 'desc': 'Mu opioid receptor'},
    '4DJH': {'name': 'KOR', 'chain': 'A', 'desc': 'Kappa opioid receptor'},
    '5TGZ': {'name': 'CB1', 'chain': 'A', 'desc': 'Cannabinoid receptor 1'},
}


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB from RCSB."""
    pdb_path = output_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"  Downloading {pdb_id}...")
        urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def check_residue_complete(res_name: str, atoms: set) -> bool:
    """Check if a residue has all required heavy atoms."""
    if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
        required = REQUIRED_ATOMS.get('HIS', set())
    else:
        required = REQUIRED_ATOMS.get(res_name, set())
    return required.issubset(atoms) if required else True


def parse_and_extract_complete_residues(pdb_path: Path, target_chain: str) -> dict:
    """
    Parse PDB and extract only complete residues from target chain.
    Returns dict with residue data.
    """
    residues = {}
    
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            atom_name = line[12:16].strip()
            altloc = line[16]
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = int(line[22:26].strip())
            
            # Filter
            if chain != target_chain:
                continue
            if altloc not in [' ', 'A']:
                continue
            if res_name not in STANDARD_AA:
                continue
            
            key = res_num
            if key not in residues:
                residues[key] = {'resname': res_name, 'atoms': set(), 'lines': []}
            
            residues[key]['atoms'].add(atom_name)
            # Clean the line
            clean_line = line[:16] + ' ' + line[17:]  # Remove altloc
            residues[key]['lines'].append(clean_line)
    
    # Filter to only complete residues
    complete = {}
    for res_num, data in residues.items():
        if check_residue_complete(data['resname'], data['atoms']):
            complete[res_num] = data
    
    return complete


def find_longest_continuous_stretch(residues: dict, min_gap: int = 2) -> list:
    """
    Find the longest continuous stretch of residues.
    Allows small gaps (1-2 residues) for flexibility.
    """
    if not residues:
        return []
    
    sorted_nums = sorted(residues.keys())
    
    # Build stretches
    stretches = []
    current = [sorted_nums[0]]
    
    for i in range(1, len(sorted_nums)):
        gap = sorted_nums[i] - sorted_nums[i-1]
        if gap <= min_gap + 1:  # Allow gap
            current.append(sorted_nums[i])
        else:
            stretches.append(current)
            current = [sorted_nums[i]]
    stretches.append(current)
    
    # Return longest
    longest = max(stretches, key=len)
    
    # Filter to only include residues we have
    return [n for n in longest if n in residues]


def add_terminal_oxt(lines: list) -> list:
    """
    Add OXT atom to the C-terminal residue.
    This creates a proper carboxylate terminus.
    """
    # Find last residue's C and O atoms
    last_c_line = None
    last_o_line = None
    last_res_num = None
    
    for line in reversed(lines):
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            res_num = int(line[22:26])
            
            if last_res_num is None:
                last_res_num = res_num
            
            if res_num != last_res_num:
                break
                
            if atom_name == 'C':
                last_c_line = line
            elif atom_name == 'O':
                last_o_line = line
    
    if last_c_line and last_o_line:
        # Get coordinates of C and O
        c_x = float(last_c_line[30:38])
        c_y = float(last_c_line[38:46])
        c_z = float(last_c_line[46:54])
        
        o_x = float(last_o_line[30:38])
        o_y = float(last_o_line[38:46])
        o_z = float(last_o_line[46:54])
        
        # OXT is roughly opposite O relative to C
        oxt_x = 2 * c_x - o_x
        oxt_y = 2 * c_y - o_y
        oxt_z = 2 * c_z - o_z
        
        # Create OXT line based on O line template
        atom_num = len([l for l in lines if l.startswith('ATOM')]) + 1
        oxt_line = f"ATOM  {atom_num:5d}  OXT {last_o_line[17:26]}" + \
                   f"{oxt_x:8.3f}{oxt_y:8.3f}{oxt_z:8.3f}" + \
                   last_o_line[54:]
        
        # Insert before TER/END
        new_lines = []
        for line in lines:
            if line.startswith('TER') or line.startswith('END'):
                new_lines.append(oxt_line)
            new_lines.append(line)
        
        return new_lines
    
    return lines


def write_clean_pdb(residues: dict, res_nums: list, output_path: Path,
                    fix_histidines: bool = True, add_oxt: bool = True) -> int:
    """Write clean PDB with proper formatting."""
    lines = []
    atom_num = 1
    new_res_num = 1  # Renumber from 1
    
    for orig_num in res_nums:
        if orig_num not in residues:
            continue
        
        data = residues[orig_num]
        
        for line in data['lines']:
            # Renumber atoms and residues
            new_line = f"ATOM  {atom_num:5d}" + line[11:22]
            new_line += f"{new_res_num:4d}" + line[26:]
            
            # Fix histidine naming for Amber
            if fix_histidines:
                new_line = new_line.replace(' HIS ', ' HID ')
            
            lines.append(new_line)
            atom_num += 1
        
        new_res_num += 1
    
    lines.append('TER\n')
    lines.append('END\n')
    
    # Add OXT to C-terminus
    if add_oxt:
        lines = add_terminal_oxt(lines)
    
    with open(output_path, 'w') as f:
        f.writelines(lines)
    
    return len(res_nums)


def prepare_with_openmm(input_pdb: Path, output_dir: Path, name: str,
                        solvate: bool = True) -> dict:
    """
    Prepare structure using OpenMM with proper terminal handling.
    """
    result = {
        'success': False,
        'input': str(input_pdb),
        'output_files': {}
    }
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load
        pdb = app.PDBFile(str(input_pdb))
        initial_atoms = pdb.topology.getNumAtoms()
        initial_res = pdb.topology.getNumResidues()
        
        print(f"    Loaded: {initial_atoms} atoms, {initial_res} residues")
        
        # Force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Add hydrogens - this should now work with proper termini
        print(f"    Adding hydrogens...")
        modeller.addHydrogens(forcefield, pH=7.4)
        
        print(f"    After H: {modeller.topology.getNumAtoms()} atoms")
        
        # Save protein
        protein_pdb = output_dir / f"{name}_hydrogenated.pdb"
        with open(protein_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        result['output_files']['protein'] = str(protein_pdb)
        
        if solvate:
            print(f"    Adding solvent...")
            modeller.addSolvent(
                forcefield, model='tip3p',
                padding=1.0*unit.nanometers,
                neutralize=True
            )
            print(f"    Solvated: {modeller.topology.getNumAtoms()} atoms")
            
            solvated_pdb = output_dir / f"{name}_solvated.pdb"
            with open(solvated_pdb, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            result['output_files']['solvated'] = str(solvated_pdb)
        
        # Validate by creating system
        print(f"    Creating system...")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME if solvate else app.NoCutoff,
            nonbondedCutoff=1.0*unit.nanometers if solvate else None,
            constraints=app.HBonds
        )
        
        result['success'] = True
        result['atoms'] = modeller.topology.getNumAtoms()
        result['residues'] = modeller.topology.getNumResidues()
        
        print(f"    ✓ System validated!")
        
    except Exception as e:
        result['error'] = str(e)
        print(f"    ✗ Error: {str(e)[:80]}")
    
    return result


def run_minimization(solvated_pdb: Path, output_dir: Path,
                     min_steps: int = 1000, eq_steps: int = 5000) -> dict:
    """Run energy minimization and brief equilibration."""
    result = {'success': False}
    
    try:
        # Load
        pdb = app.PDBFile(str(solvated_pdb))
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # System
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        # Integrator
        integrator = mm.LangevinMiddleIntegrator(
            310*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds
        )
        
        # Platform
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            platform = mm.Platform.getPlatformByName('CPU')
        
        # Simulation
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        
        # Initial energy
        state = simulation.context.getState(getEnergy=True)
        init_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Initial E: {init_E:.1f} kJ/mol")
        
        # Minimize
        print(f"      Minimizing ({min_steps} steps)...")
        simulation.minimizeEnergy(maxIterations=min_steps)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        min_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Minimized E: {min_E:.1f} kJ/mol")
        
        # Save minimized
        minimized_pdb = output_dir / 'minimized.pdb'
        with open(minimized_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Brief equilibration
        print(f"      Equilibrating ({eq_steps} steps = {eq_steps*2/1000:.1f} ps)...")
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
        simulation.step(eq_steps)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        eq_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Equilibrated E: {eq_E:.1f} kJ/mol")
        
        # Save equilibrated
        equilibrated_pdb = output_dir / 'equilibrated.pdb'
        with open(equilibrated_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        result['success'] = True
        result['initial_energy'] = init_E
        result['minimized_energy'] = min_E
        result['equilibrated_energy'] = eq_E
        result['minimized_pdb'] = str(minimized_pdb)
        result['equilibrated_pdb'] = str(equilibrated_pdb)
        result['atoms'] = simulation.topology.getNumAtoms()
        
    except Exception as e:
        result['error'] = str(e)
        print(f"      ✗ Error: {str(e)[:80]}")
    
    return result


def process_receptor(pdb_id: str, info: dict, structures_dir: Path,
                    output_dir: Path) -> dict:
    """Process a single receptor."""
    result = {
        'pdb_id': pdb_id,
        'name': info['name'],
        'description': info['desc'],
        'success': False
    }
    
    print(f"\n{'='*70}")
    print(f" {pdb_id}: {info['name']} - {info['desc']}")
    print(f"{'='*70}")
    
    # Download
    raw_pdb = download_pdb(pdb_id, structures_dir)
    
    # Parse and extract complete residues
    print(f"  Extracting chain {info['chain']}...")
    residues = parse_and_extract_complete_residues(raw_pdb, info['chain'])
    print(f"  Complete residues: {len(residues)}")
    
    if len(residues) < 20:
        result['error'] = f"Only {len(residues)} complete residues"
        return result
    
    # Find longest stretch
    print(f"  Finding continuous stretch...")
    stretch = find_longest_continuous_stretch(residues)
    print(f"  Longest stretch: {len(stretch)} residues (res {min(stretch)}-{max(stretch)})")
    
    if len(stretch) < 20:
        result['error'] = f"Longest stretch only {len(stretch)} residues"
        return result
    
    # Write clean PDB
    receptor_dir = output_dir / info['name']
    receptor_dir.mkdir(parents=True, exist_ok=True)
    
    clean_pdb = receptor_dir / f"{info['name']}_clean.pdb"
    n_residues = write_clean_pdb(residues, stretch, clean_pdb)
    print(f"  Wrote {n_residues} residues to {clean_pdb.name}")
    
    result['clean_pdb'] = str(clean_pdb)
    result['residues'] = n_residues
    result['original_range'] = f"{min(stretch)}-{max(stretch)}"
    
    # Prepare with OpenMM
    print(f"\n  OpenMM preparation:")
    prep_result = prepare_with_openmm(clean_pdb, receptor_dir, info['name'], solvate=True)
    result['preparation'] = prep_result
    
    if not prep_result['success']:
        return result
    
    # Run minimization
    print(f"\n  Running minimization and equilibration:")
    min_result = run_minimization(
        Path(prep_result['output_files']['solvated']),
        receptor_dir,
        min_steps=1000,
        eq_steps=5000
    )
    result['minimization'] = min_result
    
    result['success'] = min_result['success']
    
    if result['success']:
        print(f"\n  ✓ {info['name']} READY FOR SIMULATION!")
        print(f"    Final structure: {min_result.get('equilibrated_pdb', '')}")
        print(f"    Total atoms: {min_result.get('atoms', 0)}")
    
    return result


def main():
    """Main entry point."""
    base_dir = Path(__file__).parent
    structures_dir = base_dir / 'structures'
    output_dir = base_dir / 'prepared_receptors'
    
    structures_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print(" RECEPTOR STRUCTURE PREPARATION WITH TERMINAL CAPPING")
    print("=" * 70)
    
    results = {}
    successful = []
    
    for pdb_id, info in RECEPTORS.items():
        result = process_receptor(pdb_id, info, structures_dir, output_dir)
        results[pdb_id] = result
        
        if result['success']:
            successful.append(info['name'])
    
    # Final summary
    print("\n" + "=" * 70)
    print(" FINAL SUMMARY")
    print("=" * 70)
    
    for pdb_id, result in results.items():
        name = result['name']
        if result['success']:
            min_data = result.get('minimization', {})
            print(f"  ✓ {name}: {result['residues']} residues, "
                  f"{min_data.get('atoms', 0)} atoms, "
                  f"E={min_data.get('equilibrated_energy', 0):.0f} kJ/mol")
        else:
            print(f"  ✗ {name}: {result.get('error', result.get('preparation', {}).get('error', 'Unknown error'))[:50]}")
    
    print(f"\n  Successfully prepared: {len(successful)}/{len(RECEPTORS)} receptors")
    
    if successful:
        print(f"\n  Ready for MD simulations:")
        for name in successful:
            print(f"    - {output_dir / name / 'equilibrated.pdb'}")
    
    # Save results
    results_file = base_dir / 'receptor_preparation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n  Results saved to: {results_file}")
    
    return results


if __name__ == '__main__':
    main()
