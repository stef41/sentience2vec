#!/usr/bin/env python3
"""
Complete Receptor Structure Preparation for MD Simulations
===========================================================
Adds proper terminal groups (OXT for C-terminus) for OpenMM compatibility.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import math
import warnings
warnings.filterwarnings('ignore')

# Standard amino acids
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'
}

# Required heavy atoms
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
    """Check if residue has all required heavy atoms."""
    if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
        required = REQUIRED_ATOMS.get('HIS', set())
    else:
        required = REQUIRED_ATOMS.get(res_name, set())
    return required.issubset(atoms) if required else True


def parse_pdb_chain(pdb_path: Path, target_chain: str) -> dict:
    """Parse PDB and extract complete residues from target chain."""
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
            
            if chain != target_chain:
                continue
            if altloc not in [' ', 'A']:
                continue
            if res_name not in STANDARD_AA:
                continue
            
            if res_num not in residues:
                residues[res_num] = {'resname': res_name, 'atoms': {}}
            
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occupancy = float(line[54:60]) if len(line) > 60 else 1.0
            bfactor = float(line[60:66]) if len(line) > 66 else 0.0
            element = line[76:78].strip() if len(line) > 78 else atom_name[0]
            
            residues[res_num]['atoms'][atom_name] = {
                'x': x, 'y': y, 'z': z,
                'occupancy': occupancy,
                'bfactor': bfactor,
                'element': element
            }
    
    # Filter to complete residues
    complete = {}
    for res_num, data in residues.items():
        atom_names = set(data['atoms'].keys())
        if check_residue_complete(data['resname'], atom_names):
            complete[res_num] = data
    
    return complete


def find_longest_stretch(residues: dict, max_gap: int = 2) -> list:
    """Find longest continuous stretch allowing small gaps."""
    if not residues:
        return []
    
    sorted_nums = sorted(residues.keys())
    stretches = []
    current = [sorted_nums[0]]
    
    for i in range(1, len(sorted_nums)):
        gap = sorted_nums[i] - sorted_nums[i-1]
        if gap <= max_gap + 1:
            current.append(sorted_nums[i])
        else:
            stretches.append(current)
            current = [sorted_nums[i]]
    stretches.append(current)
    
    longest = max(stretches, key=len)
    return [n for n in longest if n in residues]


def calculate_oxt_position(c_coords: dict, o_coords: dict, ca_coords: dict) -> dict:
    """
    Calculate OXT position for proper carboxylate geometry.
    OXT should be roughly symmetric to O about the C-CA axis.
    """
    # Vectors
    c = [c_coords['x'], c_coords['y'], c_coords['z']]
    o = [o_coords['x'], o_coords['y'], o_coords['z']]
    ca = [ca_coords['x'], ca_coords['y'], ca_coords['z']]
    
    # C-CA vector
    ca_c = [c[i] - ca[i] for i in range(3)]
    
    # C-O vector
    c_o = [o[i] - c[i] for i in range(3)]
    
    # Normalize CA-C
    ca_c_len = math.sqrt(sum(x*x for x in ca_c))
    ca_c_norm = [x/ca_c_len for x in ca_c]
    
    # Project C-O onto CA-C plane, then reflect
    # Simple approach: rotate C-O by ~120 degrees around CA-C axis
    # Use cross product to get perpendicular direction
    
    # Get perpendicular in the plane
    dot = sum(c_o[i] * ca_c_norm[i] for i in range(3))
    perp = [c_o[i] - dot * ca_c_norm[i] for i in range(3)]
    perp_len = math.sqrt(sum(x*x for x in perp))
    if perp_len > 0.1:
        perp = [x/perp_len for x in perp]
    
    # Cross product for third axis
    cross = [
        ca_c_norm[1]*perp[2] - ca_c_norm[2]*perp[1],
        ca_c_norm[2]*perp[0] - ca_c_norm[0]*perp[2],
        ca_c_norm[0]*perp[1] - ca_c_norm[1]*perp[0]
    ]
    
    # OXT at ~120 degrees from O
    # cos(120) = -0.5, sin(120) = 0.866
    co_len = math.sqrt(sum(x*x for x in c_o))  # ~1.25 Å for C=O
    
    oxt = [
        c[i] + co_len * (-0.5 * perp[i] + 0.866 * cross[i])
        for i in range(3)
    ]
    
    return {
        'x': oxt[0], 'y': oxt[1], 'z': oxt[2],
        'occupancy': 1.0,
        'bfactor': o_coords['bfactor'],
        'element': 'O'
    }


def write_pdb_with_termini(residues: dict, res_nums: list, output_path: Path) -> int:
    """Write clean PDB with proper terminal groups."""
    lines = []
    atom_num = 1
    new_res_num = 1
    
    last_res_data = None
    last_res_num_new = None
    
    for orig_num in res_nums:
        if orig_num not in residues:
            continue
        
        data = residues[orig_num]
        res_name = data['resname']
        
        # Fix histidine naming for Amber
        if res_name == 'HIS':
            res_name = 'HID'
        
        for atom_name, coords in data['atoms'].items():
            line = f"ATOM  {atom_num:5d} {atom_name:^4s} {res_name:3s} A{new_res_num:4d}    "
            line += f"{coords['x']:8.3f}{coords['y']:8.3f}{coords['z']:8.3f}"
            line += f"{coords['occupancy']:6.2f}{coords['bfactor']:6.2f}"
            line += f"          {coords['element']:>2s}\n"
            
            lines.append(line)
            atom_num += 1
        
        last_res_data = data
        last_res_num_new = new_res_num
        new_res_num += 1
    
    # Add OXT to C-terminus if we have the necessary atoms
    if last_res_data and 'C' in last_res_data['atoms'] and 'O' in last_res_data['atoms']:
        c_coords = last_res_data['atoms']['C']
        o_coords = last_res_data['atoms']['O']
        ca_coords = last_res_data['atoms'].get('CA')
        
        if ca_coords:
            oxt_coords = calculate_oxt_position(c_coords, o_coords, ca_coords)
            res_name = last_res_data['resname']
            if res_name == 'HIS':
                res_name = 'HID'
            
            oxt_line = f"ATOM  {atom_num:5d}  OXT {res_name:3s} A{last_res_num_new:4d}    "
            oxt_line += f"{oxt_coords['x']:8.3f}{oxt_coords['y']:8.3f}{oxt_coords['z']:8.3f}"
            oxt_line += f"{oxt_coords['occupancy']:6.2f}{oxt_coords['bfactor']:6.2f}"
            oxt_line += f"           O\n"
            
            lines.append(oxt_line)
            atom_num += 1
    
    lines.append('TER\n')
    lines.append('END\n')
    
    with open(output_path, 'w') as f:
        f.writelines(lines)
    
    return len(res_nums)


def prepare_openmm(input_pdb: Path, output_dir: Path, name: str) -> dict:
    """Prepare structure with OpenMM."""
    result = {'success': False, 'output_files': {}}
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        pdb = app.PDBFile(str(input_pdb))
        print(f"    Loaded: {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues")
        
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        print(f"    Adding hydrogens...")
        modeller.addHydrogens(forcefield, pH=7.4)
        print(f"    After H: {modeller.topology.getNumAtoms()} atoms")
        
        # Save protein
        protein_pdb = output_dir / f"{name}_prepared.pdb"
        with open(protein_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        result['output_files']['protein'] = str(protein_pdb)
        
        # Add solvent
        print(f"    Adding solvent...")
        modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometers, neutralize=True)
        print(f"    Solvated: {modeller.topology.getNumAtoms()} atoms")
        
        solvated_pdb = output_dir / f"{name}_solvated.pdb"
        with open(solvated_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        result['output_files']['solvated'] = str(solvated_pdb)
        
        # Validate
        print(f"    Validating system...")
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        result['success'] = True
        result['atoms'] = modeller.topology.getNumAtoms()
        print(f"    ✓ System created!")
        
    except Exception as e:
        result['error'] = str(e)
        print(f"    ✗ Error: {str(e)[:80]}")
    
    return result


def run_minimization(solvated_pdb: Path, output_dir: Path) -> dict:
    """Run energy minimization and equilibration."""
    result = {'success': False}
    
    try:
        pdb = app.PDBFile(str(solvated_pdb))
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
        
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            platform = mm.Platform.getPlatformByName('CPU')
        
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        
        state = simulation.context.getState(getEnergy=True)
        init_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Initial E: {init_E:.1f} kJ/mol")
        
        print(f"      Minimizing...")
        simulation.minimizeEnergy(maxIterations=1000)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        min_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Minimized E: {min_E:.1f} kJ/mol")
        
        minimized_pdb = output_dir / 'minimized.pdb'
        with open(minimized_pdb, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        print(f"      Equilibrating (10 ps)...")
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
        simulation.step(5000)
        
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        eq_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"      Equilibrated E: {eq_E:.1f} kJ/mol")
        
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


def process_receptor(pdb_id: str, info: dict, structures_dir: Path, output_dir: Path) -> dict:
    """Process a single receptor."""
    result = {'pdb_id': pdb_id, 'name': info['name'], 'success': False}
    
    print(f"\n{'='*70}")
    print(f" {pdb_id}: {info['name']} - {info['desc']}")
    print(f"{'='*70}")
    
    raw_pdb = download_pdb(pdb_id, structures_dir)
    
    print(f"  Parsing chain {info['chain']}...")
    residues = parse_pdb_chain(raw_pdb, info['chain'])
    print(f"  Complete residues: {len(residues)}")
    
    if len(residues) < 30:
        result['error'] = f"Only {len(residues)} complete residues"
        return result
    
    stretch = find_longest_stretch(residues)
    print(f"  Longest stretch: {len(stretch)} residues")
    
    if len(stretch) < 30:
        result['error'] = f"Stretch too short ({len(stretch)} residues)"
        return result
    
    receptor_dir = output_dir / info['name']
    receptor_dir.mkdir(parents=True, exist_ok=True)
    
    clean_pdb = receptor_dir / f"{info['name']}_clean.pdb"
    n_res = write_pdb_with_termini(residues, stretch, clean_pdb)
    print(f"  Wrote {n_res} residues (with OXT)")
    result['residues'] = n_res
    result['clean_pdb'] = str(clean_pdb)
    
    print(f"\n  OpenMM preparation:")
    prep = prepare_openmm(clean_pdb, receptor_dir, info['name'])
    result['preparation'] = prep
    
    if not prep['success']:
        return result
    
    print(f"\n  Minimization and equilibration:")
    min_result = run_minimization(Path(prep['output_files']['solvated']), receptor_dir)
    result['minimization'] = min_result
    result['success'] = min_result['success']
    
    if result['success']:
        print(f"\n  ✓ {info['name']} READY!")
    
    return result


def main():
    """Main."""
    base_dir = Path(__file__).parent
    structures_dir = base_dir / 'structures'
    output_dir = base_dir / 'prepared_receptors'
    
    structures_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print(" RECEPTOR STRUCTURE PREPARATION WITH PROPER TERMINI")
    print("=" * 70)
    
    results = {}
    successful = []
    
    for pdb_id, info in RECEPTORS.items():
        result = process_receptor(pdb_id, info, structures_dir, output_dir)
        results[pdb_id] = result
        if result['success']:
            successful.append(info['name'])
    
    print("\n" + "=" * 70)
    print(" SUMMARY")
    print("=" * 70)
    
    for pdb_id, result in results.items():
        name = result['name']
        if result['success']:
            min_data = result.get('minimization', {})
            print(f"  ✓ {name}: {result['residues']} res, {min_data.get('atoms', 0)} atoms")
        else:
            err = result.get('error', result.get('preparation', {}).get('error', 'Unknown'))
            print(f"  ✗ {name}: {err[:60]}")
    
    print(f"\n  Success: {len(successful)}/{len(RECEPTORS)}")
    
    with open(base_dir / 'preparation_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    return results


if __name__ == '__main__':
    main()
