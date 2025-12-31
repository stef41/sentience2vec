#!/usr/bin/env python3
"""
Comprehensive Receptor and Ligand Preparation Pipeline

Downloads all 20 neuroreceptor structures from PDB and prepares all 91 ligand
3D structures for molecular dynamics simulations.
"""

import os
import sys
import json
import math
import time
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import our databases
from receptor_pdb_database import RECEPTOR_PDB_DATABASE
from ligand_smiles_database import LIGAND_DATABASE

# Try to import optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("Warning: RDKit not available. Ligand preparation will be skipped.")

try:
    import openmm.app as app
    import openmm as mm
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    print("Warning: OpenMM not available. Receptor preparation will be skipped.")

# Constants
BASE_DIR = Path(__file__).parent
STRUCTURES_DIR = BASE_DIR / "structures"
PREPARED_RECEPTORS_DIR = BASE_DIR / "prepared_receptors"
LIGANDS_DIR = BASE_DIR / "ligands"

# Standard amino acids
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX'
}

# Required atoms for each amino acid (for completeness checking)
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


# =============================================================================
# RECEPTOR PREPARATION FUNCTIONS
# =============================================================================

def download_pdb(pdb_id, output_path):
    """Download a PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        print(f"  Failed to download {pdb_id}: {e}")
        return False


def check_residue_complete(res_name, atoms):
    """Check if a residue has all required atoms."""
    if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
        required = REQUIRED_ATOMS.get('HIS', set())
    else:
        required = REQUIRED_ATOMS.get(res_name, set())
    return required.issubset(atoms) if required else True


def parse_pdb_chain(pdb_path, chain_id):
    """Parse a PDB file and extract complete residues from specified chain."""
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
            
            # Filter by chain and standard amino acids
            if chain != chain_id or altloc not in [' ', 'A'] or res_name not in STANDARD_AA:
                continue
            
            if res_num not in residues:
                residues[res_num] = {'resname': res_name, 'atoms': {}}
            
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ = float(line[54:60]) if len(line) > 60 else 1.0
            bf = float(line[60:66]) if len(line) > 66 else 0.0
            elem = line[76:78].strip() if len(line) > 78 else atom_name[0]
            
            residues[res_num]['atoms'][atom_name] = {
                'x': x, 'y': y, 'z': z, 'occ': occ, 'bf': bf, 'elem': elem
            }
    
    # Keep only complete residues
    return {
        k: v for k, v in residues.items()
        if check_residue_complete(v['resname'], set(v['atoms'].keys()))
    }


def find_longest_stretch(residues, max_gap=2):
    """Find the longest continuous stretch of residues."""
    if not residues:
        return []
    
    nums = sorted(residues.keys())
    stretches = []
    current = [nums[0]]
    
    for i in range(1, len(nums)):
        if nums[i] - nums[i-1] <= max_gap + 1:
            current.append(nums[i])
        else:
            stretches.append(current)
            current = [nums[i]]
    stretches.append(current)
    
    longest = max(stretches, key=len)
    return [n for n in longest if n in residues]


def calculate_oxt_position(c_atom, o_atom, ca_atom):
    """Calculate OXT position for C-terminal."""
    # CA -> C vector
    ca_c = [c_atom['x'] - ca_atom['x'],
            c_atom['y'] - ca_atom['y'],
            c_atom['z'] - ca_atom['z']]
    ca_c_len = math.sqrt(sum(x*x for x in ca_c))
    ca_c_norm = [x / ca_c_len for x in ca_c]
    
    # C -> O vector
    c_o = [o_atom['x'] - c_atom['x'],
           o_atom['y'] - c_atom['y'],
           o_atom['z'] - c_atom['z']]
    
    # Project out CA-C component
    dot = sum(c_o[i] * ca_c_norm[i] for i in range(3))
    perp = [c_o[i] - dot * ca_c_norm[i] for i in range(3)]
    perp_len = math.sqrt(sum(x*x for x in perp))
    
    if perp_len > 0.1:
        perp = [x / perp_len for x in perp]
    
    # Cross product for 120 degree rotation
    cross = [
        ca_c_norm[1] * perp[2] - ca_c_norm[2] * perp[1],
        ca_c_norm[2] * perp[0] - ca_c_norm[0] * perp[2],
        ca_c_norm[0] * perp[1] - ca_c_norm[1] * perp[0]
    ]
    
    co_len = math.sqrt(sum(x*x for x in c_o))
    
    return {
        'x': c_atom['x'] + co_len * (-0.5 * perp[0] + 0.866 * cross[0]),
        'y': c_atom['y'] + co_len * (-0.5 * perp[1] + 0.866 * cross[1]),
        'z': c_atom['z'] + co_len * (-0.5 * perp[2] + 0.866 * cross[2]),
        'occ': 1.0,
        'bf': o_atom['bf'],
        'elem': 'O'
    }


def write_clean_pdb(residues, residue_nums, output_path):
    """Write a clean PDB file with proper termini."""
    lines = []
    atom_num = 1
    res_num_out = 1
    
    last_residue_data = None
    last_res_num_out = None
    
    for orig_num in residue_nums:
        if orig_num not in residues:
            continue
        
        data = residues[orig_num]
        res_name = data['resname']
        
        # Convert HIS to HID for Amber
        if res_name == 'HIS':
            res_name = 'HID'
        
        for atom_name, coords in data['atoms'].items():
            line = f"ATOM  {atom_num:5d} {atom_name:^4s} {res_name:3s} A{res_num_out:4d}    "
            line += f"{coords['x']:8.3f}{coords['y']:8.3f}{coords['z']:8.3f}"
            line += f"{coords['occ']:6.2f}{coords['bf']:6.2f}          {coords['elem']:>2s}\n"
            lines.append(line)
            atom_num += 1
        
        last_residue_data = data
        last_res_num_out = res_num_out
        res_num_out += 1
    
    # Add OXT to C-terminus
    if last_residue_data:
        atoms = last_residue_data['atoms']
        if 'C' in atoms and 'O' in atoms and 'CA' in atoms:
            oxt = calculate_oxt_position(atoms['C'], atoms['O'], atoms['CA'])
            res_name = last_residue_data['resname']
            if res_name == 'HIS':
                res_name = 'HID'
            
            line = f"ATOM  {atom_num:5d}  OXT {res_name:3s} A{last_res_num_out:4d}    "
            line += f"{oxt['x']:8.3f}{oxt['y']:8.3f}{oxt['z']:8.3f}"
            line += f"{oxt['occ']:6.2f}{oxt['bf']:6.2f}           O\n"
            lines.append(line)
    
    lines.append("TER\n")
    lines.append("END\n")
    
    with open(output_path, 'w') as f:
        f.writelines(lines)
    
    return len(residue_nums)


def prepare_receptor_openmm(clean_pdb_path, output_dir):
    """Prepare receptor with OpenMM: add H, solvate, minimize, equilibrate."""
    if not HAS_OPENMM:
        return False
    
    try:
        # Load structure
        pdb = app.PDBFile(str(clean_pdb_path))
        
        # Setup force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        
        # Create modeller and add hydrogens
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)
        
        # Save with hydrogens
        with open(output_dir / "with_hydrogens.pdb", 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
        # Add solvent
        modeller.addSolvent(forcefield, model='tip3p', padding=0.8*unit.nanometers)
        
        # Save solvated
        with open(output_dir / "solvated.pdb", 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
        # Create system
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        # Setup integrator and simulation
        integrator = mm.LangevinMiddleIntegrator(
            310*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds
        )
        
        # Use OpenCL if available
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            platform = mm.Platform.getPlatformByName('CPU')
        
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=1000)
        state = simulation.context.getState(getPositions=True)
        
        with open(output_dir / "minimized.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Short equilibration
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
        simulation.step(5000)  # 10 ps
        
        state = simulation.context.getState(getPositions=True)
        with open(output_dir / "equilibrated.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        return True
        
    except Exception as e:
        print(f"  OpenMM preparation failed: {e}")
        return False


def prepare_all_receptors():
    """Download and prepare all receptors."""
    print("=" * 70)
    print("RECEPTOR PREPARATION")
    print("=" * 70)
    
    STRUCTURES_DIR.mkdir(exist_ok=True)
    PREPARED_RECEPTORS_DIR.mkdir(exist_ok=True)
    
    results = {}
    
    for receptor_name, info in RECEPTOR_PDB_DATABASE.items():
        pdb_id = info['pdb_id']
        chain = info['chain']
        
        print(f"\n[{receptor_name}]")
        print(f"  PDB: {pdb_id}, Chain: {chain}")
        
        # Create output directory
        output_dir = PREPARED_RECEPTORS_DIR / receptor_name
        output_dir.mkdir(exist_ok=True)
        
        # Check if already prepared
        equilibrated = output_dir / "equilibrated.pdb"
        if equilibrated.exists():
            print(f"  ✓ Already prepared")
            results[receptor_name] = 'cached'
            continue
        
        # Download PDB
        pdb_path = STRUCTURES_DIR / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            print(f"  Downloading {pdb_id}...", end=' ')
            if download_pdb(pdb_id, pdb_path):
                print("✓")
            else:
                results[receptor_name] = 'download_failed'
                continue
        
        # Parse and extract clean chain
        print(f"  Extracting chain {chain}...", end=' ')
        residues = parse_pdb_chain(pdb_path, chain)
        if not residues:
            print(f"✗ No residues found")
            results[receptor_name] = 'no_residues'
            continue
        
        stretch = find_longest_stretch(residues)
        if len(stretch) < 30:
            print(f"✗ Only {len(stretch)} residues")
            results[receptor_name] = 'too_short'
            continue
        
        clean_pdb = output_dir / f"{receptor_name}_clean.pdb"
        n_res = write_clean_pdb(residues, stretch, clean_pdb)
        print(f"✓ {n_res} residues")
        
        # OpenMM preparation
        if HAS_OPENMM:
            print(f"  OpenMM preparation...", end=' ')
            if prepare_receptor_openmm(clean_pdb, output_dir):
                print("✓")
                results[receptor_name] = 'success'
            else:
                results[receptor_name] = 'openmm_failed'
        else:
            results[receptor_name] = 'no_openmm'
    
    return results


# =============================================================================
# LIGAND PREPARATION FUNCTIONS
# =============================================================================

def prepare_ligand_3d(name, smiles, output_dir):
    """Generate 3D structure from SMILES using RDKit."""
    if not HAS_RDKIT:
        return False
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # Try with different parameters
            result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
            if result != 0:
                return False
        
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        
        # Create output directory
        ligand_dir = output_dir / name
        ligand_dir.mkdir(parents=True, exist_ok=True)
        
        # Save as PDB
        Chem.MolToPDBFile(mol, str(ligand_dir / f"{name}.pdb"))
        
        # Save as SDF (for better atom typing)
        writer = Chem.SDWriter(str(ligand_dir / f"{name}.sdf"))
        writer.write(mol)
        writer.close()
        
        # Save as MOL2 format string (for docking)
        # Note: RDKit doesn't directly write MOL2, but SDF is usually sufficient
        
        return True
        
    except Exception as e:
        print(f"  Failed to prepare {name}: {e}")
        return False


def prepare_all_ligands():
    """Prepare 3D structures for all ligands."""
    print("\n" + "=" * 70)
    print("LIGAND PREPARATION")
    print("=" * 70)
    
    if not HAS_RDKIT:
        print("RDKit not available - skipping ligand preparation")
        return {}
    
    LIGANDS_DIR.mkdir(exist_ok=True)
    
    results = {}
    total = len(LIGAND_DATABASE)
    
    for i, (name, data) in enumerate(LIGAND_DATABASE.items(), 1):
        smiles = data['smiles']
        
        # Check if already prepared
        ligand_dir = LIGANDS_DIR / name
        pdb_file = ligand_dir / f"{name}.pdb"
        
        if pdb_file.exists():
            results[name] = 'cached'
            continue
        
        print(f"  [{i}/{total}] {name}...", end=' ')
        
        if prepare_ligand_3d(name, smiles, LIGANDS_DIR):
            print("✓")
            results[name] = 'success'
        else:
            print("✗")
            results[name] = 'failed'
    
    return results


# =============================================================================
# MAIN
# =============================================================================

def print_summary(receptor_results, ligand_results):
    """Print preparation summary."""
    print("\n" + "=" * 70)
    print("PREPARATION SUMMARY")
    print("=" * 70)
    
    # Receptor summary
    receptor_success = sum(1 for v in receptor_results.values() if v in ['success', 'cached'])
    print(f"\n📦 Receptors: {receptor_success}/{len(receptor_results)} prepared")
    
    for name, status in sorted(receptor_results.items()):
        symbol = "✓" if status in ['success', 'cached'] else "✗"
        print(f"  {symbol} {name}: {status}")
    
    # Ligand summary
    ligand_success = sum(1 for v in ligand_results.values() if v in ['success', 'cached'])
    print(f"\n💊 Ligands: {ligand_success}/{len(ligand_results)} prepared")
    
    failed = [name for name, status in ligand_results.items() if status == 'failed']
    if failed:
        print(f"  Failed: {', '.join(failed[:10])}" + ("..." if len(failed) > 10 else ""))
    
    # Save results
    results = {
        'receptors': receptor_results,
        'ligands': ligand_results,
        'summary': {
            'total_receptors': len(receptor_results),
            'receptors_prepared': receptor_success,
            'total_ligands': len(ligand_results),
            'ligands_prepared': ligand_success
        }
    }
    
    with open(BASE_DIR / 'preparation_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to preparation_results.json")


def main():
    """Main entry point."""
    print("=" * 70)
    print("COMPREHENSIVE MD SIMULATION PREPARATION")
    print("All Neuroreceptors + All Psychoactive Ligands")
    print("=" * 70)
    print(f"\nTarget: {len(RECEPTOR_PDB_DATABASE)} receptors, {len(LIGAND_DATABASE)} ligands")
    
    # Prepare receptors
    receptor_results = prepare_all_receptors()
    
    # Prepare ligands
    ligand_results = prepare_all_ligands()
    
    # Print summary
    print_summary(receptor_results, ligand_results)


if __name__ == '__main__':
    main()
