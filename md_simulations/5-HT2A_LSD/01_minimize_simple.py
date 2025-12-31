#!/usr/bin/env python3
"""
Simplified Energy Minimization for MD Simulations
Uses OpenMM's built-in PDB fixer capabilities
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import sys

# Configuration
SYSTEM_DIR = Path(__file__).parent
SYSTEM_NAME = SYSTEM_DIR.name

# Map receptor to PDB
RECEPTOR_PDB = {
    "5-HT2A": "6WHA", "5-HT2B": "5TVN", "D2": "6CM4",
    "MOR": "5C1M", "KOR": "4DJH", "NMDA": "4PE5", "CB1": "5TGZ"
}

# Get receptor name
parts = SYSTEM_NAME.split('_')
RECEPTOR_NAME = parts[0]
LIGAND_NAME = '_'.join(parts[1:])

print("="*60)
print(f"ENERGY MINIMIZATION: {SYSTEM_NAME}")
print("="*60)

# Paths
structures_dir = SYSTEM_DIR.parent / "structures"
pdb_id = RECEPTOR_PDB.get(RECEPTOR_NAME)

# Find the cleaned PDB
pdb_path = structures_dir / f"{pdb_id}_fixed.pdb"
if not pdb_path.exists():
    pdb_path = structures_dir / f"{pdb_id}.pdb"

print(f"\nReceptor: {RECEPTOR_NAME} (PDB: {pdb_id})")
print(f"Ligand: {LIGAND_NAME}")
print(f"PDB file: {pdb_path}")

# Load PDB
print("\nLoading structure...")
pdb = app.PDBFile(str(pdb_path))
print(f"  Initial atoms: {pdb.topology.getNumAtoms()}")

# Create modeller
modeller = app.Modeller(pdb.topology, pdb.positions)

# Remove non-protein atoms first to simplify 
print("Cleaning structure...")
chains_to_keep = set()
for chain in modeller.topology.chains():
    # Keep chains that have protein residues
    for res in chain.residues():
        if res.name in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                        'HIS', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
            chains_to_keep.add(chain.id)
            break

# Only keep the first chain to simplify
first_chain = list(modeller.topology.chains())[0]
print(f"  Using chain: {first_chain.id}")

# Delete all atoms not in first chain
to_delete = []
for chain in modeller.topology.chains():
    if chain.id != first_chain.id:
        to_delete.extend(list(chain.residues()))

# Also delete any non-standard residues
for residue in first_chain.residues():
    if residue.name not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                            'HIS', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 
                            'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
        to_delete.append(residue)

modeller.delete(to_delete)
print(f"  Atoms after cleaning: {modeller.topology.getNumAtoms()}")

# Set up force field
print("\nSetting up simulation...")
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Add missing hydrogens
print("Adding missing atoms and hydrogens...")
try:
    modeller.addHydrogens(forcefield)
    print(f"  Atoms after hydrogens: {modeller.topology.getNumAtoms()}")
except Exception as e:
    print(f"  Warning: Could not add hydrogens: {e}")
    print("  Proceeding without adding hydrogens...")

# Add solvent (smaller box for faster test)
print("Adding solvent...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=0.8 * unit.nanometers,  # Smaller padding for faster test
    ionicStrength=0.15 * unit.molar
)
print(f"  Total atoms: {modeller.topology.getNumAtoms()}")

# Create system
print("\nCreating OpenMM system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds
)

# Integrator
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picoseconds,
    0.002 * unit.picoseconds
)

# Platform
platform_name = "CPU"
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    platform_name = "CUDA"
    properties = {'CudaPrecision': 'mixed'}
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        platform_name = "OpenCL"
        properties = {}
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {}
print(f"Using platform: {platform_name}")

# Create simulation
simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties if properties else {}
)
simulation.context.setPositions(modeller.positions)

# Initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"\nInitial energy: {initial_energy}")

# Minimize
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=1000)

# Final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"Final energy: {final_energy}")
print(f"Energy reduced by: {initial_energy - final_energy}")

# Save
output_pdb = SYSTEM_DIR / "minimized.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"\n✓ Saved minimized structure: {output_pdb}")
