#!/usr/bin/env python3
"""
Energy Minimization for 5-HT2A_DMT
"""
import warnings
warnings.filterwarnings('ignore')

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer
from pathlib import Path
import sys
import time

def progress(msg):
    """Print timestamped progress message."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

print("="*60)
print("ENERGY MINIMIZATION: 5-HT2A_DMT")
print("="*60)

# Load and fix receptor with PDBFixer
progress("Loading receptor...")
fixer = PDBFixer("/data/users/zacharie/whc/consciousness_study/md_simulations/structures/6WHA.pdb")
progress(f"  Loaded: {fixer.topology.getNumAtoms()} atoms")

progress("Fixing missing atoms...")
fixer.findMissingResidues()
fixer.missingResidues = {}  # Skip missing residues for speed
fixer.findMissingAtoms()
fixer.addMissingAtoms()

progress("Adding hydrogens...")
fixer.addMissingHydrogens(7.0)

progress("Removing heterogens...")
fixer.removeHeterogens(keepWater=False)
progress(f"  Protein atoms: {fixer.topology.getNumAtoms()}")

# Set up force field
progress("Setting up force field...")
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create modeller from fixed structure
modeller = app.Modeller(fixer.topology, fixer.positions)

# Add solvent with progress
progress("Adding solvent (0.8nm padding)...")
start = time.time()
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=0.8 * unit.nanometers,
    ionicStrength=0.15 * unit.molar
)
progress(f"  Solvated: {modeller.topology.getNumAtoms()} atoms ({time.time()-start:.1f}s)")

# Create system
progress("Creating system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds
)

# Set up integrator (doesn't matter much for minimization)
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picoseconds,
    0.002 * unit.picoseconds
)

# Get best platform - use CPU since GPUs may be busy
progress("Using CPU platform (32 threads)...")
platform = mm.Platform.getPlatformByName('CPU')
properties = {'Threads': '32'}

# Create simulation
progress("Creating simulation...")
simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties
)
simulation.context.setPositions(modeller.positions)

# Get initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
progress(f"Initial energy: {initial_energy}")

# Minimize with progress
progress("Minimizing energy (5000 steps)...")
start = time.time()
simulation.minimizeEnergy(maxIterations=5000)

# Get final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
progress(f"Final energy: {final_energy} ({time.time()-start:.1f}s)")
progress(f"Energy reduced by: {initial_energy - final_energy}")

# Save minimized structure
output_pdb = "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_DMT/minimized.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(
        simulation.topology,
        state.getPositions(),
        f
    )
progress(f"Saved: {output_pdb}")
print("\n✓ Minimization complete!")
