#!/usr/bin/env python3
"""
Energy Minimization for KOR_salvinorin_A
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import sys

print("="*60)
print("ENERGY MINIMIZATION: KOR_salvinorin_A")
print("="*60)

# Load receptor
print("\nLoading receptor...")
receptor = app.PDBFile("/data/users/zacharie/whc/consciousness_study/md_simulations/structures/4DJH.pdb")
print(f"  Atoms: {receptor.topology.getNumAtoms()}")

# Set up force field
print("\nSetting up force field...")
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create modeller
modeller = app.Modeller(receptor.topology, receptor.positions)

# Add hydrogens if missing
print("Adding hydrogens...")
modeller.addHydrogens(forcefield)

# Add solvent
print("Adding solvent (this may take a minute)...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0 * unit.nanometers,
    ionicStrength=0.15 * unit.molar
)
print(f"  Total atoms after solvation: {modeller.topology.getNumAtoms()}")

# Create system
print("\nCreating system...")
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

# Get best platform
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    print("Using CUDA platform")
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {}
        print("Using OpenCL platform")
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {}
        print("Using CPU platform")

# Create simulation
simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties if properties else {}
)
simulation.context.setPositions(modeller.positions)

# Get initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"\nInitial energy: {initial_energy}")

# Minimize
print("\nMinimizing energy (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

# Get final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"Final energy: {final_energy}")
print(f"Energy reduced by: {initial_energy - final_energy}")

# Save minimized structure
output_pdb = "/data/users/zacharie/whc/consciousness_study/md_simulations/KOR_salvinorin_A/minimized.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(
        simulation.topology,
        state.getPositions(),
        f
    )
print(f"\nSaved minimized structure to: {output_pdb}")
print("\n✓ Minimization complete!")
