#!/usr/bin/env python3
"""
Equilibration (NVT + NPT) for 5-HT2A_mescaline
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path

print("="*60)
print("EQUILIBRATION: 5-HT2A_mescaline")
print("="*60)

# Load minimized structure
print("\nLoading minimized structure...")
pdb = app.PDBFile("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/minimized.pdb")

# Force field
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system with barostat for NPT
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds
)

# Add barostat for pressure control
barostat = mm.MonteCarloBarostat(
    1.0 * unit.atmospheres,
    310 * unit.kelvin
)
system.addForce(barostat)

# Integrator
integrator = mm.LangevinMiddleIntegrator(
    310 * unit.kelvin,
    1.0 / unit.picoseconds,
    0.002 * unit.picoseconds
)

# Platform
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {}
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {}

# Create simulation
simulation = app.Simulation(
    pdb.topology,
    system,
    integrator,
    platform,
    properties if properties else {}
)
simulation.context.setPositions(pdb.positions)

# Initialize velocities at target temperature
print("\nInitializing velocities at 310 K...")
simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)

# Add reporters
simulation.reporters.append(
    app.StateDataReporter(
        '/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/equilibration.log',
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        density=True,
        speed=True
    )
)

# NVT equilibration (100 ps)
print("\nRunning NVT equilibration (100 ps)...")
nvt_steps = 50000  # 100 ps at 2 fs timestep
simulation.step(nvt_steps)

# NPT equilibration (400 ps)
print("Running NPT equilibration (400 ps)...")
npt_steps = 200000  # 400 ps
simulation.step(npt_steps)

# Save equilibrated state
state = simulation.context.getState(getPositions=True, getVelocities=True)
with open("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/equilibrated.pdb", 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

# Save checkpoint
simulation.saveCheckpoint("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/equilibrated.chk")

print("\n✓ Equilibration complete!")
print(f"  Saved equilibrated structure: /data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/equilibrated.pdb")
print(f"  Saved checkpoint: /data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline/equilibrated.chk")
