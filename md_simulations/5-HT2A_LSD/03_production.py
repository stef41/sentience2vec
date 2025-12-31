#!/usr/bin/env python3
"""
Production MD Simulation for 5-HT2A_LSD
Duration: 10.0 ns
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import time

print("="*60)
print("PRODUCTION MD: 5-HT2A_LSD")
print("Duration: 10.0 ns")
print("="*60)

# Load equilibrated structure
print("\nLoading equilibrated structure...")
pdb = app.PDBFile("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/equilibrated.pdb")

# Force field
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    hydrogenMass=4.0 * unit.amu  # HMR for longer timestep
)

# Barostat
barostat = mm.MonteCarloBarostat(
    1.0 * unit.atmospheres,
    310 * unit.kelvin
)
system.addForce(barostat)

# Integrator with 4 fs timestep (HMR enabled)
integrator = mm.LangevinMiddleIntegrator(
    310 * unit.kelvin,
    1.0 / unit.picoseconds,
    0.004 * unit.picoseconds  # 4 fs with HMR
)

# Platform
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
        print("Using CPU platform (will be slow)")

# Create simulation
simulation = app.Simulation(
    pdb.topology,
    system,
    integrator,
    platform,
    properties if properties else {}
)
simulation.context.setPositions(pdb.positions)

# Load checkpoint if available
checkpoint_file = "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/equilibrated.chk"
try:
    simulation.loadCheckpoint(checkpoint_file)
    print("Loaded checkpoint from equilibration")
except:
    simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)
    print("Starting with fresh velocities")

# Production reporters
simulation.reporters.append(
    app.DCDReporter('/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/trajectory.dcd', 5000)
)
simulation.reporters.append(
    app.StateDataReporter(
        '/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/production.log',
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        speed=True
    )
)
simulation.reporters.append(
    app.CheckpointReporter('/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/production.chk', 50000)
)

# Run production
production_steps = 5000000
print(f"\nRunning {production_steps} steps (10.0 ns)...")
start_time = time.time()

# Run in chunks to report progress
chunk_size = 25000  # ~50 ps per chunk
n_chunks = production_steps // chunk_size

for i in range(n_chunks):
    simulation.step(chunk_size)
    elapsed = time.time() - start_time
    progress = (i + 1) / n_chunks * 100
    time_ns = (i + 1) * chunk_size * 0.004 / 1000  # Convert to ns
    print(f"  Progress: {progress:.1f}% ({time_ns:.2f} ns) - Elapsed: {elapsed/60:.1f} min")

# Final steps
remaining = production_steps % chunk_size
if remaining > 0:
    simulation.step(remaining)

# Save final state
simulation.saveState("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/final_state.xml")
state = simulation.context.getState(getPositions=True)
with open("/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/final.pdb", 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

total_time = time.time() - start_time
print(f"\n✓ Production complete!")
print(f"  Total wall time: {total_time/3600:.2f} hours")
print(f"  Output: /data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD/trajectory.dcd")
