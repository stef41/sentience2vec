#!/usr/bin/env python3
"""
Launch All Molecular Dynamics Simulations
==========================================
Sets up and runs MD simulations for all receptor-ligand combinations
relevant to consciousness research.

Receptors: 5-HT2A, 5-HT2B, D2, MOR, KOR, NMDA, CB1
Ligands: LSD, psilocin, DMT, mescaline, MDMA, ketamine, salvinorin_A, THC
"""

import os
import sys
import json
import urllib.request
import time
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# ============================================================
# CONFIGURATION
# ============================================================

WORK_DIR = Path(__file__).parent.parent / "md_simulations"
WORK_DIR.mkdir(parents=True, exist_ok=True)

# Receptor-ligand pairs relevant for consciousness research
SIMULATION_PAIRS = [
    # Psychedelic binding
    ("5-HT2A", "LSD"),
    ("5-HT2A", "psilocin"),
    ("5-HT2A", "DMT"),
    ("5-HT2A", "mescaline"),
    
    # Entactogen binding
    ("5-HT2A", "MDMA"),
    ("5-HT2B", "MDMA"),
    
    # Dissociative binding
    ("NMDA", "ketamine"),
    
    # Kappa opioid (salvia)
    ("KOR", "salvinorin_A"),
    
    # Cannabinoid
    ("CB1", "THC"),
    
    # Dopaminergic (stimulant comparison)
    ("D2", "LSD"),  # LSD has D2 affinity
]

# Receptor PDB structures
RECEPTOR_STRUCTURES = {
    "5-HT2A": {
        "name": "Serotonin 2A Receptor",
        "pdb_id": "6WHA",  # LSD-bound structure
        "binding_site": [151, 155, 159, 207, 229, 230, 234, 313, 336, 340, 366]
    },
    "5-HT2B": {
        "name": "Serotonin 2B Receptor",
        "pdb_id": "5TVN",  # Ergotamine-bound
        "binding_site": [134, 135, 136, 139, 207, 217, 313, 317, 340, 344]
    },
    "D2": {
        "name": "Dopamine D2 Receptor",
        "pdb_id": "6CM4",
        "binding_site": [114, 115, 118, 119, 194, 197, 379, 380, 383, 386, 389]
    },
    "MOR": {
        "name": "Mu Opioid Receptor",
        "pdb_id": "5C1M",
        "binding_site": [148, 149, 202, 230, 233, 289, 293, 322, 325, 328]
    },
    "KOR": {
        "name": "Kappa Opioid Receptor", 
        "pdb_id": "4DJH",
        "binding_site": [138, 139, 141, 192, 219, 223, 280, 283, 284, 312]
    },
    "NMDA": {
        "name": "NMDA Receptor (GluN1-GluN2B)",
        "pdb_id": "4PE5",
        "binding_site": [518, 531, 532, 786, 788, 790, 810, 811, 813]
    },
    "CB1": {
        "name": "Cannabinoid CB1 Receptor",
        "pdb_id": "5TGZ",
        "binding_site": [178, 182, 192, 193, 196, 274, 277, 278, 353, 356, 360]
    }
}

# Ligand SMILES
LIGAND_SMILES = {
    "LSD": "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34",
    "psilocin": "CN(C)CCc1c[nH]c2ccc(O)cc12",
    "DMT": "CN(C)CCc1c[nH]c2ccccc12",
    "mescaline": "COc1cc(CCN)cc(OC)c1OC",
    "MDMA": "CC(NC)Cc1ccc2OCOc2c1",
    "ketamine": "CNC1(c2ccccc2Cl)CCCCC1=O",
    "salvinorin_A": "CC(=O)O[C@@H]1C[C@@H]2[C@]3(C)CC[C@@H]4C(C)(C)[C@H](OC(=O)c5ccoc5)C(=O)O[C@@]4(C)[C@@H]3C(=O)[C@@H]([C@@H]2C)[C@H]1OC(C)=O",
    "THC": "CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=C[C@H]21"
}


# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def download_pdb(pdb_id: str, output_path: Path) -> bool:
    """Download PDB structure from RCSB"""
    if output_path.exists():
        print(f"  ✓ PDB {pdb_id} already exists")
        return True
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        print(f"  Downloading {pdb_id}...", end=" ", flush=True)
        urllib.request.urlretrieve(url, output_path)
        print("✓")
        return True
    except Exception as e:
        print(f"✗ ({e})")
        return False


def prepare_ligand_3d(ligand_name: str, smiles: str, output_dir: Path) -> Optional[Path]:
    """Generate 3D structure from SMILES using RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        output_dir.mkdir(parents=True, exist_ok=True)
        mol_path = output_dir / f"{ligand_name}.mol"
        pdb_path = output_dir / f"{ligand_name}.pdb"
        
        if pdb_path.exists():
            print(f"  ✓ Ligand {ligand_name} already prepared")
            return pdb_path
        
        print(f"  Preparing {ligand_name}...", end=" ", flush=True)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("✗ (invalid SMILES)")
            return None
        
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result == -1:
            # Try with different parameters
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        if result == -1:
            print("✗ (embedding failed)")
            return None
        
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Save files
        Chem.MolToMolFile(mol, str(mol_path))
        Chem.MolToPDBFile(mol, str(pdb_path))
        
        print("✓")
        return pdb_path
        
    except ImportError:
        print("✗ (RDKit not available)")
        return None
    except Exception as e:
        print(f"✗ ({e})")
        return None


def check_openmm_available() -> Tuple[bool, str]:
    """Check if OpenMM is available and what platform"""
    try:
        import openmm as mm
        platforms = []
        for i in range(mm.Platform.getNumPlatforms()):
            plat = mm.Platform.getPlatform(i)
            platforms.append(plat.getName())
        
        # Prefer CUDA > OpenCL > CPU
        if "CUDA" in platforms:
            return True, "CUDA"
        elif "OpenCL" in platforms:
            return True, "OpenCL"
        elif "CPU" in platforms:
            return True, "CPU"
        else:
            return True, "Reference"
    except ImportError:
        return False, "Not installed"


# ============================================================
# SIMULATION SETUP
# ============================================================

def generate_minimization_script(system_dir: Path, receptor_pdb: str, ligand_pdb: str) -> Path:
    """Generate energy minimization script"""
    script = f'''#!/usr/bin/env python3
"""
Energy Minimization for {system_dir.name}
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import sys

print("="*60)
print("ENERGY MINIMIZATION: {system_dir.name}")
print("="*60)

# Load receptor
print("\\nLoading receptor...")
receptor = app.PDBFile("{receptor_pdb}")
print(f"  Atoms: {{receptor.topology.getNumAtoms()}}")

# Set up force field
print("\\nSetting up force field...")
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
print(f"  Total atoms after solvation: {{modeller.topology.getNumAtoms()}}")

# Create system
print("\\nCreating system...")
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
    properties = {{'CudaPrecision': 'mixed'}}
    print("Using CUDA platform")
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {{}}
        print("Using OpenCL platform")
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {{}}
        print("Using CPU platform")

# Create simulation
simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties if properties else {{}}
)
simulation.context.setPositions(modeller.positions)

# Get initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"\\nInitial energy: {{initial_energy}}")

# Minimize
print("\\nMinimizing energy (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

# Get final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"Final energy: {{final_energy}}")
print(f"Energy reduced by: {{initial_energy - final_energy}}")

# Save minimized structure
output_pdb = "{system_dir}/minimized.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(
        simulation.topology,
        state.getPositions(),
        f
    )
print(f"\\nSaved minimized structure to: {{output_pdb}}")
print("\\n✓ Minimization complete!")
'''
    
    script_path = system_dir / "01_minimize.py"
    script_path.write_text(script)
    return script_path


def generate_equilibration_script(system_dir: Path) -> Path:
    """Generate NVT/NPT equilibration script"""
    script = f'''#!/usr/bin/env python3
"""
Equilibration (NVT + NPT) for {system_dir.name}
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path

print("="*60)
print("EQUILIBRATION: {system_dir.name}")
print("="*60)

# Load minimized structure
print("\\nLoading minimized structure...")
pdb = app.PDBFile("{system_dir}/minimized.pdb")

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
    properties = {{'CudaPrecision': 'mixed'}}
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {{}}
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {{}}

# Create simulation
simulation = app.Simulation(
    pdb.topology,
    system,
    integrator,
    platform,
    properties if properties else {{}}
)
simulation.context.setPositions(pdb.positions)

# Initialize velocities at target temperature
print("\\nInitializing velocities at 310 K...")
simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)

# Add reporters
simulation.reporters.append(
    app.StateDataReporter(
        '{system_dir}/equilibration.log',
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
print("\\nRunning NVT equilibration (100 ps)...")
nvt_steps = 50000  # 100 ps at 2 fs timestep
simulation.step(nvt_steps)

# NPT equilibration (400 ps)
print("Running NPT equilibration (400 ps)...")
npt_steps = 200000  # 400 ps
simulation.step(npt_steps)

# Save equilibrated state
state = simulation.context.getState(getPositions=True, getVelocities=True)
with open("{system_dir}/equilibrated.pdb", 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

# Save checkpoint
simulation.saveCheckpoint("{system_dir}/equilibrated.chk")

print("\\n✓ Equilibration complete!")
print(f"  Saved equilibrated structure: {system_dir}/equilibrated.pdb")
print(f"  Saved checkpoint: {system_dir}/equilibrated.chk")
'''
    
    script_path = system_dir / "02_equilibrate.py"
    script_path.write_text(script)
    return script_path


def generate_production_script(system_dir: Path, production_ns: float = 10.0) -> Path:
    """Generate production MD script"""
    production_steps = int(production_ns * 500000)  # ns to steps at 2fs
    output_freq = 5000  # every 10 ps
    
    script = f'''#!/usr/bin/env python3
"""
Production MD Simulation for {system_dir.name}
Duration: {production_ns} ns
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import time

print("="*60)
print("PRODUCTION MD: {system_dir.name}")
print("Duration: {production_ns} ns")
print("="*60)

# Load equilibrated structure
print("\\nLoading equilibrated structure...")
pdb = app.PDBFile("{system_dir}/equilibrated.pdb")

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
    properties = {{'CudaPrecision': 'mixed'}}
    print("Using CUDA platform")
except:
    try:
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {{}}
        print("Using OpenCL platform")
    except:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {{}}
        print("Using CPU platform (will be slow)")

# Create simulation
simulation = app.Simulation(
    pdb.topology,
    system,
    integrator,
    platform,
    properties if properties else {{}}
)
simulation.context.setPositions(pdb.positions)

# Load checkpoint if available
checkpoint_file = "{system_dir}/equilibrated.chk"
try:
    simulation.loadCheckpoint(checkpoint_file)
    print("Loaded checkpoint from equilibration")
except:
    simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)
    print("Starting with fresh velocities")

# Production reporters
simulation.reporters.append(
    app.DCDReporter('{system_dir}/trajectory.dcd', {output_freq})
)
simulation.reporters.append(
    app.StateDataReporter(
        '{system_dir}/production.log',
        {output_freq},
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
    app.CheckpointReporter('{system_dir}/production.chk', {output_freq * 10})
)

# Run production
production_steps = {production_steps}
print(f"\\nRunning {{production_steps}} steps ({production_ns} ns)...")
start_time = time.time()

# Run in chunks to report progress
chunk_size = 25000  # ~50 ps per chunk
n_chunks = production_steps // chunk_size

for i in range(n_chunks):
    simulation.step(chunk_size)
    elapsed = time.time() - start_time
    progress = (i + 1) / n_chunks * 100
    time_ns = (i + 1) * chunk_size * 0.004 / 1000  # Convert to ns
    print(f"  Progress: {{progress:.1f}}% ({{time_ns:.2f}} ns) - Elapsed: {{elapsed/60:.1f}} min")

# Final steps
remaining = production_steps % chunk_size
if remaining > 0:
    simulation.step(remaining)

# Save final state
simulation.saveState("{system_dir}/final_state.xml")
state = simulation.context.getState(getPositions=True)
with open("{system_dir}/final.pdb", 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

total_time = time.time() - start_time
print(f"\\n✓ Production complete!")
print(f"  Total wall time: {{total_time/3600:.2f}} hours")
print(f"  Output: {system_dir}/trajectory.dcd")
'''
    
    script_path = system_dir / "03_production.py"
    script_path.write_text(script)
    return script_path


def generate_analysis_script(system_dir: Path, binding_residues: List[int]) -> Path:
    """Generate trajectory analysis script"""
    script = f'''#!/usr/bin/env python3
"""
Trajectory Analysis for {system_dir.name}
Analyzes binding interactions and structural stability
"""
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, distances
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json

print("="*60)
print("TRAJECTORY ANALYSIS: {system_dir.name}")
print("="*60)

# Load trajectory
print("\\nLoading trajectory...")
try:
    u = mda.Universe(
        "{system_dir}/minimized.pdb",  # Topology
        "{system_dir}/trajectory.dcd"   # Trajectory
    )
    print(f"  Frames: {{len(u.trajectory)}}")
    print(f"  Atoms: {{u.atoms.n_atoms}}")
except Exception as e:
    print(f"Error loading trajectory: {{e}}")
    print("Using final structure for static analysis...")
    u = mda.Universe("{system_dir}/final.pdb")

# Define selections
protein = u.select_atoms("protein")
backbone = u.select_atoms("protein and backbone")
binding_site = u.select_atoms("protein and resid {' '.join(map(str, binding_residues))}")

print(f"\\nSelections:")
print(f"  Protein atoms: {{protein.n_atoms}}")
print(f"  Backbone atoms: {{backbone.n_atoms}}")
print(f"  Binding site atoms: {{binding_site.n_atoms}}")
print(f"  Binding residues: {binding_residues}")

# ========================================
# RMSD Analysis
# ========================================
print("\\nCalculating RMSD...")

if hasattr(u.trajectory, '__len__') and len(u.trajectory) > 1:
    # Align trajectory to first frame
    aligner = align.AlignTraj(u, u, select='backbone', in_memory=True)
    aligner.run()
    
    # Calculate RMSD
    rmsd_data = rms.RMSD(u, u, select='backbone', ref_frame=0)
    rmsd_data.run()
    
    rmsd_values = rmsd_data.results.rmsd[:, 2]  # RMSD column
    time_ps = rmsd_data.results.rmsd[:, 1]      # Time column
    
    print(f"  Average RMSD: {{rmsd_values.mean():.2f}} Å")
    print(f"  Max RMSD: {{rmsd_values.max():.2f}} Å")
    print(f"  Final RMSD: {{rmsd_values[-1]:.2f}} Å")
    
    # ========================================
    # Binding Site RMSD
    # ========================================
    print("\\nCalculating binding site RMSD...")
    bs_rmsd = rms.RMSD(u, u, select='resid {' '.join(map(str, binding_residues))} and backbone', ref_frame=0)
    bs_rmsd.run()
    bs_rmsd_values = bs_rmsd.results.rmsd[:, 2]
    print(f"  Average: {{bs_rmsd_values.mean():.2f}} Å")
    
    # ========================================
    # Plotting
    # ========================================
    print("\\nGenerating plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Backbone RMSD
    ax1 = axes[0, 0]
    ax1.plot(time_ps/1000, rmsd_values)
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("RMSD (Å)")
    ax1.set_title("Backbone RMSD")
    ax1.axhline(y=rmsd_values.mean(), color='r', linestyle='--', alpha=0.5)
    
    # Binding site RMSD
    ax2 = axes[0, 1]
    ax2.plot(time_ps/1000, bs_rmsd_values)
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("RMSD (Å)")
    ax2.set_title("Binding Site RMSD")
    
    # RMSD distribution
    ax3 = axes[1, 0]
    ax3.hist(rmsd_values, bins=50, alpha=0.7, label="Backbone")
    ax3.hist(bs_rmsd_values, bins=50, alpha=0.7, label="Binding Site")
    ax3.set_xlabel("RMSD (Å)")
    ax3.set_ylabel("Count")
    ax3.set_title("RMSD Distribution")
    ax3.legend()
    
    # Energy from log file (if available)
    ax4 = axes[1, 1]
    try:
        import pandas as pd
        log_data = pd.read_csv("{system_dir}/production.log", sep=',', skiprows=1, 
                               names=['step', 'time', 'PE', 'KE', 'TE', 'T', 'density', 'speed'])
        ax4.plot(log_data['time']/1000, log_data['TE'])
        ax4.set_xlabel("Time (ns)")
        ax4.set_ylabel("Total Energy (kJ/mol)")
        ax4.set_title("Total Energy")
    except:
        ax4.text(0.5, 0.5, "Energy data not available", ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig("{system_dir}/analysis_results.png", dpi=150)
    print(f"  Saved: {system_dir}/analysis_results.png")
    
    # ========================================
    # Save Results
    # ========================================
    results = {{
        "system": "{system_dir.name}",
        "n_frames": len(u.trajectory),
        "backbone_rmsd_mean": float(rmsd_values.mean()),
        "backbone_rmsd_std": float(rmsd_values.std()),
        "backbone_rmsd_max": float(rmsd_values.max()),
        "binding_site_rmsd_mean": float(bs_rmsd_values.mean()),
        "binding_site_rmsd_std": float(bs_rmsd_values.std()),
        "binding_residues": {binding_residues}
    }}
    
else:
    print("  Single frame - static analysis only")
    results = {{
        "system": "{system_dir.name}",
        "n_frames": 1,
        "binding_residues": {binding_residues}
    }}

with open("{system_dir}/analysis_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\\n✓ Analysis complete!")
print(f"  Results: {system_dir}/analysis_results.json")
'''
    
    script_path = system_dir / "04_analyze.py"
    script_path.write_text(script)
    return script_path


def generate_master_runner(systems: List[Path]) -> Path:
    """Generate master script to run all simulations"""
    system_dirs = [str(s) for s in systems]
    
    script = f'''#!/usr/bin/env python3
"""
Master Runner for All MD Simulations
=====================================
Runs minimization, equilibration, and production for all systems.
"""
import subprocess
import sys
from pathlib import Path
import time
import json

SYSTEMS = {json.dumps(system_dirs, indent=4)}

def run_phase(system_dir: str, phase: str) -> bool:
    """Run a simulation phase for a system"""
    scripts = {{
        "minimize": "01_minimize.py",
        "equilibrate": "02_equilibrate.py",
        "production": "03_production.py",
        "analyze": "04_analyze.py"
    }}
    
    script = Path(system_dir) / scripts[phase]
    if not script.exists():
        print(f"  ✗ Script not found: {{script}}")
        return False
    
    print(f"\\n>>> Running {{phase}} for {{Path(system_dir).name}}...")
    try:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=system_dir,
            capture_output=True,
            text=True,
            timeout=3600 * 24  # 24 hour timeout
        )
        
        if result.returncode == 0:
            print(f"  ✓ {{phase}} complete")
            return True
        else:
            print(f"  ✗ {{phase}} failed:")
            print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ✗ {{phase}} timed out")
        return False
    except Exception as e:
        print(f"  ✗ {{phase}} error: {{e}}")
        return False


def main():
    print("="*70)
    print("MOLECULAR DYNAMICS SIMULATION RUNNER")
    print("="*70)
    print(f"\\nSystems to simulate: {{len(SYSTEMS)}}")
    
    results = {{}}
    
    for system_dir in SYSTEMS:
        system_name = Path(system_dir).name
        print(f"\\n" + "="*70)
        print(f"SYSTEM: {{system_name}}")
        print("="*70)
        
        results[system_name] = {{}}
        
        # Run each phase
        for phase in ["minimize", "equilibrate", "production", "analyze"]:
            success = run_phase(system_dir, phase)
            results[system_name][phase] = "success" if success else "failed"
            
            if not success and phase in ["minimize", "equilibrate"]:
                print(f"  Skipping remaining phases for {{system_name}}")
                break
    
    # Summary
    print("\\n" + "="*70)
    print("SIMULATION SUMMARY")
    print("="*70)
    
    for system, phases in results.items():
        status = "✓" if all(v == "success" for v in phases.values()) else "✗"
        print(f"  {{status}} {{system}}: {{phases}}")
    
    # Save results
    with open("{WORK_DIR}/simulation_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\\nResults saved to: {WORK_DIR}/simulation_results.json")


if __name__ == "__main__":
    main()
'''
    
    runner_path = WORK_DIR / "run_all_simulations.py"
    runner_path.write_text(script)
    return runner_path


# ============================================================
# MAIN SETUP
# ============================================================

def setup_all_simulations():
    """Set up all receptor-ligand simulations"""
    print("="*70)
    print("MOLECULAR DYNAMICS SIMULATION SETUP")
    print("="*70)
    
    # Check OpenMM
    openmm_ok, platform = check_openmm_available()
    print(f"\nOpenMM: {'✓ Available' if openmm_ok else '✗ Not available'} (Platform: {platform})")
    
    if not openmm_ok:
        print("\nInstall OpenMM with: pip install openmm")
        return
    
    # Create directories
    structures_dir = WORK_DIR / "structures"
    ligands_dir = WORK_DIR / "ligands"
    structures_dir.mkdir(parents=True, exist_ok=True)
    ligands_dir.mkdir(parents=True, exist_ok=True)
    
    # Download receptor structures
    print("\n" + "-"*50)
    print("DOWNLOADING RECEPTOR STRUCTURES")
    print("-"*50)
    
    for receptor_name, info in RECEPTOR_STRUCTURES.items():
        pdb_path = structures_dir / f"{info['pdb_id']}.pdb"
        download_pdb(info['pdb_id'], pdb_path)
    
    # Prepare ligands
    print("\n" + "-"*50)
    print("PREPARING LIGANDS")
    print("-"*50)
    
    ligand_files = {}
    for ligand_name, smiles in LIGAND_SMILES.items():
        lig_dir = ligands_dir / ligand_name
        result = prepare_ligand_3d(ligand_name, smiles, lig_dir)
        if result:
            ligand_files[ligand_name] = result
    
    # Set up simulation systems
    print("\n" + "-"*50)
    print("SETTING UP SIMULATION SYSTEMS")
    print("-"*50)
    
    system_dirs = []
    
    for receptor_name, ligand_name in SIMULATION_PAIRS:
        print(f"\n  {receptor_name} + {ligand_name}")
        
        receptor_info = RECEPTOR_STRUCTURES.get(receptor_name)
        if not receptor_info:
            print(f"    ✗ Unknown receptor: {receptor_name}")
            continue
        
        if ligand_name not in ligand_files:
            print(f"    ✗ Ligand not prepared: {ligand_name}")
            continue
        
        system_dir = WORK_DIR / f"{receptor_name}_{ligand_name}"
        system_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy/link structure files
        receptor_pdb = structures_dir / f"{receptor_info['pdb_id']}.pdb"
        
        # Generate scripts
        generate_minimization_script(system_dir, str(receptor_pdb), str(ligand_files[ligand_name]))
        generate_equilibration_script(system_dir)
        generate_production_script(system_dir, production_ns=10.0)  # 10 ns production
        generate_analysis_script(system_dir, receptor_info['binding_site'])
        
        system_dirs.append(system_dir)
        print(f"    ✓ Generated scripts in {system_dir}")
    
    # Generate master runner
    print("\n" + "-"*50)
    print("GENERATING MASTER RUNNER")
    print("-"*50)
    runner_path = generate_master_runner(system_dirs)
    print(f"  ✓ {runner_path}")
    
    # Summary
    print("\n" + "="*70)
    print("SETUP COMPLETE")
    print("="*70)
    print(f"""
Systems prepared: {len(system_dirs)}
Work directory: {WORK_DIR}

To run all simulations:
  python {runner_path}

Or run individual systems:
  cd {WORK_DIR}/<system>
  python 01_minimize.py
  python 02_equilibrate.py  
  python 03_production.py
  python 04_analyze.py

Estimated time per system:
  - Minimization: ~5-10 min
  - Equilibration: ~30-60 min
  - Production (10 ns): ~4-8 hours (GPU) / 24-48 hours (CPU)
  
Note: For longer simulations, modify production_ns in 03_production.py
""")
    
    return system_dirs


if __name__ == "__main__":
    setup_all_simulations()
