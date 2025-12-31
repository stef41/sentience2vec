"""
Molecular Dynamics Simulation Setup for Receptor-Ligand Studies
Uses OpenMM/GROMACS for atomistic simulations
"""

import os
import json
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from enum import Enum


class SimulationType(Enum):
    EQUILIBRATION = "equilibration"
    PRODUCTION = "production"
    MINIMIZATION = "minimization"
    HEATING = "heating"
    COOLING = "cooling"


class ForceField(Enum):
    AMBER_FF14SB = "amber14-all.xml"        # Protein
    AMBER_GAFF2 = "gaff2.xml"               # Small molecules
    CHARMM36 = "charmm36.xml"               # Alternative
    TIP3P = "tip3p.xml"                     # Water model
    TIP4PEW = "tip4pew.xml"                 # Better water model


@dataclass
class SimulationParameters:
    """Parameters for MD simulation"""
    temperature_k: float = 310.0            # Body temperature
    pressure_atm: float = 1.0
    timestep_fs: float = 2.0                # 2 femtoseconds
    total_time_ns: float = 100.0            # 100 nanoseconds
    
    # Output settings
    output_frequency_ps: float = 10.0       # Save every 10 ps
    
    # System setup
    box_padding_nm: float = 1.2             # Water padding
    ionic_strength_molar: float = 0.15      # Physiological salt
    
    # Constraints
    hydrogen_mass_repartitioning: bool = True  # Allow 4fs timestep
    constrain_bonds: str = "all"


@dataclass
class ReceptorStructure:
    """Receptor structural information"""
    name: str
    pdb_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    alphafold_id: Optional[str] = None
    
    # Key residues
    binding_site_residues: List[int] = field(default_factory=list)
    active_site_residues: List[int] = field(default_factory=list)
    
    # Structure file paths
    pdb_file: Optional[str] = None
    topology_file: Optional[str] = None


# Key receptor structures for consciousness research
RECEPTOR_STRUCTURES = {
    "5-HT2A": ReceptorStructure(
        name="Serotonin 2A Receptor",
        pdb_id="6WHA",          # LSD-bound structure!
        uniprot_id="P28223",
        binding_site_residues=[151, 155, 159, 207, 229, 230, 234, 313, 336, 340, 366],
        # Key residues: D155 (salt bridge), S159 (H-bond), F340 (pi-stacking)
    ),
    
    "5-HT2B": ReceptorStructure(
        name="Serotonin 2B Receptor", 
        pdb_id="5TVN",          # Ergotamine-bound
        uniprot_id="P41595",
        binding_site_residues=[134, 135, 136, 139, 207, 217, 313, 317, 340, 344]
    ),
    
    "D2": ReceptorStructure(
        name="Dopamine D2 Receptor",
        pdb_id="6CM4",          # Risperidone-bound
        uniprot_id="P14416",
        binding_site_residues=[114, 115, 118, 119, 194, 197, 379, 380, 383, 386, 389]
    ),
    
    "MOR": ReceptorStructure(
        name="Mu Opioid Receptor",
        pdb_id="5C1M",          # Activated structure
        uniprot_id="P35372",
        binding_site_residues=[148, 149, 202, 230, 233, 289, 293, 322, 325, 328]
    ),
    
    "KOR": ReceptorStructure(
        name="Kappa Opioid Receptor",
        pdb_id="4DJH",          # JDTic-bound
        uniprot_id="P41145",
        binding_site_residues=[138, 139, 141, 192, 219, 223, 280, 283, 284, 312]
    ),
    
    "NMDA": ReceptorStructure(
        name="NMDA Receptor (GluN1-GluN2B)",
        pdb_id="4PE5",          # Ifenprodil-bound
        uniprot_id="Q00959",    # GluN1
        binding_site_residues=[518, 531, 532, 786, 788, 790, 810, 811, 813]
    ),
    
    "CB1": ReceptorStructure(
        name="Cannabinoid CB1 Receptor",
        pdb_id="5TGZ",          # THC analogue-bound
        uniprot_id="P21554",
        binding_site_residues=[178, 182, 192, 193, 196, 274, 277, 278, 353, 356, 360]
    )
}


# Common psychoactive ligand SMILES
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


class MDSimulationSetup:
    """
    Set up molecular dynamics simulations for receptor-ligand binding studies
    """
    
    def __init__(self, work_dir: str = "./md_simulations"):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.receptor_structures = RECEPTOR_STRUCTURES
        self.ligand_smiles = LIGAND_SMILES
    
    def fetch_receptor_structure(
        self,
        receptor_name: str,
        source: str = "rcsb"
    ) -> Optional[str]:
        """
        Fetch receptor structure from database
        Returns path to downloaded PDB file
        """
        structure = self.receptor_structures.get(receptor_name)
        if not structure or not structure.pdb_id:
            return None
        
        receptor_dir = self.work_dir / receptor_name
        receptor_dir.mkdir(exist_ok=True)
        
        pdb_path = receptor_dir / f"{structure.pdb_id}.pdb"
        
        if not pdb_path.exists():
            # Download from RCSB PDB
            import urllib.request
            url = f"https://files.rcsb.org/download/{structure.pdb_id}.pdb"
            
            try:
                urllib.request.urlretrieve(url, pdb_path)
                print(f"Downloaded {structure.pdb_id} to {pdb_path}")
            except Exception as e:
                print(f"Error downloading {structure.pdb_id}: {e}")
                return None
        
        return str(pdb_path)
    
    def prepare_ligand(
        self,
        ligand_name: str,
        output_dir: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Prepare ligand for simulation:
        1. Generate 3D coordinates from SMILES
        2. Assign partial charges
        3. Generate GAFF2 parameters
        
        Returns dict with file paths
        """
        smiles = self.ligand_smiles.get(ligand_name)
        if not smiles:
            raise ValueError(f"Unknown ligand: {ligand_name}")
        
        if output_dir is None:
            output_dir = self.work_dir / "ligands" / ligand_name
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # This would use RDKit for coordinate generation
        script = f'''
# Ligand preparation script for {ligand_name}
# Requires: rdkit, openmm, openff-toolkit

from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField

# Generate 3D structure from SMILES
smiles = "{smiles}"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

# Save to files
Chem.MolToMolFile(mol, "{output_dir}/ligand.mol")
Chem.MolToPDBFile(mol, "{output_dir}/ligand.pdb")

# Generate parameters with OpenFF
offmol = Molecule.from_rdkit(mol)
ff = ForceField("openff-2.0.0.offxml")
# topology = offmol.to_topology()
# system = ff.create_openmm_system(topology)

print(f"Ligand prepared: {output_dir}")
'''
        
        script_path = output_dir / "prepare_ligand.py"
        with open(script_path, "w") as f:
            f.write(script)
        
        return {
            "smiles": smiles,
            "script": str(script_path),
            "output_dir": str(output_dir),
            "mol_file": str(output_dir / "ligand.mol"),
            "pdb_file": str(output_dir / "ligand.pdb")
        }
    
    def setup_system(
        self,
        receptor_name: str,
        ligand_name: str,
        params: SimulationParameters = None
    ) -> Dict[str, str]:
        """
        Set up complete receptor-ligand system for simulation
        """
        if params is None:
            params = SimulationParameters()
        
        system_dir = self.work_dir / f"{receptor_name}_{ligand_name}"
        system_dir.mkdir(exist_ok=True)
        
        receptor = self.receptor_structures.get(receptor_name)
        
        # Generate OpenMM setup script
        script = f'''
"""
OpenMM Simulation Setup: {receptor_name} + {ligand_name}
Generated automatically for consciousness research
"""

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import GAFFTemplateGenerator
import numpy as np

# ========================================
# System Parameters
# ========================================
TEMPERATURE = {params.temperature_k} * unit.kelvin
PRESSURE = {params.pressure_atm} * unit.atmospheres
TIMESTEP = {params.timestep_fs} * unit.femtoseconds
PRODUCTION_TIME = {params.total_time_ns} * unit.nanoseconds
OUTPUT_FREQ_PS = {params.output_frequency_ps}
BOX_PADDING = {params.box_padding_nm} * unit.nanometers
IONIC_STRENGTH = {params.ionic_strength_molar} * unit.molar

# ========================================
# Load Structures
# ========================================

# Load receptor
print("Loading receptor structure...")
receptor_pdb = app.PDBFile("{receptor.pdb_id}.pdb")

# Load ligand (prepared separately)
print("Loading ligand...")
ligand_mol = app.PDBFile("../ligands/{ligand_name}/ligand.pdb")

# ========================================
# Build System
# ========================================

# Create modeller to combine structures
modeller = app.Modeller(receptor_pdb.topology, receptor_pdb.positions)
modeller.add(ligand_mol.topology, ligand_mol.positions)

# Add solvent
print("Adding solvent...")
modeller.addSolvent(
    forcefield=None,  # Will use TIP3P
    model='tip3p',
    padding=BOX_PADDING,
    ionicStrength=IONIC_STRENGTH
)

# Set up force field with GAFF for ligand
forcefield = app.ForceField(
    'amber14-all.xml',
    'amber14/tip3p.xml'
)

# Add GAFF parameters for ligand
gaff = GAFFTemplateGenerator(
    molecules='{LIGAND_SMILES.get(ligand_name, "")}'
)
forcefield.registerTemplateGenerator(gaff.generator)

# Create system
print("Creating system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    hydrogenMass=4.0 * unit.amu if {params.hydrogen_mass_repartitioning} else None
)

# Add barostat for NPT
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

# ========================================
# Integrator and Simulation
# ========================================

integrator = mm.LangevinMiddleIntegrator(
    TEMPERATURE,
    1.0 / unit.picoseconds,  # Friction
    TIMESTEP
)

# Create simulation
platform = mm.Platform.getPlatformByName('CUDA')
properties = {{'CudaPrecision': 'mixed'}}

simulation = app.Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties
)
simulation.context.setPositions(modeller.positions)

# ========================================
# Minimization
# ========================================

print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=5000)

# ========================================
# Equilibration
# ========================================

print("Equilibrating (heating)...")
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

# Report during equilibration
eq_reporters = [
    app.StateDataReporter(
        'equilibration.log', 1000,
        step=True, time=True, potentialEnergy=True,
        temperature=True, density=True
    )
]
for r in eq_reporters:
    simulation.reporters.append(r)

# 1 ns equilibration
eq_steps = int(1.0 * unit.nanoseconds / TIMESTEP)
simulation.step(eq_steps)
simulation.reporters.clear()

# ========================================
# Production
# ========================================

print("Running production simulation...")

# Set up production reporters
output_freq = int(OUTPUT_FREQ_PS * unit.picoseconds / TIMESTEP)

simulation.reporters.append(
    app.DCDReporter('trajectory.dcd', output_freq)
)
simulation.reporters.append(
    app.StateDataReporter(
        'production.log', output_freq,
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, totalEnergy=True,
        temperature=True, density=True, speed=True
    )
)
simulation.reporters.append(
    app.CheckpointReporter('checkpoint.chk', output_freq * 10)
)

# Run production
prod_steps = int(PRODUCTION_TIME / TIMESTEP)
simulation.step(prod_steps)

# Save final state
simulation.saveState('final_state.xml')
print("Simulation complete!")

# ========================================
# Key Binding Site Residues to Analyze
# ========================================
BINDING_SITE_RESIDUES = {receptor.binding_site_residues}

# These residues should be analyzed for:
# - Hydrogen bond occupancy
# - Salt bridge distances  
# - Pi-stacking interactions
# - RMSD of binding pose
'''
        
        script_path = system_dir / "run_simulation.py"
        with open(script_path, "w") as f:
            f.write(script)
        
        # Also generate GROMACS input files
        self._generate_gromacs_mdp(system_dir, params)
        
        return {
            "system_dir": str(system_dir),
            "openmm_script": str(script_path),
            "receptor_pdb": receptor.pdb_id,
            "binding_site_residues": receptor.binding_site_residues
        }
    
    def _generate_gromacs_mdp(
        self,
        output_dir: Path,
        params: SimulationParameters
    ) -> None:
        """Generate GROMACS .mdp parameter files"""
        
        # Energy minimization
        em_mdp = f"""; Energy Minimization
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

; Output
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Neighbor searching
cutoff-scheme = Verlet
nstlist     = 10
ns_type     = grid
pbc         = xyz
rlist       = 1.0

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0

; VdW
vdwtype     = Cut-off
rvdw        = 1.0
"""
        
        # Production run
        md_mdp = f"""; Production MD
integrator  = md
dt          = {params.timestep_fs / 1000}  ; ps
nsteps      = {int(params.total_time_ns * 1e6 / params.timestep_fs)}

; Output control
nstxout-compressed = {int(params.output_frequency_ps * 1000 / params.timestep_fs)}
nstenergy   = {int(params.output_frequency_ps * 1000 / params.timestep_fs)}
nstlog      = {int(params.output_frequency_ps * 1000 / params.timestep_fs)}

; Bond constraints
constraint_algorithm = LINCS
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
nstlist     = 20
ns_type     = grid
pbc         = xyz
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2
pme_order   = 4
fourierspacing = 0.16

; VdW
vdwtype     = Cut-off
rvdw        = 1.2
DispCorr    = EnerPres

; Temperature coupling
tcoupl      = V-rescale
tc_grps     = Protein_LIG Water_and_ions
tau_t       = 0.1 0.1
ref_t       = {params.temperature_k} {params.temperature_k}

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = {params.pressure_atm}
compressibility = 4.5e-5
refcoord_scaling = com

; Velocity generation
gen_vel     = no
continuation = yes
"""
        
        (output_dir / "em.mdp").write_text(em_mdp)
        (output_dir / "md_prod.mdp").write_text(md_mdp)
    
    def setup_docking(
        self,
        receptor_name: str,
        ligand_name: str
    ) -> Dict[str, str]:
        """
        Set up molecular docking with AutoDock Vina
        Useful for initial pose generation before MD
        """
        receptor = self.receptor_structures.get(receptor_name)
        if not receptor:
            raise ValueError(f"Unknown receptor: {receptor_name}")
        
        docking_dir = self.work_dir / "docking" / f"{receptor_name}_{ligand_name}"
        docking_dir.mkdir(parents=True, exist_ok=True)
        
        # Calculate box center from binding site residues
        # (In practice, would be computed from actual coordinates)
        
        vina_config = f"""# AutoDock Vina configuration
# {receptor_name} + {ligand_name}

receptor = {receptor.pdb_id}_prepared.pdbqt
ligand = {ligand_name}.pdbqt

out = {ligand_name}_docked.pdbqt
log = docking.log

# Box centered on binding site
# Residues: {receptor.binding_site_residues}
center_x = 0.0
center_y = 0.0  
center_z = 0.0
size_x = 25
size_y = 25
size_z = 25

exhaustiveness = 32
num_modes = 20
energy_range = 5
"""
        
        config_path = docking_dir / "vina_config.txt"
        config_path.write_text(vina_config)
        
        # Preparation script
        prep_script = f'''#!/bin/bash
# Prepare structures for docking

# Receptor preparation (requires MGLTools)
prepare_receptor4.py -r {receptor.pdb_id}.pdb -o {receptor.pdb_id}_prepared.pdbqt

# Ligand preparation
obabel ligand.mol -O {ligand_name}.pdbqt -h

# Run docking
vina --config vina_config.txt

# Extract best pose
vina_split --input {ligand_name}_docked.pdbqt
'''
        
        script_path = docking_dir / "run_docking.sh"
        script_path.write_text(prep_script)
        
        return {
            "docking_dir": str(docking_dir),
            "config_file": str(config_path),
            "run_script": str(script_path)
        }


class TrajectoryAnalyzer:
    """Analyze MD trajectories for binding interactions"""
    
    def __init__(self, trajectory_path: str, topology_path: str):
        self.trajectory_path = trajectory_path
        self.topology_path = topology_path
    
    def generate_analysis_script(self, binding_residues: List[int]) -> str:
        """Generate MDAnalysis script for trajectory analysis"""
        
        script = f'''
"""
Trajectory Analysis for Receptor-Ligand Binding
Uses MDAnalysis library
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances, hydrogen_bonds
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory
u = mda.Universe("{self.topology_path}", "{self.trajectory_path}")

# Define selections
receptor = u.select_atoms("protein")
ligand = u.select_atoms("resname LIG")  # Adjust resname as needed
binding_site = u.select_atoms("protein and resid {' '.join(map(str, binding_residues))}")

# ========================================
# RMSD Analysis
# ========================================

print("Calculating RMSD...")
rmsd_analysis = rms.RMSD(
    u,
    select="backbone",  # Protein backbone
    ref_frame=0
)
rmsd_analysis.run()

# Ligand RMSD
ligand_rmsd = rms.RMSD(
    u,
    select="resname LIG",
    ref_frame=0
)
ligand_rmsd.run()

# ========================================
# Binding Site Contacts
# ========================================

print("Analyzing binding contacts...")
contact_distances = []

for ts in u.trajectory:
    distances_frame = distances.distance_array(
        ligand.positions,
        binding_site.positions
    )
    min_dist = distances_frame.min()
    contact_distances.append(min_dist)

contact_distances = np.array(contact_distances)

# ========================================
# Hydrogen Bond Analysis
# ========================================

print("Analyzing hydrogen bonds...")
hbonds = hydrogen_bonds.HydrogenBondAnalysis(
    u,
    donors_sel="resname LIG",
    acceptors_sel=f"protein and resid {' '.join(map(str, binding_residues))}",
    d_h_cutoff=1.2,
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=150.0
)
hbonds.run()

hbond_count = hbonds.count_by_time()

# ========================================
# Binding Free Energy (MM-GBSA style)
# ========================================
# Note: This is a simplified version
# For accurate free energies, use proper FEP/TI methods

print("Estimating binding interactions...")

# This would require energy calculation plugins
# interaction_energies = calculate_interaction_energy(receptor, ligand)

# ========================================
# Output Results
# ========================================

print("\\n=== Analysis Results ===")
print(f"Simulation time: {{u.trajectory.n_frames * u.trajectory.dt / 1000:.1f}} ns")
print(f"Average RMSD: {{rmsd_analysis.results.rmsd[:, 2].mean():.2f}} Å")
print(f"Ligand RMSD: {{ligand_rmsd.results.rmsd[:, 2].mean():.2f}} Å")
print(f"Average contact distance: {{contact_distances.mean():.2f}} Å")
print(f"Total H-bonds: {{hbond_count.sum():.0f}}")
print(f"Avg H-bonds per frame: {{hbond_count.mean():.2f}}")

# ========================================
# Plotting
# ========================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Protein RMSD
ax1 = axes[0, 0]
time_ps = rmsd_analysis.results.rmsd[:, 1]
ax1.plot(time_ps/1000, rmsd_analysis.results.rmsd[:, 2])
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_title("Protein Backbone RMSD")

# Ligand RMSD
ax2 = axes[0, 1]
ax2.plot(time_ps/1000, ligand_rmsd.results.rmsd[:, 2])
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("RMSD (Å)")
ax2.set_title("Ligand RMSD")

# Contact distance
ax3 = axes[1, 0]
ax3.plot(time_ps/1000, contact_distances)
ax3.set_xlabel("Time (ns)")
ax3.set_ylabel("Min Distance (Å)")
ax3.set_title("Ligand-Binding Site Distance")
ax3.axhline(y=3.5, color='r', linestyle='--', label="Contact threshold")
ax3.legend()

# H-bond count
ax4 = axes[1, 1]
ax4.plot(time_ps/1000, hbond_count)
ax4.set_xlabel("Time (ns)")
ax4.set_ylabel("H-bond Count")
ax4.set_title("Hydrogen Bonds Over Time")

plt.tight_layout()
plt.savefig("binding_analysis.png", dpi=150)
plt.show()

print("\\nAnalysis complete! Saved to binding_analysis.png")
'''
        return script


if __name__ == "__main__":
    # Example: Set up LSD-5HT2A simulation
    setup = MDSimulationSetup("./md_simulations")
    
    print("=== Molecular Dynamics Setup ===\n")
    
    # 1. Prepare ligand
    print("1. Preparing LSD ligand...")
    ligand_info = setup.prepare_ligand("LSD")
    print(f"   SMILES: {ligand_info['smiles'][:50]}...")
    print(f"   Output: {ligand_info['output_dir']}")
    
    # 2. Set up complete system
    print("\n2. Setting up 5-HT2A + LSD system...")
    params = SimulationParameters(
        temperature_k=310.0,
        total_time_ns=100.0,
        output_frequency_ps=10.0
    )
    system_info = setup.setup_system("5-HT2A", "LSD", params)
    print(f"   System directory: {system_info['system_dir']}")
    print(f"   Receptor PDB: {system_info['receptor_pdb']}")
    print(f"   Binding site residues: {system_info['binding_site_residues']}")
    
    # 3. Set up docking
    print("\n3. Setting up molecular docking...")
    docking_info = setup.setup_docking("5-HT2A", "LSD")
    print(f"   Docking directory: {docking_info['docking_dir']}")
    
    # 4. Generate analysis script
    print("\n4. Generating trajectory analysis script...")
    analyzer = TrajectoryAnalyzer(
        "trajectory.dcd",
        "system.pdb"
    )
    analysis_script = analyzer.generate_analysis_script(
        RECEPTOR_STRUCTURES["5-HT2A"].binding_site_residues
    )
    print("   Analysis script ready for post-simulation processing")
    
    print("\n=== Available Receptor Structures ===")
    for name, struct in RECEPTOR_STRUCTURES.items():
        print(f"  {name}: PDB {struct.pdb_id} ({struct.name})")
