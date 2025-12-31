#!/usr/bin/env python3
"""
MD Simulation using AlphaFold Structures
=========================================
AlphaFold structures are complete and well-formatted, making them easier to use
for molecular dynamics simulations than experimental PDB structures.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import urllib.request
import json
import sys

# AlphaFold UniProt IDs for key receptors
ALPHAFOLD_RECEPTORS = {
    "5-HT2A": "P28223",  # Serotonin 2A receptor
    "5-HT2B": "P41595",  # Serotonin 2B receptor
    "D2": "P14416",      # Dopamine D2 receptor
    "MOR": "P35372",     # Mu opioid receptor
    "KOR": "P41145",     # Kappa opioid receptor
    "CB1": "P21554",     # Cannabinoid CB1 receptor
}

# Configuration
SYSTEM_DIR = Path(__file__).parent
SYSTEM_NAME = SYSTEM_DIR.name
parts = SYSTEM_NAME.split('_')
RECEPTOR_NAME = parts[0]
LIGAND_NAME = '_'.join(parts[1:])

print("="*60)
print(f"MD SIMULATION SETUP: {SYSTEM_NAME}")
print("="*60)


def download_alphafold_structure(uniprot_id: str, output_dir: Path) -> Path:
    """Download AlphaFold structure"""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    output_path = output_dir / f"AF-{uniprot_id}.pdb"
    
    if output_path.exists():
        print(f"  AlphaFold structure exists: {output_path}")
        return output_path
    
    print(f"  Downloading AlphaFold structure for {uniprot_id}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"  ✓ Downloaded: {output_path}")
        return output_path
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return None


def prepare_receptor_for_md(pdb_path: Path, forcefield: app.ForceField) -> app.Modeller:
    """
    Prepare receptor structure for MD simulation.
    AlphaFold structures are generally well-formed but may need:
    - Truncation to remove disordered regions
    - Addition of hydrogens
    - Solvation
    """
    print("\nPreparing receptor...")
    
    # Load PDB
    pdb = app.PDBFile(str(pdb_path))
    print(f"  Initial atoms: {pdb.topology.getNumAtoms()}")
    
    # Create modeller
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # AlphaFold structures are typically a single chain
    # Keep only residues with good confidence (pLDDT > 70)
    # For now, we'll use the full structure
    
    # Add hydrogens at physiological pH
    print("  Adding hydrogens (pH 7.4)...")
    modeller.addHydrogens(forcefield, pH=7.4)
    print(f"  Atoms after hydrogens: {modeller.topology.getNumAtoms()}")
    
    return modeller


def run_minimization():
    """Run energy minimization"""
    
    # Get receptor info
    uniprot_id = ALPHAFOLD_RECEPTORS.get(RECEPTOR_NAME)
    if not uniprot_id:
        print(f"ERROR: No AlphaFold ID for receptor: {RECEPTOR_NAME}")
        sys.exit(1)
    
    print(f"\nReceptor: {RECEPTOR_NAME}")
    print(f"UniProt ID: {uniprot_id}")
    print(f"Ligand: {LIGAND_NAME}")
    
    # Create structures directory
    structures_dir = SYSTEM_DIR.parent / "alphafold_structures"
    structures_dir.mkdir(exist_ok=True)
    
    # Download AlphaFold structure
    pdb_path = download_alphafold_structure(uniprot_id, structures_dir)
    if not pdb_path:
        sys.exit(1)
    
    # Set up force field
    print("\nSetting up force field...")
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    
    # Prepare receptor
    try:
        modeller = prepare_receptor_for_md(pdb_path, forcefield)
    except Exception as e:
        print(f"ERROR preparing receptor: {e}")
        sys.exit(1)
    
    # Add solvent
    print("\nAdding solvent...")
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=1.0 * unit.nanometers,
        ionicStrength=0.15 * unit.molar
    )
    print(f"  Total atoms: {modeller.topology.getNumAtoms()}")
    
    # Create system
    print("\nCreating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds
    )
    
    # Integrator
    integrator = mm.LangevinMiddleIntegrator(
        310 * unit.kelvin,  # Body temperature
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
    print(f"Platform: {platform_name}")
    
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
    print("Running energy minimization (2000 steps)...")
    simulation.minimizeEnergy(maxIterations=2000)
    
    # Final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy()
    print(f"Final energy: {final_energy}")
    print(f"Energy change: {(final_energy - initial_energy).value_in_unit(unit.kilojoules_per_mole):.1f} kJ/mol")
    
    # Save minimized structure
    output_pdb = SYSTEM_DIR / "minimized.pdb"
    with open(output_pdb, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print(f"\n✓ Saved: {output_pdb}")
    
    # Save system info
    info = {
        "receptor": RECEPTOR_NAME,
        "uniprot_id": uniprot_id,
        "ligand": LIGAND_NAME,
        "platform": platform_name,
        "initial_energy_kj": initial_energy.value_in_unit(unit.kilojoules_per_mole),
        "final_energy_kj": final_energy.value_in_unit(unit.kilojoules_per_mole),
        "total_atoms": modeller.topology.getNumAtoms()
    }
    with open(SYSTEM_DIR / "system_info.json", 'w') as f:
        json.dump(info, f, indent=2)
    
    print("\n✓ Minimization complete!")
    return simulation


if __name__ == "__main__":
    run_minimization()
