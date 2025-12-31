#!/usr/bin/env python3
"""
Prepare ALL receptors with PDB structures.
Downloads, cleans, adds hydrogens, and minimizes receptor structures.
"""

import os
import sys
import json
import urllib.request
from pathlib import Path

# Import the full receptor database
sys.path.insert(0, str(Path(__file__).parent))
from receptor_pdb_database_full import RECEPTOR_PDB_DATABASE as FULL_RECEPTOR_DATABASE

BASE_DIR = Path("/data/users/zacharie/whc/consciousness_study/md_simulations")
OUTPUT_DIR = BASE_DIR / "receptors_full"
OUTPUT_DIR.mkdir(exist_ok=True)

def download_pdb(pdb_id, outpath):
    """Download PDB from RCSB"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, outpath)
        return outpath.exists() and outpath.stat().st_size > 1000
    except Exception as e:
        print(f"    Download error: {e}")
        return False

def clean_receptor(pdb_path, output_path):
    """Clean receptor - remove ligands, waters, keep only protein"""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        
        fixer = PDBFixer(filename=str(pdb_path))
        
        # Find and remove non-standard residues (ligands)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)  # Remove ligands and water
        
        # Add missing atoms
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        # Add hydrogens
        fixer.addMissingHydrogens(7.4)  # pH 7.4
        
        # Write cleaned structure
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        
        return True
    except Exception as e:
        print(f"    Cleaning error: {e}")
        return False

def minimize_receptor(pdb_path, output_path, steps=500):
    """Quick energy minimization"""
    try:
        from openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
        from openmm import LangevinMiddleIntegrator, Platform
        from openmm.unit import kelvin, nanometer, picoseconds
        import openmm
        
        # Load
        pdb = PDBFile(str(pdb_path))
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        # Create modeller (no solvation for speed)
        modeller = Modeller(pdb.topology, pdb.positions)
        
        # Create system with implicit solvent (faster)
        try:
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=PME,
                nonbondedCutoff=1.0*nanometer,
                constraints=HBonds
            )
        except:
            # If PME fails, use no cutoff (vacuum)
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=openmm.app.NoCutoff,
                constraints=HBonds
            )
        
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 0.002*picoseconds)
        
        # Try GPU first
        try:
            platform = Platform.getPlatformByName('OpenCL')
            simulation = openmm.app.Simulation(modeller.topology, system, integrator, platform)
        except:
            platform = Platform.getPlatformByName('CPU')
            simulation = openmm.app.Simulation(modeller.topology, system, integrator, platform)
        
        simulation.context.setPositions(modeller.positions)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=steps)
        
        # Save
        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(output_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, positions, f)
        
        return True
    except Exception as e:
        print(f"    Minimization error: {e}")
        return False

def prepare_receptor(name, pdb_id, receptor_dir):
    """Full receptor preparation pipeline"""
    receptor_dir.mkdir(exist_ok=True)
    
    # Paths
    raw_pdb = receptor_dir / f"{name}_raw.pdb"
    cleaned_pdb = receptor_dir / f"{name}_cleaned.pdb"
    minimized_pdb = receptor_dir / f"{name}_minimized.pdb"
    
    # Check if already done
    if minimized_pdb.exists():
        return "cached"
    
    # Download
    if not raw_pdb.exists():
        if not download_pdb(pdb_id, raw_pdb):
            return "download_failed"
    
    # Clean
    if not cleaned_pdb.exists():
        if not clean_receptor(raw_pdb, cleaned_pdb):
            # Try without cleaning (just download)
            cleaned_pdb = raw_pdb
    
    # Minimize (quick, 500 steps)
    if not minimize_receptor(cleaned_pdb, minimized_pdb, steps=500):
        # If minimization fails, keep cleaned version
        if cleaned_pdb.exists() and not minimized_pdb.exists():
            import shutil
            shutil.copy(cleaned_pdb, minimized_pdb)
            return "cleaned_only"
    
    return "success"

def main():
    print("=" * 70)
    print("FULL RECEPTOR PREPARATION")
    print("=" * 70)
    
    # Count receptors with PDB
    receptors_with_pdb = []
    for name, info in FULL_RECEPTOR_DATABASE.items():
        if info.get("pdb_id") and info["pdb_id"] is not None:
            receptors_with_pdb.append((name, info))
    
    print(f"\n📊 Total receptors with PDB: {len(receptors_with_pdb)}")
    print("=" * 70)
    
    results = {"success": [], "cached": [], "cleaned_only": [], "failed": []}
    
    for i, (name, info) in enumerate(receptors_with_pdb, 1):
        pdb_id = info["pdb_id"]
        safe_name = name.replace("/", "_").replace(" ", "_")
        receptor_dir = OUTPUT_DIR / safe_name
        
        print(f"  [{i}/{len(receptors_with_pdb)}] {name} ({pdb_id})...", end=" ", flush=True)
        
        status = prepare_receptor(safe_name, pdb_id, receptor_dir)
        
        if status == "cached":
            print("(cached)")
            results["cached"].append(name)
        elif status == "success":
            print("✓")
            results["success"].append(name)
        elif status == "cleaned_only":
            print("✓ (cleaned only)")
            results["cleaned_only"].append(name)
        else:
            print("✗")
            results["failed"].append(name)
    
    # Summary
    print("\n" + "=" * 70)
    print("RECEPTOR PREPARATION SUMMARY")
    print("=" * 70)
    print(f"  New: {len(results['success'])}")
    print(f"  Cached: {len(results['cached'])}")
    print(f"  Cleaned only: {len(results['cleaned_only'])}")
    print(f"  Failed: {len(results['failed'])}")
    print(f"  TOTAL PREPARED: {len(results['success']) + len(results['cached']) + len(results['cleaned_only'])}")
    
    if results["failed"]:
        print(f"\nFailed receptors:")
        for r in results["failed"]:
            print(f"  - {r}")
    
    # Save results
    with open(BASE_DIR / "receptor_preparation_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✓ Results saved to receptor_preparation_results.json")

if __name__ == "__main__":
    main()
