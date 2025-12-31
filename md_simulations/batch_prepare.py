#!/usr/bin/env python3
"""
Batch prepare all receptor structures for MD simulations.
Runs each receptor preparation sequentially.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pathlib import Path
import json
import time
import warnings
warnings.filterwarnings('ignore')


RECEPTORS = {
    '5-HT2A': 'prepared_receptors/5-HT2A/5-HT2A_clean.pdb',
    '5-HT2B': 'prepared_receptors/5-HT2B/5-HT2B_clean.pdb',
    'D2': 'prepared_receptors/D2/D2_clean.pdb',
    'MOR': 'prepared_receptors/MOR/MOR_clean.pdb',
    'KOR': 'prepared_receptors/KOR/KOR_clean.pdb',
    'CB1': 'prepared_receptors/CB1/CB1_clean.pdb',
}


def prepare_receptor(name: str, pdb_path: Path) -> dict:
    """Prepare a single receptor for MD."""
    result = {'name': name, 'success': False}
    output_dir = pdb_path.parent
    
    print(f"\n{'='*60}")
    print(f" {name}")
    print(f"{'='*60}")
    
    try:
        # Load
        start = time.time()
        pdb = app.PDBFile(str(pdb_path))
        print(f"  Loaded: {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues")
        
        # Force field
        ff = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Add hydrogens
        print(f"  Adding hydrogens...")
        modeller.addHydrogens(ff, pH=7.4)
        print(f"  After H: {modeller.topology.getNumAtoms()} atoms")
        
        # Save protein
        with open(output_dir / f'{name}_prepared.pdb', 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
        # Solvate with smaller padding for speed
        print(f"  Adding solvent...")
        modeller.addSolvent(ff, model='tip3p', padding=0.8*unit.nanometers, neutralize=True)
        print(f"  Solvated: {modeller.topology.getNumAtoms()} atoms")
        
        with open(output_dir / f'{name}_solvated.pdb', 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
        # Create system
        print(f"  Creating system...")
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds
        )
        
        integrator = mm.LangevinMiddleIntegrator(
            310*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds
        )
        
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            platform = mm.Platform.getPlatformByName('CPU')
        
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        
        # Minimize
        state = simulation.context.getState(getEnergy=True)
        init_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Initial E: {init_E:.1f} kJ/mol")
        
        print(f"  Minimizing...")
        simulation.minimizeEnergy(maxIterations=1000)
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        min_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Minimized E: {min_E:.1f} kJ/mol")
        
        with open(output_dir / 'minimized.pdb', 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Equilibrate
        print(f"  Equilibrating (10 ps)...")
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
        simulation.step(5000)
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        eq_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Equilibrated E: {eq_E:.1f} kJ/mol")
        
        with open(output_dir / 'equilibrated.pdb', 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        elapsed = time.time() - start
        
        result['success'] = True
        result['atoms'] = simulation.topology.getNumAtoms()
        result['residues'] = pdb.topology.getNumResidues()
        result['initial_energy'] = init_E
        result['minimized_energy'] = min_E
        result['equilibrated_energy'] = eq_E
        result['time'] = elapsed
        result['output_dir'] = str(output_dir)
        
        print(f"  ✓ {name} done in {elapsed:.1f}s")
        
    except Exception as e:
        result['error'] = str(e)
        print(f"  ✗ Error: {e}")
    
    return result


def main():
    base_dir = Path(__file__).parent
    
    print("=" * 70)
    print(" BATCH RECEPTOR PREPARATION FOR MD SIMULATIONS")
    print("=" * 70)
    
    results = {}
    successful = []
    
    for name, rel_path in RECEPTORS.items():
        pdb_path = base_dir / rel_path
        if pdb_path.exists():
            result = prepare_receptor(name, pdb_path)
            results[name] = result
            if result['success']:
                successful.append(name)
        else:
            print(f"\n⚠️  {name}: File not found: {pdb_path}")
            results[name] = {'name': name, 'success': False, 'error': 'File not found'}
    
    # Summary
    print("\n" + "=" * 70)
    print(" FINAL SUMMARY")
    print("=" * 70)
    
    for name, result in results.items():
        if result['success']:
            print(f"  ✓ {name}: {result['residues']} res, {result['atoms']} atoms, "
                  f"E={result['equilibrated_energy']:.0f} kJ/mol ({result['time']:.1f}s)")
        else:
            print(f"  ✗ {name}: {result.get('error', 'Unknown error')[:50]}")
    
    print(f"\n  Successfully prepared: {len(successful)}/{len(RECEPTORS)}")
    
    # Save results
    with open(base_dir / 'batch_preparation_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    return results


if __name__ == '__main__':
    main()
