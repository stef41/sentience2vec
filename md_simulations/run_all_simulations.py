#!/usr/bin/env python3
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

SYSTEMS = [
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_LSD",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_psilocin",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_DMT",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_mescaline",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2A_MDMA",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/5-HT2B_MDMA",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/NMDA_ketamine",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/KOR_salvinorin_A",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/CB1_THC",
    "/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD"
]

def run_phase(system_dir: str, phase: str) -> bool:
    """Run a simulation phase for a system"""
    scripts = {
        "minimize": "01_minimize.py",
        "equilibrate": "02_equilibrate.py",
        "production": "03_production.py",
        "analyze": "04_analyze.py"
    }
    
    script = Path(system_dir) / scripts[phase]
    if not script.exists():
        print(f"  ✗ Script not found: {script}")
        return False
    
    print(f"\n>>> Running {phase} for {Path(system_dir).name}...")
    try:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=system_dir,
            capture_output=True,
            text=True,
            timeout=3600 * 24  # 24 hour timeout
        )
        
        if result.returncode == 0:
            print(f"  ✓ {phase} complete")
            return True
        else:
            print(f"  ✗ {phase} failed:")
            print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ✗ {phase} timed out")
        return False
    except Exception as e:
        print(f"  ✗ {phase} error: {e}")
        return False


def main():
    print("="*70)
    print("MOLECULAR DYNAMICS SIMULATION RUNNER")
    print("="*70)
    print(f"\nSystems to simulate: {len(SYSTEMS)}")
    
    results = {}
    
    for system_dir in SYSTEMS:
        system_name = Path(system_dir).name
        print(f"\n" + "="*70)
        print(f"SYSTEM: {system_name}")
        print("="*70)
        
        results[system_name] = {}
        
        # Run each phase
        for phase in ["minimize", "equilibrate", "production", "analyze"]:
            success = run_phase(system_dir, phase)
            results[system_name][phase] = "success" if success else "failed"
            
            if not success and phase in ["minimize", "equilibrate"]:
                print(f"  Skipping remaining phases for {system_name}")
                break
    
    # Summary
    print("\n" + "="*70)
    print("SIMULATION SUMMARY")
    print("="*70)
    
    for system, phases in results.items():
        status = "✓" if all(v == "success" for v in phases.values()) else "✗"
        print(f"  {status} {system}: {phases}")
    
    # Save results
    with open("/data/users/zacharie/whc/consciousness_study/md_simulations/simulation_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: /data/users/zacharie/whc/consciousness_study/md_simulations/simulation_results.json")


if __name__ == "__main__":
    main()
