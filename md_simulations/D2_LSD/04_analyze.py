#!/usr/bin/env python3
"""
Trajectory Analysis for D2_LSD
Analyzes binding interactions and structural stability
"""
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, distances
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json

print("="*60)
print("TRAJECTORY ANALYSIS: D2_LSD")
print("="*60)

# Load trajectory
print("\nLoading trajectory...")
try:
    u = mda.Universe(
        "/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/minimized.pdb",  # Topology
        "/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/trajectory.dcd"   # Trajectory
    )
    print(f"  Frames: {len(u.trajectory)}")
    print(f"  Atoms: {u.atoms.n_atoms}")
except Exception as e:
    print(f"Error loading trajectory: {e}")
    print("Using final structure for static analysis...")
    u = mda.Universe("/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/final.pdb")

# Define selections
protein = u.select_atoms("protein")
backbone = u.select_atoms("protein and backbone")
binding_site = u.select_atoms("protein and resid 114 115 118 119 194 197 379 380 383 386 389")

print(f"\nSelections:")
print(f"  Protein atoms: {protein.n_atoms}")
print(f"  Backbone atoms: {backbone.n_atoms}")
print(f"  Binding site atoms: {binding_site.n_atoms}")
print(f"  Binding residues: [114, 115, 118, 119, 194, 197, 379, 380, 383, 386, 389]")

# ========================================
# RMSD Analysis
# ========================================
print("\nCalculating RMSD...")

if hasattr(u.trajectory, '__len__') and len(u.trajectory) > 1:
    # Align trajectory to first frame
    aligner = align.AlignTraj(u, u, select='backbone', in_memory=True)
    aligner.run()
    
    # Calculate RMSD
    rmsd_data = rms.RMSD(u, u, select='backbone', ref_frame=0)
    rmsd_data.run()
    
    rmsd_values = rmsd_data.results.rmsd[:, 2]  # RMSD column
    time_ps = rmsd_data.results.rmsd[:, 1]      # Time column
    
    print(f"  Average RMSD: {rmsd_values.mean():.2f} Å")
    print(f"  Max RMSD: {rmsd_values.max():.2f} Å")
    print(f"  Final RMSD: {rmsd_values[-1]:.2f} Å")
    
    # ========================================
    # Binding Site RMSD
    # ========================================
    print("\nCalculating binding site RMSD...")
    bs_rmsd = rms.RMSD(u, u, select='resid 114 115 118 119 194 197 379 380 383 386 389 and backbone', ref_frame=0)
    bs_rmsd.run()
    bs_rmsd_values = bs_rmsd.results.rmsd[:, 2]
    print(f"  Average: {bs_rmsd_values.mean():.2f} Å")
    
    # ========================================
    # Plotting
    # ========================================
    print("\nGenerating plots...")
    
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
        log_data = pd.read_csv("/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/production.log", sep=',', skiprows=1, 
                               names=['step', 'time', 'PE', 'KE', 'TE', 'T', 'density', 'speed'])
        ax4.plot(log_data['time']/1000, log_data['TE'])
        ax4.set_xlabel("Time (ns)")
        ax4.set_ylabel("Total Energy (kJ/mol)")
        ax4.set_title("Total Energy")
    except:
        ax4.text(0.5, 0.5, "Energy data not available", ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig("/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/analysis_results.png", dpi=150)
    print(f"  Saved: /data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/analysis_results.png")
    
    # ========================================
    # Save Results
    # ========================================
    results = {
        "system": "D2_LSD",
        "n_frames": len(u.trajectory),
        "backbone_rmsd_mean": float(rmsd_values.mean()),
        "backbone_rmsd_std": float(rmsd_values.std()),
        "backbone_rmsd_max": float(rmsd_values.max()),
        "binding_site_rmsd_mean": float(bs_rmsd_values.mean()),
        "binding_site_rmsd_std": float(bs_rmsd_values.std()),
        "binding_residues": [114, 115, 118, 119, 194, 197, 379, 380, 383, 386, 389]
    }
    
else:
    print("  Single frame - static analysis only")
    results = {
        "system": "D2_LSD",
        "n_frames": 1,
        "binding_residues": [114, 115, 118, 119, 194, 197, 379, 380, 383, 386, 389]
    }

with open("/data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/analysis_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n✓ Analysis complete!")
print(f"  Results: /data/users/zacharie/whc/consciousness_study/md_simulations/D2_LSD/analysis_results.json")
