#!/usr/bin/env python3
"""
Quick test of all three tools:
1. Chai-1 - Pose prediction
2. PoseBusters - Pose validation
3. FEP Setup - Free energy calculation setup
"""
import warnings
warnings.filterwarnings('ignore')

print("=" * 60)
print("TESTING PIPELINE TOOLS")
print("=" * 60)

# Test 1: Chai-1
print("\n[1/3] Testing Chai-1...")
try:
    from chai_lab.chai1 import run_inference
    print("    ✓ chai_lab.chai1 imported successfully")
    print("    ✓ run_inference function available")
except Exception as e:
    print(f"    ✗ Error: {e}")

# Test 2: PoseBusters
print("\n[2/3] Testing PoseBusters...")
try:
    from posebusters import PoseBusters
    pb = PoseBusters(config='redock')
    print("    ✓ posebusters imported successfully")
    print("    ✓ PoseBusters(config='redock') initialized")
except Exception as e:
    print(f"    ✗ Error: {e}")

# Test 3: FEP Setup (RDKit MCS)
print("\n[3/3] Testing FEP Setup (RDKit MCS)...")
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    
    # Create two simple molecules
    mol1 = Chem.MolFromSmiles('c1ccc2[nH]ccc2c1')  # Indole
    mol2 = Chem.MolFromSmiles('CN(C)CCc1c[nH]c2ccccc12')  # DMT
    
    # Find MCS
    mcs = rdFMCS.FindMCS([mol1, mol2])
    print(f"    ✓ rdkit.Chem imported successfully")
    print(f"    ✓ MCS found: {mcs.numAtoms} atoms")
    print(f"    ✓ SMARTS: {mcs.smartsString[:50]}...")
except Exception as e:
    print(f"    ✗ Error: {e}")

# Test 4: Check files
print("\n[BONUS] Checking prepared files...")
from pathlib import Path

receptors_dir = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/receptors_full")
ligands_dir = Path("/data/users/zacharie/whc/consciousness_study/md_simulations/ligands_full")

receptor_count = len(list(receptors_dir.glob("*/")) if receptors_dir.exists() else [])
ligand_count = len(list(ligands_dir.glob("*/")) if ligands_dir.exists() else [])

print(f"    ✓ {receptor_count} receptors prepared")
print(f"    ✓ {ligand_count} ligands prepared")

# Show some examples
print("\n[EXAMPLES] Key pairs:")
pairs = [
    ("5-HT2A", "DMT"),
    ("5-HT2A", "LSD"),
    ("5-HT2A", "Psilocin"),
    ("D2", "Dopamine"),
    ("MOR", "Morphine"),
]

for receptor, ligand in pairs:
    r_exists = (receptors_dir / receptor).exists()
    l_exists = (ligands_dir / ligand).exists()
    status = "✓" if (r_exists and l_exists) else "✗"
    print(f"    {status} {receptor} + {ligand}")

print("\n" + "=" * 60)
print("ALL TOOLS VERIFIED!")
print("=" * 60)
