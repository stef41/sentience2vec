#!/usr/bin/env python3
"""
Full Docking Pipeline for Novel Molecules

Integrates:
1. SMILES → 3D structure generation (RDKit)
2. Chai-1 pose prediction
3. OpenMM energy minimization
4. Binding energy scoring

With caching at every step to avoid redundant computation.
"""

import os
import sys
import json
import hashlib
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Try to import dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

OPENMM_AVAILABLE = False
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    OPENMM_AVAILABLE = True
except ImportError:
    pass

CHAI_AVAILABLE = False
try:
    from chai_lab.chai1 import run_inference
    CHAI_AVAILABLE = True
except ImportError:
    pass


# Base directories
BASE_DIR = Path(__file__).parent.parent / "md_simulations"
CACHE_DIR = Path(__file__).parent / "cache"
RECEPTORS_DIR = BASE_DIR / "receptors_full"


@dataclass
class DockingResult:
    """Result of docking a molecule to a receptor."""
    smiles: str
    receptor: str
    binding_energy_kJ: float = 0.0
    initial_energy_kJ: float = 0.0
    minimized_energy_kJ: float = 0.0
    pose_confidence: float = 0.0
    success: bool = False
    error: Optional[str] = None
    cached: bool = False
    timestamp: str = ""
    
    def to_dict(self) -> Dict:
        return asdict(self)


class DockingCache:
    """
    Persistent cache for docking results.
    
    Uses SMILES hash + receptor name as key.
    Stores results in JSON files for persistence.
    """
    
    def __init__(self, cache_dir: Path = None):
        self.cache_dir = cache_dir or CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.memory_cache: Dict[str, DockingResult] = {}
        self._load_cache_index()
    
    def _load_cache_index(self):
        """Load cache index from disk."""
        index_file = self.cache_dir / "cache_index.json"
        if index_file.exists():
            with open(index_file, 'r') as f:
                self.index = json.load(f)
        else:
            self.index = {}
    
    def _save_cache_index(self):
        """Save cache index to disk."""
        index_file = self.cache_dir / "cache_index.json"
        with open(index_file, 'w') as f:
            json.dump(self.index, f, indent=2)
    
    def _get_key(self, smiles: str, receptor: str) -> str:
        """Generate cache key from SMILES and receptor."""
        # Canonicalize SMILES if possible
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Hash for filename safety
        content = f"{smiles}:{receptor}"
        return hashlib.md5(content.encode()).hexdigest()
    
    def get(self, smiles: str, receptor: str) -> Optional[DockingResult]:
        """Get cached result if available."""
        key = self._get_key(smiles, receptor)
        
        # Check memory cache first
        if key in self.memory_cache:
            result = self.memory_cache[key]
            result.cached = True
            return result
        
        # Check disk cache
        cache_file = self.cache_dir / f"{key}.json"
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    data = json.load(f)
                result = DockingResult(**data)
                result.cached = True
                self.memory_cache[key] = result
                return result
            except Exception:
                pass
        
        return None
    
    def put(self, result: DockingResult):
        """Store result in cache."""
        key = self._get_key(result.smiles, result.receptor)
        
        # Store in memory
        self.memory_cache[key] = result
        
        # Store on disk
        cache_file = self.cache_dir / f"{key}.json"
        result.timestamp = datetime.now().isoformat()
        with open(cache_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        
        # Update index
        self.index[key] = {
            'smiles': result.smiles,
            'receptor': result.receptor,
            'energy': result.binding_energy_kJ,
            'success': result.success,
            'timestamp': result.timestamp
        }
        self._save_cache_index()
    
    def get_stats(self) -> Dict:
        """Get cache statistics."""
        return {
            'total_entries': len(self.index),
            'memory_cached': len(self.memory_cache),
            'cache_dir': str(self.cache_dir)
        }
    
    def clear(self):
        """Clear all cache."""
        self.memory_cache.clear()
        self.index.clear()
        for f in self.cache_dir.glob("*.json"):
            f.unlink()


def smiles_to_3d(smiles: str, output_path: Path) -> bool:
    """Convert SMILES to 3D SDF file."""
    if not RDKIT_AVAILABLE:
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            # Try with random coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        
        # Write SDF
        output_path.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        
        return True
    except Exception as e:
        print(f"    Error generating 3D: {e}")
        return False


def extract_sequence_from_pdb(pdb_path: Path) -> Optional[str]:
    """Extract amino acid sequence from PDB file."""
    AA_MAP = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'CYX': 'C'
    }
    
    residues = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    chain = line[21]
                    if chain == 'A' and res_name in AA_MAP:
                        residues[res_num] = AA_MAP[res_name]
    except Exception:
        return None
    
    if not residues:
        return None
    
    return ''.join(residues[k] for k in sorted(residues.keys()))


def run_chai1_docking(receptor_pdb: Path, ligand_smiles: str, 
                      output_dir: Path, gpu_id: int = 0) -> Tuple[Optional[Path], float]:
    """
    Run Chai-1 pose prediction for a ligand.
    
    Returns: (pose_path, confidence_score)
    """
    if not CHAI_AVAILABLE:
        return None, 0.0
    
    # Get receptor sequence
    sequence = extract_sequence_from_pdb(receptor_pdb)
    if not sequence:
        return None, 0.0
    
    # Prepare FASTA input
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    fasta_content = f">protein|name=receptor\n{sequence}\n>ligand|name=ligand\n{ligand_smiles}\n"
    
    fasta_path = output_dir.parent / f"{output_dir.name}_input.fasta"
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    # Clear and create output directory
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Run Chai-1
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=output_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device=f"cuda:{gpu_id}",
            use_esm_embeddings=True,
        )
        
        if candidates and len(candidates.cif_paths) > 0:
            # Get confidence from aggregate scores
            confidence = 0.5  # Default
            if hasattr(candidates, 'aggregate_scores'):
                scores = candidates.aggregate_scores
                if len(scores) > 0:
                    confidence = float(scores[0].get('ptm', 0.5))
            
            return candidates.cif_paths[0], confidence
    
    except Exception as e:
        print(f"    Chai-1 error: {str(e)[:80]}")
    
    return None, 0.0


def run_openmm_minimization(cif_path: Path, output_dir: Path, 
                            gpu_id: int = 0) -> Tuple[float, float]:
    """
    Run OpenMM energy minimization on predicted pose.
    
    Returns: (initial_energy, final_energy) in kJ/mol
    """
    if not OPENMM_AVAILABLE:
        return 0.0, 0.0
    
    try:
        from pdbfixer import PDBFixer
        
        # Load and fix structure
        fixer = PDBFixer(str(cif_path))
        
        # Keep only standard amino acids
        aa_names = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
                    'TYR', 'VAL', 'HIE', 'HID', 'HIP', 'CYX'}
        
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        
        if fixer.topology.getNumAtoms() < 100:
            return 0.0, 0.0
        
        # Set up force field with implicit solvent
        forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds
        )
        
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picoseconds,
            0.002 * unit.picoseconds
        )
        
        # Use GPU
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'DeviceIndex': str(gpu_id)}
        
        simulation = app.Simulation(
            fixer.topology, system, integrator, platform, properties
        )
        simulation.context.setPositions(fixer.positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Get final energy
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Save minimized structure
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "minimized.pdb", 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        return initial_energy, final_energy
    
    except Exception as e:
        print(f"    OpenMM error: {str(e)[:80]}")
        return 0.0, 0.0


class DockingPipeline:
    """
    Full docking pipeline with caching.
    
    Workflow:
    1. Check cache for existing results
    2. Generate 3D structure from SMILES
    3. Run Chai-1 pose prediction
    4. Run OpenMM energy minimization
    5. Store results in cache
    """
    
    def __init__(self, cache_dir: Path = None, gpu_id: int = 0):
        self.cache = DockingCache(cache_dir)
        self.gpu_id = gpu_id
        self.work_dir = Path(tempfile.mkdtemp(prefix="docking_"))
        
        # Load available receptors
        self.receptors = self._find_receptors()
    
    def _find_receptors(self) -> Dict[str, Path]:
        """Find all prepared receptor PDB files."""
        receptors = {}
        
        if RECEPTORS_DIR.exists():
            for receptor_dir in RECEPTORS_DIR.iterdir():
                if receptor_dir.is_dir():
                    pdb_file = receptor_dir / f"{receptor_dir.name}_prepared.pdb"
                    if pdb_file.exists():
                        receptors[receptor_dir.name] = pdb_file
        
        return receptors
    
    def get_available_receptors(self) -> List[str]:
        """Get list of available receptor names."""
        return list(self.receptors.keys())
    
    def dock_molecule(self, smiles: str, receptor: str, 
                      use_cache: bool = True) -> DockingResult:
        """
        Dock a molecule to a receptor.
        
        Args:
            smiles: SMILES string of the molecule
            receptor: Name of the receptor (e.g., "5-HT2A")
            use_cache: Whether to use cached results
        
        Returns:
            DockingResult with binding energy and metadata
        """
        # Check cache first
        if use_cache:
            cached = self.cache.get(smiles, receptor)
            if cached is not None:
                return cached
        
        result = DockingResult(
            smiles=smiles,
            receptor=receptor,
            timestamp=datetime.now().isoformat()
        )
        
        # Find receptor PDB
        receptor_pdb = self.receptors.get(receptor)
        if receptor_pdb is None:
            result.error = f"Receptor {receptor} not found"
            return result
        
        # Create working directory for this docking
        mol_hash = hashlib.md5(smiles.encode()).hexdigest()[:8]
        work_subdir = self.work_dir / f"{receptor}_{mol_hash}"
        work_subdir.mkdir(parents=True, exist_ok=True)
        
        # Step 1: Generate 3D structure
        ligand_sdf = work_subdir / "ligand.sdf"
        if not smiles_to_3d(smiles, ligand_sdf):
            result.error = "Failed to generate 3D structure"
            self.cache.put(result)
            return result
        
        # Step 2: Run Chai-1 pose prediction
        chai_dir = work_subdir / "chai1"
        pose_path, confidence = run_chai1_docking(
            receptor_pdb, smiles, chai_dir, self.gpu_id
        )
        
        result.pose_confidence = confidence
        
        if pose_path is None or not pose_path.exists():
            # Fall back to estimated energy based on structure
            result.binding_energy_kJ = self._estimate_binding_energy(smiles)
            result.error = "Chai-1 prediction failed, using estimate"
            self.cache.put(result)
            return result
        
        # Step 3: Run OpenMM minimization
        md_dir = work_subdir / "md"
        initial_e, final_e = run_openmm_minimization(pose_path, md_dir, self.gpu_id)
        
        result.initial_energy_kJ = initial_e
        result.minimized_energy_kJ = final_e
        result.binding_energy_kJ = final_e  # Use minimized energy as binding proxy
        result.success = True
        
        # Cache result
        self.cache.put(result)
        
        return result
    
    def _estimate_binding_energy(self, smiles: str) -> float:
        """
        Estimate binding energy from molecular properties.
        Used as fallback when full pipeline fails.
        """
        if not RDKIT_AVAILABLE:
            return -50000.0  # Default estimate
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return -50000.0
            
            # Simple estimate based on molecular properties
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            
            # Rough estimate (more negative = stronger binding)
            # Based on typical drug-like properties
            energy = -40000.0  # Base
            energy -= mw * 50  # Size contribution
            energy -= logp * 2000  # Lipophilicity
            energy -= (hbd + hba) * 1000  # H-bonding
            energy += tpsa * 100  # Polar area penalty
            
            return energy
        except Exception:
            return -50000.0
    
    def dock_to_all_receptors(self, smiles: str, 
                              receptors: List[str] = None) -> Dict[str, DockingResult]:
        """
        Dock a molecule to multiple receptors.
        
        Returns dict of receptor -> DockingResult
        """
        if receptors is None:
            receptors = list(self.receptors.keys())
        
        results = {}
        for receptor in receptors:
            print(f"  Docking to {receptor}...")
            results[receptor] = self.dock_molecule(smiles, receptor)
        
        return results
    
    def get_binding_profile(self, smiles: str, 
                           receptors: List[str] = None) -> Dict[str, float]:
        """
        Get binding energy profile for a molecule.
        
        Returns dict of receptor -> binding energy (kJ/mol)
        """
        dock_results = self.dock_to_all_receptors(smiles, receptors)
        
        return {
            receptor: result.binding_energy_kJ 
            for receptor, result in dock_results.items()
        }
    
    def cleanup(self):
        """Clean up temporary files."""
        if self.work_dir.exists():
            shutil.rmtree(self.work_dir)
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.cleanup()


def main():
    """Test the docking pipeline."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run docking pipeline")
    parser.add_argument("smiles", help="SMILES string to dock")
    parser.add_argument("--receptor", "-r", default="5-HT2A",
                       help="Receptor to dock to")
    parser.add_argument("--all", "-a", action="store_true",
                       help="Dock to all available receptors")
    parser.add_argument("--gpu", "-g", type=int, default=0,
                       help="GPU ID to use")
    parser.add_argument("--no-cache", action="store_true",
                       help="Don't use cache")
    
    args = parser.parse_args()
    
    print(f"Docking Pipeline")
    print(f"================")
    print(f"SMILES: {args.smiles}")
    print(f"RDKit available: {RDKIT_AVAILABLE}")
    print(f"OpenMM available: {OPENMM_AVAILABLE}")
    print(f"Chai-1 available: {CHAI_AVAILABLE}")
    
    with DockingPipeline(gpu_id=args.gpu) as pipeline:
        print(f"\nAvailable receptors: {len(pipeline.receptors)}")
        
        cache_stats = pipeline.cache.get_stats()
        print(f"Cache entries: {cache_stats['total_entries']}")
        
        if args.all:
            print(f"\nDocking to all receptors...")
            results = pipeline.dock_to_all_receptors(
                args.smiles, 
            )
            
            print(f"\n{'Receptor':<20} {'Energy (kJ/mol)':<20} {'Status'}")
            print("-" * 60)
            for receptor, result in sorted(results.items(), key=lambda x: x[1].binding_energy_kJ):
                status = "✓ cached" if result.cached else ("✓" if result.success else f"✗ {result.error}")
                print(f"{receptor:<20} {result.binding_energy_kJ:<20.0f} {status}")
        else:
            print(f"\nDocking to {args.receptor}...")
            result = pipeline.dock_molecule(
                args.smiles, 
                args.receptor, 
                use_cache=not args.no_cache
            )
            
            print(f"\nResult:")
            print(f"  Success: {result.success}")
            print(f"  Cached: {result.cached}")
            print(f"  Binding Energy: {result.binding_energy_kJ:.0f} kJ/mol")
            print(f"  Initial Energy: {result.initial_energy_kJ:.0f} kJ/mol")
            print(f"  Minimized Energy: {result.minimized_energy_kJ:.0f} kJ/mol")
            print(f"  Pose Confidence: {result.pose_confidence:.3f}")
            if result.error:
                print(f"  Error: {result.error}")


if __name__ == "__main__":
    main()
