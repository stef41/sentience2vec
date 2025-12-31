# Molecule Generator from Docking Profile

Generate novel molecules that match a target receptor binding profile.

## Overview

Given a desired binding profile (e.g., "strong 5-HT2A agonist, weak D2 binding"), this tool generates valid molecular structures that are predicted to have similar properties.

## Approach

1. **Profile Matching**: Find existing molecules with similar binding profiles
2. **Fragment Analysis**: Extract common structural fragments from matched molecules
3. **Molecule Generation**: Combine fragments and apply mutations to create novel structures
4. **Validation**: Ensure generated molecules are chemically valid (RDKit)
5. **Scoring**: Rank by predicted profile similarity and drug-likeness

## Usage

```bash
# Generate molecules for a target profile
python generate_molecule.py --profile "5-HT2A:-90000,D2:-50000,5-HT2B:-70000"

# Generate molecules similar to a reference drug
python generate_molecule.py --like "LSD" --modify "increase D2"

# Interactive mode
python generate_molecule.py --interactive
```

## Requirements

```bash
pip install rdkit numpy scipy
```

## Output

- SMILES strings for generated molecules
- Predicted binding profiles
- Drug-likeness scores (Lipinski, QED)
- 2D structure images
