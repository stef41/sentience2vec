# Full Consciousness Study Structure Preparation Summary

## Overview

Comprehensive preparation of ALL neuroreceptor structures and ALL psychoactive 
compounds from PsychonautWiki for molecular dynamics simulations.

## Ligand Structures

| Category | Count |
|----------|-------|
| Total compounds in database | 320 |
| **Compounds with 3D structures** | **301** |
| Failed (generic categories/plants) | 19 |

### Ligand Sources
- **PubChem 3D**: Direct download of conformer SDF files
- **RDKit**: SMILES → 3D conversion for compounds not in PubChem

### Classes Covered
- Psychedelics (LSD, psilocin, DMT, mescaline, 2C-x, NBOMe series)
- Dissociatives (ketamine, PCP, MXE, DXM)
- Stimulants (amphetamines, cathinones, cocaine)
- Depressants (benzodiazepines, barbiturates, opioids)
- Cannabinoids (THC, CBD, synthetic cannabinoids)
- Nootropics and supplements
- Entactogens (MDMA, MDA)
- Deliriants (scopolamine, diphenhydramine)

## Receptor Structures

| Category | Count |
|----------|-------|
| Total receptor subtypes | 84 |
| **Receptors with PDB structures** | **76** |
| No crystal structure available | 8 |

### Receptor Families
| Family | Subtypes | With PDB |
|--------|----------|----------|
| Serotonin (5-HT) | 14 | 12 |
| Glutamate (NMDA/AMPA/mGluR) | 14 | 12 |
| Adrenergic (α/β) | 9 | 7 |
| Dopamine (D1-D5) | 5 | 5 |
| Muscarinic (M1-M5) | 5 | 5 |
| GABA (A/B) | 5 | 5 |
| Opioid (μ/κ/δ/NOP) | 4 | 4 |
| Histamine (H1-H4) | 4 | 4 |
| Nicotinic (nAChR) | 4 | 4 |
| Adenosine (A1-A3) | 4 | 4 |
| Cannabinoid (CB1/CB2) | 2 | 2 |
| Sigma (σ1/σ2) | 2 | 2 |
| Melatonin (MT1/MT2) | 2 | 2 |
| Orexin (OX1/OX2) | 2 | 2 |
| TRP Channels (TRPV1/TRPM8) | 2 | 2 |
| Transporters (SERT/DAT/NET/VMAT2) | 4 | 4 |
| Others (TAAR1, GlyR) | 2 | 2 |

## Directory Structure

\`\`\`
md_simulations/
├── ligands_full/           # 301 compounds with 3D structures
│   ├── LSD/
│   │   ├── LSD.sdf
│   │   └── LSD.pdb
│   ├── Ketamine/
│   ├── MDMA/
│   └── ...
├── receptors_full/         # 76 receptor structures
│   ├── 5-HT1A/
│   │   └── 5-HT1A_raw.pdb
│   ├── 5-HT2A/
│   ├── D2/
│   └── ...
├── ligands/                # Original 91 curated ligands
├── prepared_receptors/     # Original 20 equilibrated receptors
└── receptor_pdb_database_full.py  # Full receptor database
\`\`\`

## Total Possible Combinations

- **301 ligands × 76 receptors = 22,876 potential receptor-ligand pairs**

This provides comprehensive coverage for studying:
- Classical psychedelic binding (5-HT2A, 5-HT1A)
- Dissociative mechanisms (NMDA, σ receptors)
- Opioid system interactions
- Cannabinoid signaling
- Dopamine/norepinephrine/serotonin transporter inhibition
- GABA-ergic effects
- Novel compound predictions

## Files Generated

- \`full_preparation_results.json\` - Ligand preparation status
- \`receptor_preparation_results.json\` - Receptor preparation status
- \`receptor_pdb_database_full.py\` - Complete receptor → PDB mapping
- \`prepare_all_full.py\` - Full ligand preparation script
- \`prepare_all_receptors.py\` - Full receptor preparation script

## Next Steps for MD Simulations

1. **Docking**: Use AutoDock Vina for all pairs
2. **MD Setup**: Solvate receptor-ligand complexes with OpenMM
3. **Production**: Run ns-scale simulations on GPU
4. **Analysis**: Calculate binding free energies, interaction fingerprints

---
Generated: $(date)
