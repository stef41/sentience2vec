# Consciousness Study Framework

A multi-modal research framework for mapping subjective experiences to brain regions using:
- **Phenomenological data** from Erowid, PsychonautWiki, and TripSit
- **Local LLM** for NLP/effect extraction
- **Molecular dynamics** for receptor binding analysis

## Data Sources

| Source | API/Access | Data Type |
|--------|------------|-----------|
| PsychonautWiki | GraphQL API (`api.psychonautwiki.org`) | Substances, effects, dosages, pharmacology |
| TripSit | REST API (`tripbot.tripsit.me/api/tripsit/getDrug`) | Factsheets, combinations, dosing |
| Erowid | Web scraping (respect ToS) | ~40k+ experience reports |

## Project Structure

```
consciousness_study/
├── data/
│   ├── raw/                    # Raw scraped/API data
│   ├── processed/              # Cleaned datasets
│   └── knowledge_graph/        # Graph database exports
├── src/
│   ├── data_collection/        # API clients & scrapers
│   ├── nlp/                    # LLM-based effect extraction
│   ├── pharmacology/           # Receptor binding analysis
│   ├── molecular_dynamics/     # MD simulation scripts
│   └── analysis/               # Statistical modeling
├── notebooks/                  # Jupyter notebooks for exploration
└── configs/                    # Configuration files
```

## Key Mappings

### Subjective Effect Index (PsychonautWiki)
- **Visual**: geometry, drifting, color enhancement, hallucinations
- **Cognitive**: ego dissolution, time distortion, thought loops
- **Physical**: stimulation, sedation, body high
- **Transpersonal**: unity, interconnectedness, mystical states

### Receptor Systems → Brain Regions → Experiences
| Receptor | Primary Regions | Associated Qualia |
|----------|-----------------|-------------------|
| 5-HT2A | V1, PFC, claustrum | Visual geometry, ego dissolution |
| 5-HT1A | Raphe, hippocampus | Anxiolysis, introspection |
| D2 | Striatum, NAcc | Euphoria, motivation |
| NMDA | Hippocampus, cortex | Dissociation, k-hole |
| κ-opioid | Insula, limbic | Dysphoria, depersonalization |
| CB1 | Hippocampus, cerebellum | Time dilation, memory effects |

## Installation

```bash
pip install -r requirements.txt
```

## Usage

1. **Collect data**: `python -m src.data_collection.main`
2. **Extract effects**: `python -m src.nlp.extract_effects`
3. **Build knowledge graph**: `python -m src.analysis.build_graph`
4. **Run MD simulations**: See `src/molecular_dynamics/README.md`
Receptors with detailed density maps:

