"""
Receptor Binding Data and Brain Region Mapping
Links pharmacology to neuroanatomy to phenomenology
"""

import json
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict, field
from enum import Enum
import pandas as pd
import numpy as np


class ReceptorActivity(Enum):
    """Type of receptor interaction"""
    FULL_AGONIST = "full_agonist"
    PARTIAL_AGONIST = "partial_agonist"
    AGONIST = "agonist"  # Generic agonist (when efficacy unknown)
    ANTAGONIST = "antagonist"
    INVERSE_AGONIST = "inverse_agonist"
    ALLOSTERIC_MODULATOR = "allosteric_modulator"
    REUPTAKE_INHIBITOR = "reuptake_inhibitor"
    RELEASING_AGENT = "releasing_agent"


@dataclass
class ReceptorBinding:
    """Receptor binding profile for a substance"""
    receptor: str              # e.g., "5-HT2A", "D2", "NMDA"
    ki_nm: Optional[float]     # Binding affinity (nM) - lower = stronger
    activity: ReceptorActivity
    ec50_nm: Optional[float] = None  # Functional potency
    efficacy: Optional[float] = None  # 0-1, relative to reference agonist
    source: Optional[str] = None


@dataclass
class SubstancePharmacology:
    """Complete pharmacological profile"""
    name: str
    bindings: List[ReceptorBinding]
    primary_targets: List[str]  # Main receptors responsible for effects
    half_life_hours: Optional[float] = None
    bioavailability: Optional[Dict[str, float]] = None  # By route
    metabolism: Optional[str] = None
    active_metabolites: List[str] = field(default_factory=list)


# Major neurotransmitter receptors and their brain distributions
# Based on autoradiography, PET imaging, and immunohistochemistry studies
# Density values in fmol/mg protein where available (from human post-mortem and PET studies)

RECEPTOR_BRAIN_REGIONS = {
    # ==========================================================================
    # SEROTONIN (5-HT) RECEPTORS
    # Source: Pazos & Palacios 1985, Hall et al. 2000, Beliveau et al. 2017
    # ==========================================================================
    "5-HT1A": {
        "density_map": {
            # Highest density regions (>100 fmol/mg)
            "dorsal_raphe_nucleus": {"density_fmol_mg": 450, "cell_type": "somatodendritic_autoreceptor"},
            "median_raphe_nucleus": {"density_fmol_mg": 380, "cell_type": "somatodendritic_autoreceptor"},
            "hippocampus_CA1": {"density_fmol_mg": 180, "cell_type": "pyramidal_dendrites"},
            "hippocampus_dentate_gyrus": {"density_fmol_mg": 165, "cell_type": "granule_cells"},
            "entorhinal_cortex_layer_II": {"density_fmol_mg": 145, "cell_type": "stellate_cells"},
            "septum_lateral": {"density_fmol_mg": 130, "cell_type": "GABAergic_neurons"},
            "amygdala_basolateral": {"density_fmol_mg": 120, "cell_type": "principal_neurons"},
            # Moderate density (40-100 fmol/mg)
            "prefrontal_cortex_layer_II": {"density_fmol_mg": 85, "cell_type": "pyramidal_apical_dendrites"},
            "anterior_cingulate_cortex": {"density_fmol_mg": 75, "cell_type": "layer_II_III_pyramidal"},
            "temporal_cortex": {"density_fmol_mg": 65, "cell_type": "pyramidal_neurons"},
            "hypothalamus_paraventricular": {"density_fmol_mg": 55, "cell_type": "neuroendocrine"},
            # Low density (<40 fmol/mg)
            "striatum": {"density_fmol_mg": 15, "cell_type": "sparse"},
            "cerebellum": {"density_fmol_mg": 8, "cell_type": "minimal"},
            "thalamus": {"density_fmol_mg": 25, "cell_type": "relay_neurons"},
        },
        "pet_tracer": "[11C]WAY-100635",
        "function": "Gi/o_coupled_inhibitory",
        "subcellular": "somatodendritic",
        "associated_effects": ["anxiolysis", "introspection", "mood_elevation", "sexual_dysfunction"]
    },
    
    "5-HT2A": {
        "density_map": {
            # Highest density - Layer V pyramidal neurons (>30 fmol/mg)
            "prefrontal_cortex_layer_V": {"density_fmol_mg": 48, "cell_type": "pyramidal_apical_dendrites"},
            "anterior_cingulate_BA24": {"density_fmol_mg": 45, "cell_type": "layer_V_pyramidal"},
            "posterior_cingulate_BA23": {"density_fmol_mg": 42, "cell_type": "layer_V_pyramidal"},
            "claustrum": {"density_fmol_mg": 55, "cell_type": "principal_neurons"},
            "visual_cortex_V1_layer_IV": {"density_fmol_mg": 38, "cell_type": "spiny_stellate"},
            "visual_cortex_V2": {"density_fmol_mg": 35, "cell_type": "layer_IV_V"},
            "temporal_cortex_BA21": {"density_fmol_mg": 32, "cell_type": "pyramidal"},
            "parietal_cortex_BA7": {"density_fmol_mg": 30, "cell_type": "layer_V"},
            "insula_anterior": {"density_fmol_mg": 36, "cell_type": "von_Economo_neurons"},
            # Moderate density (15-30 fmol/mg)
            "hippocampus_CA3": {"density_fmol_mg": 22, "cell_type": "pyramidal"},
            "amygdala_basolateral": {"density_fmol_mg": 25, "cell_type": "principal"},
            "entorhinal_cortex": {"density_fmol_mg": 28, "cell_type": "layer_II_stellate"},
            "nucleus_accumbens_shell": {"density_fmol_mg": 18, "cell_type": "medium_spiny"},
            "caudate_nucleus": {"density_fmol_mg": 15, "cell_type": "medium_spiny"},
            # Low density (<15 fmol/mg)
            "thalamus_pulvinar": {"density_fmol_mg": 12, "cell_type": "relay"},
            "cerebellum": {"density_fmol_mg": 3, "cell_type": "minimal"},
            "brainstem": {"density_fmol_mg": 8, "cell_type": "sparse"},
        },
        "pet_tracer": "[18F]altanserin, [11C]MDL100907",
        "function": "Gq_coupled_excitatory",
        "subcellular": "apical_dendrites_layer_V",
        "associated_effects": ["visual_geometry", "ego_dissolution", "mystical_experience", "synesthesia", "thought_disorder"]
    },
    
    "5-HT2B": {
        "density_map": {
            # Primarily peripheral - minimal CNS expression
            "heart_valves": {"density_fmol_mg": 250, "cell_type": "valvular_interstitial"},
            "stomach_fundus": {"density_fmol_mg": 180, "cell_type": "smooth_muscle"},
            "intestine": {"density_fmol_mg": 120, "cell_type": "myenteric_plexus"},
            "pulmonary_artery": {"density_fmol_mg": 90, "cell_type": "smooth_muscle"},
            # Limited CNS
            "cerebellum_purkinje": {"density_fmol_mg": 12, "cell_type": "purkinje_cells"},
            "dorsal_raphe": {"density_fmol_mg": 8, "cell_type": "serotonergic"},
            "hypothalamus": {"density_fmol_mg": 6, "cell_type": "sparse"},
        },
        "pet_tracer": "none_validated",
        "function": "Gq_coupled",
        "subcellular": "postsynaptic",
        "associated_effects": ["nausea", "cardiac_valvulopathy", "pulmonary_hypertension"]
    },
    
    "5-HT2C": {
        "density_map": {
            # Very high in choroid plexus
            "choroid_plexus": {"density_fmol_mg": 850, "cell_type": "epithelial"},
            # High limbic density
            "nucleus_accumbens": {"density_fmol_mg": 65, "cell_type": "medium_spiny_neurons"},
            "substantia_nigra_pars_compacta": {"density_fmol_mg": 55, "cell_type": "dopaminergic"},
            "VTA": {"density_fmol_mg": 50, "cell_type": "dopaminergic"},
            "amygdala_central": {"density_fmol_mg": 48, "cell_type": "GABAergic"},
            "hippocampus_CA3": {"density_fmol_mg": 42, "cell_type": "pyramidal"},
            "hypothalamus_arcuate": {"density_fmol_mg": 65, "cell_type": "POMC_neurons"},
            "prefrontal_cortex": {"density_fmol_mg": 28, "cell_type": "interneurons"},
            "habenula": {"density_fmol_mg": 45, "cell_type": "principal"},
            # Moderate
            "striatum": {"density_fmol_mg": 22, "cell_type": "cholinergic_interneurons"},
            "olfactory_tubercle": {"density_fmol_mg": 35, "cell_type": "medium_spiny"},
        },
        "pet_tracer": "[11C]IMiD",
        "function": "Gq_coupled_constitutively_active",
        "subcellular": "postsynaptic_GABAergic",
        "associated_effects": ["appetite_suppression", "anxiety", "dopamine_modulation", "CSF_production"]
    },
    
    # ==========================================================================
    # DOPAMINE RECEPTORS  
    # Source: Hall et al. 1994, Kessler et al. 1993, Seeman 2006
    # ==========================================================================
    "D1": {
        "density_map": {
            # Highest in striatum
            "caudate_nucleus": {"density_fmol_mg": 380, "cell_type": "medium_spiny_direct_pathway"},
            "putamen": {"density_fmol_mg": 420, "cell_type": "medium_spiny_direct_pathway"},
            "nucleus_accumbens_core": {"density_fmol_mg": 350, "cell_type": "medium_spiny"},
            "nucleus_accumbens_shell": {"density_fmol_mg": 280, "cell_type": "medium_spiny"},
            "olfactory_tubercle": {"density_fmol_mg": 260, "cell_type": "medium_spiny"},
            # Cortical (layer V-VI)
            "prefrontal_cortex_layer_V": {"density_fmol_mg": 45, "cell_type": "pyramidal"},
            "prefrontal_cortex_layer_VI": {"density_fmol_mg": 38, "cell_type": "pyramidal"},
            "anterior_cingulate": {"density_fmol_mg": 35, "cell_type": "pyramidal"},
            "entorhinal_cortex": {"density_fmol_mg": 28, "cell_type": "stellate"},
            # Limbic
            "amygdala_basolateral": {"density_fmol_mg": 32, "cell_type": "principal"},
            "hippocampus_subiculum": {"density_fmol_mg": 25, "cell_type": "pyramidal"},
            # Low
            "substantia_nigra": {"density_fmol_mg": 15, "cell_type": "pars_reticulata"},
            "thalamus": {"density_fmol_mg": 12, "cell_type": "sparse"},
            "cerebellum": {"density_fmol_mg": 2, "cell_type": "minimal"},
        },
        "pet_tracer": "[11C]SCH23390, [11C]NNC112",
        "function": "Gs_coupled_excitatory",
        "subcellular": "postsynaptic_dendritic_spines",
        "associated_effects": ["reward_learning", "working_memory", "motor_initiation", "motivation"]
    },
    
    "D2": {
        "density_map": {
            # Very high striatal density
            "caudate_nucleus": {"density_fmol_mg": 280, "cell_type": "medium_spiny_indirect_pathway"},
            "putamen": {"density_fmol_mg": 320, "cell_type": "medium_spiny_indirect_pathway"},
            "nucleus_accumbens": {"density_fmol_mg": 240, "cell_type": "medium_spiny"},
            "olfactory_tubercle": {"density_fmol_mg": 180, "cell_type": "medium_spiny"},
            # Midbrain autoreceptors
            "substantia_nigra_pars_compacta": {"density_fmol_mg": 85, "cell_type": "dopaminergic_autoreceptor"},
            "VTA": {"density_fmol_mg": 75, "cell_type": "dopaminergic_autoreceptor"},
            # Moderate
            "globus_pallidus_external": {"density_fmol_mg": 45, "cell_type": "GABAergic"},
            "hypothalamus_arcuate": {"density_fmol_mg": 55, "cell_type": "TIDA_neurons"},
            "pituitary_anterior": {"density_fmol_mg": 120, "cell_type": "lactotrophs"},
            "amygdala": {"density_fmol_mg": 25, "cell_type": "intercalated"},
            # Low cortical
            "prefrontal_cortex": {"density_fmol_mg": 8, "cell_type": "layer_V_interneurons"},
            "temporal_cortex": {"density_fmol_mg": 6, "cell_type": "sparse"},
            "hippocampus": {"density_fmol_mg": 5, "cell_type": "minimal"},
            "thalamus_mediodorsal": {"density_fmol_mg": 18, "cell_type": "relay"},
        },
        "pet_tracer": "[11C]raclopride, [18F]fallypride",
        "function": "Gi/o_coupled_inhibitory",
        "subcellular": "postsynaptic_and_presynaptic_autoreceptor",
        "associated_effects": ["euphoria", "reward", "motor_control", "prolactin_inhibition", "psychosis_risk"]
    },
    
    "D3": {
        "density_map": {
            # Limbic-preferring
            "nucleus_accumbens_shell": {"density_fmol_mg": 85, "cell_type": "medium_spiny"},
            "islands_of_Calleja": {"density_fmol_mg": 120, "cell_type": "granule_cells"},
            "olfactory_tubercle": {"density_fmol_mg": 65, "cell_type": "medium_spiny"},
            "VTA": {"density_fmol_mg": 45, "cell_type": "dopaminergic"},
            "substantia_nigra": {"density_fmol_mg": 35, "cell_type": "dopaminergic"},
            "hippocampus": {"density_fmol_mg": 18, "cell_type": "CA1_pyramidal"},
            "amygdala": {"density_fmol_mg": 22, "cell_type": "basolateral"},
            "hypothalamus_mammillary": {"density_fmol_mg": 28, "cell_type": "projection"},
            # Very low in dorsal striatum
            "caudate": {"density_fmol_mg": 8, "cell_type": "sparse"},
            "putamen": {"density_fmol_mg": 6, "cell_type": "sparse"},
            "cortex": {"density_fmol_mg": 3, "cell_type": "minimal"},
        },
        "pet_tracer": "[11C]PHNO",
        "function": "Gi/o_coupled",
        "subcellular": "postsynaptic_limbic",
        "associated_effects": ["addiction_vulnerability", "reward_salience", "mood_regulation"]
    },
    
    # ==========================================================================
    # GLUTAMATE RECEPTORS
    # Source: Monaghan & Cotman 1985, Zilles et al. 2002
    # ==========================================================================
    "NMDA_GluN2A": {
        "density_map": {
            # High cortical expression
            "hippocampus_CA1_stratum_radiatum": {"density_fmol_mg": 280, "cell_type": "pyramidal_dendrites"},
            "hippocampus_CA3": {"density_fmol_mg": 220, "cell_type": "pyramidal"},
            "hippocampus_dentate_molecular": {"density_fmol_mg": 240, "cell_type": "granule_dendrites"},
            "cerebral_cortex_layer_II_III": {"density_fmol_mg": 180, "cell_type": "pyramidal"},
            "cerebral_cortex_layer_V": {"density_fmol_mg": 165, "cell_type": "pyramidal"},
            "amygdala_basolateral": {"density_fmol_mg": 145, "cell_type": "principal"},
            "striatum": {"density_fmol_mg": 95, "cell_type": "medium_spiny"},
            "thalamus_sensory": {"density_fmol_mg": 85, "cell_type": "relay"},
            "cerebellum_granule": {"density_fmol_mg": 45, "cell_type": "granule_cells"},
        },
        "pet_tracer": "none_clinical",
        "function": "ionotropic_Ca2+_permeable",
        "subcellular": "synaptic_PSD",
        "associated_effects": ["memory_formation", "LTP", "excitotoxicity"]
    },
    
    "NMDA_GluN2B": {
        "density_map": {
            # Enriched in forebrain
            "prefrontal_cortex": {"density_fmol_mg": 150, "cell_type": "pyramidal_extrasynaptic"},
            "hippocampus_CA1": {"density_fmol_mg": 180, "cell_type": "pyramidal"},
            "anterior_cingulate": {"density_fmol_mg": 140, "cell_type": "layer_II_III"},
            "striatum": {"density_fmol_mg": 120, "cell_type": "medium_spiny"},
            "amygdala": {"density_fmol_mg": 135, "cell_type": "principal"},
            "thalamus": {"density_fmol_mg": 75, "cell_type": "relay"},
            # Low in hindbrain
            "cerebellum": {"density_fmol_mg": 8, "cell_type": "minimal"},
            "brainstem": {"density_fmol_mg": 25, "cell_type": "sparse"},
        },
        "pet_tracer": "[18F]GE-179",
        "function": "ionotropic_extrasynaptic",
        "subcellular": "extrasynaptic_perisynaptic",
        "associated_effects": ["dissociation", "ketamine_target", "antidepressant_mechanism", "memory_suppression"]
    },
    
    "AMPA_GluA1": {
        "density_map": {
            "hippocampus_CA1": {"density_fmol_mg": 450, "cell_type": "pyramidal_spines"},
            "hippocampus_CA3": {"density_fmol_mg": 380, "cell_type": "pyramidal"},
            "cerebral_cortex_layer_II_III": {"density_fmol_mg": 320, "cell_type": "pyramidal"},
            "cerebral_cortex_layer_V": {"density_fmol_mg": 280, "cell_type": "pyramidal"},
            "amygdala": {"density_fmol_mg": 220, "cell_type": "principal"},
            "striatum": {"density_fmol_mg": 180, "cell_type": "medium_spiny"},
            "thalamus": {"density_fmol_mg": 150, "cell_type": "relay"},
            "cerebellum_purkinje": {"density_fmol_mg": 280, "cell_type": "purkinje_cells"},
        },
        "pet_tracer": "none_clinical",
        "function": "ionotropic_fast_excitatory",
        "subcellular": "postsynaptic_PSD",
        "associated_effects": ["fast_transmission", "synaptic_plasticity", "consciousness_arousal"]
    },
    
    # ==========================================================================
    # GABA RECEPTORS
    # Source: Zezula et al. 1988, Waldvogel et al. 1999
    # ==========================================================================
    "GABA_A_alpha1": {
        "density_map": {
            # Widespread high expression
            "cerebral_cortex_layer_IV": {"density_fmol_mg": 280, "cell_type": "stellate"},
            "cerebral_cortex_layer_II_III": {"density_fmol_mg": 220, "cell_type": "pyramidal"},
            "hippocampus_CA1_stratum_oriens": {"density_fmol_mg": 185, "cell_type": "interneurons"},
            "thalamus_reticular": {"density_fmol_mg": 320, "cell_type": "GABAergic"},
            "thalamus_sensory_relay": {"density_fmol_mg": 180, "cell_type": "thalamocortical"},
            "amygdala_basolateral": {"density_fmol_mg": 165, "cell_type": "principal"},
            "globus_pallidus": {"density_fmol_mg": 250, "cell_type": "GABAergic"},
            "substantia_nigra_pars_reticulata": {"density_fmol_mg": 280, "cell_type": "GABAergic"},
            "cerebellum_granule_layer": {"density_fmol_mg": 320, "cell_type": "granule_cells"},
            "inferior_colliculus": {"density_fmol_mg": 220, "cell_type": "auditory"},
        },
        "pet_tracer": "[11C]flumazenil",
        "function": "Cl-_channel_inhibitory",
        "subcellular": "synaptic",
        "associated_effects": ["sedation", "amnesia", "anticonvulsant", "muscle_relaxation"]
    },
    
    "GABA_A_alpha2": {
        "density_map": {
            # Limbic enriched - anxiolytic target
            "hippocampus_pyramidal_layer": {"density_fmol_mg": 145, "cell_type": "pyramidal_AIS"},
            "amygdala_central": {"density_fmol_mg": 120, "cell_type": "output_neurons"},
            "amygdala_basolateral": {"density_fmol_mg": 95, "cell_type": "principal"},
            "anterior_cingulate": {"density_fmol_mg": 85, "cell_type": "layer_II_III"},
            "hypothalamus": {"density_fmol_mg": 75, "cell_type": "various"},
            "spinal_cord_dorsal_horn": {"density_fmol_mg": 110, "cell_type": "interneurons"},
            # Lower elsewhere
            "thalamus": {"density_fmol_mg": 35, "cell_type": "sparse"},
            "cerebellum": {"density_fmol_mg": 15, "cell_type": "minimal"},
        },
        "pet_tracer": "none_selective",
        "function": "Cl-_channel_inhibitory",
        "subcellular": "axon_initial_segment",
        "associated_effects": ["anxiolysis", "muscle_relaxation", "minimal_sedation"]
    },
    
    # ==========================================================================
    # OPIOID RECEPTORS
    # Source: Pfeiffer et al. 1982, Frost et al. 1985, Henriksen & Bhattacharyya 2020
    # ==========================================================================
    "MOR": {
        "density_map": {
            # Pain/reward circuitry
            "periaqueductal_gray": {"density_fmol_mg": 280, "cell_type": "GABAergic"},
            "thalamus_medial": {"density_fmol_mg": 185, "cell_type": "relay"},
            "nucleus_accumbens_patch": {"density_fmol_mg": 165, "cell_type": "medium_spiny"},
            "caudate_striosomes": {"density_fmol_mg": 145, "cell_type": "medium_spiny"},
            "amygdala_central": {"density_fmol_mg": 155, "cell_type": "GABAergic"},
            "amygdala_basolateral": {"density_fmol_mg": 85, "cell_type": "principal"},
            "habenula_medial": {"density_fmol_mg": 175, "cell_type": "cholinergic"},
            "locus_coeruleus": {"density_fmol_mg": 120, "cell_type": "noradrenergic"},
            "rostral_ventromedial_medulla": {"density_fmol_mg": 195, "cell_type": "ON_OFF_cells"},
            "spinal_cord_dorsal_horn_I_II": {"density_fmol_mg": 220, "cell_type": "interneurons"},
            "nucleus_tractus_solitarius": {"density_fmol_mg": 145, "cell_type": "visceral_afferent"},
            # Cortical - low
            "prefrontal_cortex": {"density_fmol_mg": 25, "cell_type": "interneurons"},
            "anterior_cingulate": {"density_fmol_mg": 45, "cell_type": "layer_II_III"},
            "insula": {"density_fmol_mg": 55, "cell_type": "interoceptive"},
        },
        "pet_tracer": "[11C]carfentanil, [11C]diprenorphine",
        "function": "Gi/o_coupled_inhibitory",
        "subcellular": "presynaptic_and_postsynaptic",
        "associated_effects": ["euphoria", "analgesia", "respiratory_depression", "constipation", "dependence"]
    },
    
    "KOR": {
        "density_map": {
            # Unique high claustrum density
            "claustrum": {"density_fmol_mg": 320, "cell_type": "principal_neurons"},
            "nucleus_accumbens_shell": {"density_fmol_mg": 185, "cell_type": "medium_spiny"},
            "hypothalamus_paraventricular": {"density_fmol_mg": 165, "cell_type": "CRH_neurons"},
            "hypothalamus_arcuate": {"density_fmol_mg": 145, "cell_type": "POMC_neurons"},
            "amygdala_basolateral": {"density_fmol_mg": 125, "cell_type": "principal"},
            "amygdala_central": {"density_fmol_mg": 95, "cell_type": "GABAergic"},
            "bed_nucleus_stria_terminalis": {"density_fmol_mg": 135, "cell_type": "GABAergic"},
            "VTA": {"density_fmol_mg": 85, "cell_type": "dopaminergic"},
            "substantia_nigra": {"density_fmol_mg": 75, "cell_type": "dopaminergic"},
            "periaqueductal_gray": {"density_fmol_mg": 95, "cell_type": "various"},
            "hippocampus_CA3": {"density_fmol_mg": 65, "cell_type": "mossy_fiber_terminals"},
            "prefrontal_cortex": {"density_fmol_mg": 45, "cell_type": "layer_V"},
            "spinal_cord": {"density_fmol_mg": 110, "cell_type": "dorsal_horn"},
        },
        "pet_tracer": "[11C]GR103545, [11C]LY2795050",
        "function": "Gi/o_coupled_inhibitory",
        "subcellular": "presynaptic_DA_terminals",
        "associated_effects": ["dysphoria", "depersonalization", "derealization", "salvia_experience", "stress_response", "analgesia"]
    },
    
    "DOR": {
        "density_map": {
            # Complementary to MOR
            "olfactory_bulb": {"density_fmol_mg": 185, "cell_type": "mitral_cells"},
            "neocortex_layer_II_III": {"density_fmol_mg": 65, "cell_type": "pyramidal"},
            "neocortex_layer_V_VI": {"density_fmol_mg": 55, "cell_type": "pyramidal"},
            "amygdala_basolateral": {"density_fmol_mg": 95, "cell_type": "principal"},
            "striatum_matrix": {"density_fmol_mg": 85, "cell_type": "medium_spiny"},
            "nucleus_accumbens": {"density_fmol_mg": 75, "cell_type": "medium_spiny"},
            "hippocampus": {"density_fmol_mg": 45, "cell_type": "interneurons"},
            "periaqueductal_gray": {"density_fmol_mg": 65, "cell_type": "various"},
            "pontine_nuclei": {"density_fmol_mg": 120, "cell_type": "projection"},
        },
        "pet_tracer": "[11C]methylnaltrindole",
        "function": "Gi/o_coupled",
        "subcellular": "postsynaptic",
        "associated_effects": ["anxiolysis", "antidepressant", "mild_analgesia", "seizure_modulation"]
    },
    
    # ==========================================================================
    # CANNABINOID RECEPTORS
    # Source: Glass et al. 1997, Herkenham et al. 1990, Ceccarini et al. 2015
    # ==========================================================================
    "CB1": {
        "density_map": {
            # Extremely high in specific regions
            "globus_pallidus_internal": {"density_fmol_mg": 850, "cell_type": "striatopallidal_terminals"},
            "globus_pallidus_external": {"density_fmol_mg": 780, "cell_type": "striatopallidal_terminals"},
            "substantia_nigra_pars_reticulata": {"density_fmol_mg": 720, "cell_type": "striatonigral_terminals"},
            "cerebellum_molecular_layer": {"density_fmol_mg": 680, "cell_type": "parallel_fiber_terminals"},
            # High
            "hippocampus_CA1_stratum_oriens": {"density_fmol_mg": 320, "cell_type": "CCK_interneurons"},
            "hippocampus_CA3": {"density_fmol_mg": 280, "cell_type": "interneurons"},
            "hippocampus_dentate_molecular": {"density_fmol_mg": 245, "cell_type": "perforant_path"},
            "striatum": {"density_fmol_mg": 185, "cell_type": "medium_spiny_axon_collaterals"},
            "amygdala_basolateral": {"density_fmol_mg": 165, "cell_type": "CCK_interneurons"},
            "prefrontal_cortex_layer_II_III": {"density_fmol_mg": 145, "cell_type": "CCK_basket_cells"},
            "anterior_cingulate": {"density_fmol_mg": 135, "cell_type": "interneurons"},
            "insula": {"density_fmol_mg": 125, "cell_type": "interneurons"},
            # Moderate
            "hypothalamus": {"density_fmol_mg": 65, "cell_type": "various"},
            "periaqueductal_gray": {"density_fmol_mg": 55, "cell_type": "various"},
            # Low
            "thalamus": {"density_fmol_mg": 25, "cell_type": "sparse"},
            "brainstem_respiratory": {"density_fmol_mg": 15, "cell_type": "minimal"},
        },
        "pet_tracer": "[18F]MK-9470, [11C]OMAR",
        "function": "Gi/o_retrograde_messenger",
        "subcellular": "presynaptic_GABAergic_and_glutamatergic",
        "associated_effects": ["time_dilation", "memory_impairment", "appetite_stimulation", "anxiolysis", "motor_impairment"]
    },
    
    # ==========================================================================
    # ADRENERGIC RECEPTORS
    # Source: Palacios & Kuhar 1980, Aoki et al. 1994
    # ==========================================================================
    "alpha2A": {
        "density_map": {
            # Autoreceptor on noradrenergic neurons
            "locus_coeruleus": {"density_fmol_mg": 420, "cell_type": "noradrenergic_autoreceptor"},
            # Postsynaptic
            "prefrontal_cortex_layer_II": {"density_fmol_mg": 145, "cell_type": "pyramidal_spines"},
            "prefrontal_cortex_layer_V": {"density_fmol_mg": 125, "cell_type": "pyramidal"},
            "anterior_cingulate": {"density_fmol_mg": 115, "cell_type": "layer_II_III"},
            "hippocampus_CA1": {"density_fmol_mg": 95, "cell_type": "pyramidal"},
            "amygdala_basolateral": {"density_fmol_mg": 85, "cell_type": "principal"},
            "hypothalamus": {"density_fmol_mg": 105, "cell_type": "various"},
            "spinal_cord_dorsal_horn": {"density_fmol_mg": 165, "cell_type": "substantia_gelatinosa"},
            "nucleus_tractus_solitarius": {"density_fmol_mg": 135, "cell_type": "visceral"},
            "thalamus": {"density_fmol_mg": 55, "cell_type": "relay"},
            "striatum": {"density_fmol_mg": 25, "cell_type": "sparse"},
        },
        "pet_tracer": "[11C]yohimbine",
        "function": "Gi/o_coupled_inhibitory",
        "subcellular": "presynaptic_autoreceptor_and_postsynaptic",
        "associated_effects": ["sedation", "analgesia", "hypotension", "working_memory_enhancement"]
    },
    
    # ==========================================================================
    # SIGMA RECEPTORS
    # Source: Walker et al. 1990, Hashimoto & Bhattacharyya 2017
    # ==========================================================================
    "Sigma1": {
        "density_map": {
            # ER-associated chaperone - widespread
            "motor_cortex_layer_V": {"density_fmol_mg": 165, "cell_type": "pyramidal"},
            "hippocampus_pyramidal": {"density_fmol_mg": 145, "cell_type": "CA1_CA3_pyramidal"},
            "hippocampus_dentate": {"density_fmol_mg": 125, "cell_type": "granule"},
            "cerebellum_purkinje": {"density_fmol_mg": 185, "cell_type": "purkinje_cells"},
            "cerebellum_granule": {"density_fmol_mg": 135, "cell_type": "granule_cells"},
            "hypothalamus": {"density_fmol_mg": 95, "cell_type": "neuroendocrine"},
            "red_nucleus": {"density_fmol_mg": 155, "cell_type": "motor"},
            "substantia_nigra": {"density_fmol_mg": 85, "cell_type": "dopaminergic"},
            "dorsal_raphe": {"density_fmol_mg": 75, "cell_type": "serotonergic"},
            "spinal_cord_motor_neurons": {"density_fmol_mg": 175, "cell_type": "alpha_motor_neurons"},
            "prefrontal_cortex": {"density_fmol_mg": 65, "cell_type": "pyramidal"},
            "striatum": {"density_fmol_mg": 55, "cell_type": "medium_spiny"},
        },
        "pet_tracer": "[11C]SA4503, [18F]FTC-146",
        "function": "ER_chaperone_Ca2+_modulator",
        "subcellular": "ER_MAM_plasma_membrane",
        "associated_effects": ["neuroprotection", "neuroplasticity", "dissociation_mild", "antidepressant"]
    },
    
    # ==========================================================================
    # ACETYLCHOLINE RECEPTORS
    # Source: Paterson & Bhattacharyya 1999, Pimlott et al. 2004
    # ==========================================================================
    "nAChR_alpha4beta2": {
        "density_map": {
            # High thalamic density
            "thalamus_anteroventral": {"density_fmol_mg": 285, "cell_type": "relay"},
            "thalamus_lateral_geniculate": {"density_fmol_mg": 265, "cell_type": "relay"},
            "thalamus_mediodorsal": {"density_fmol_mg": 225, "cell_type": "relay"},
            "superior_colliculus": {"density_fmol_mg": 195, "cell_type": "visual"},
            "interpeduncular_nucleus": {"density_fmol_mg": 320, "cell_type": "habenular_target"},
            "medial_habenula": {"density_fmol_mg": 245, "cell_type": "cholinergic"},
            # Cortical
            "prefrontal_cortex": {"density_fmol_mg": 85, "cell_type": "layer_V"},
            "visual_cortex": {"density_fmol_mg": 75, "cell_type": "layer_IV"},
            "hippocampus": {"density_fmol_mg": 65, "cell_type": "interneurons"},
            "amygdala": {"density_fmol_mg": 55, "cell_type": "basolateral"},
            "striatum": {"density_fmol_mg": 45, "cell_type": "cholinergic_interneurons"},
            "cerebellum": {"density_fmol_mg": 15, "cell_type": "minimal"},
        },
        "pet_tracer": "[18F]flubatine, [18F]nifene",
        "function": "cation_channel_excitatory",
        "subcellular": "presynaptic_and_postsynaptic",
        "associated_effects": ["attention_enhancement", "memory_consolidation", "nicotine_addiction", "arousal"]
    },
    
    "M1_muscarinic": {
        "density_map": {
            # Cortical and striatal
            "cerebral_cortex_layer_II_III": {"density_fmol_mg": 185, "cell_type": "pyramidal"},
            "cerebral_cortex_layer_V": {"density_fmol_mg": 165, "cell_type": "pyramidal"},
            "hippocampus_CA1": {"density_fmol_mg": 195, "cell_type": "pyramidal"},
            "hippocampus_CA3": {"density_fmol_mg": 175, "cell_type": "pyramidal"},
            "striatum": {"density_fmol_mg": 225, "cell_type": "medium_spiny"},
            "amygdala": {"density_fmol_mg": 145, "cell_type": "principal"},
            "olfactory_bulb": {"density_fmol_mg": 125, "cell_type": "mitral"},
            # Low
            "thalamus": {"density_fmol_mg": 35, "cell_type": "sparse"},
            "cerebellum": {"density_fmol_mg": 12, "cell_type": "minimal"},
            "brainstem": {"density_fmol_mg": 25, "cell_type": "sparse"},
        },
        "pet_tracer": "[11C]LSN3172176",
        "function": "Gq_coupled_excitatory",
        "subcellular": "postsynaptic",
        "associated_effects": ["memory_enhancement", "psychosis", "anticholinergic_delirium", "cognition"]
    }
}


# Known pharmacology for common substances
SUBSTANCE_PHARMACOLOGY = {
    "LSD": SubstancePharmacology(
        name="LSD",
        bindings=[
            ReceptorBinding("5-HT2A", ki_nm=2.9, activity=ReceptorActivity.PARTIAL_AGONIST, efficacy=0.4),
            ReceptorBinding("5-HT2B", ki_nm=4.9, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("5-HT2C", ki_nm=5.5, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("5-HT1A", ki_nm=1.1, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("D2", ki_nm=120, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("D1", ki_nm=180, activity=ReceptorActivity.PARTIAL_AGONIST),
        ],
        primary_targets=["5-HT2A", "5-HT1A", "D2"],
        half_life_hours=3.0,
        bioavailability={"oral": 0.71, "sublingual": 0.71}
    ),
    
    "Psilocin": SubstancePharmacology(
        name="Psilocin",
        bindings=[
            ReceptorBinding("5-HT2A", ki_nm=6.0, activity=ReceptorActivity.PARTIAL_AGONIST, efficacy=0.6),
            ReceptorBinding("5-HT2C", ki_nm=97, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("5-HT1A", ki_nm=123, activity=ReceptorActivity.PARTIAL_AGONIST),
        ],
        primary_targets=["5-HT2A"],
        half_life_hours=2.5
    ),
    
    "DMT": SubstancePharmacology(
        name="DMT",
        bindings=[
            ReceptorBinding("5-HT2A", ki_nm=75, activity=ReceptorActivity.AGONIST, efficacy=0.8),
            ReceptorBinding("5-HT2C", ki_nm=360, activity=ReceptorActivity.AGONIST),
            ReceptorBinding("Sigma1", ki_nm=14000, activity=ReceptorActivity.AGONIST),
        ],
        primary_targets=["5-HT2A", "Sigma1"],
        half_life_hours=0.25  # Very short when IV/vaporized
    ),
    
    "MDMA": SubstancePharmacology(
        name="MDMA",
        bindings=[
            ReceptorBinding("SERT", ki_nm=238, activity=ReceptorActivity.RELEASING_AGENT),
            ReceptorBinding("NET", ki_nm=463, activity=ReceptorActivity.RELEASING_AGENT),
            ReceptorBinding("DAT", ki_nm=1572, activity=ReceptorActivity.RELEASING_AGENT),
            ReceptorBinding("5-HT2A", ki_nm=5115, activity=ReceptorActivity.AGONIST),
        ],
        primary_targets=["SERT", "NET"],
        half_life_hours=7.0
    ),
    
    "Ketamine": SubstancePharmacology(
        name="Ketamine",
        bindings=[
            ReceptorBinding("NMDA", ki_nm=659, activity=ReceptorActivity.ANTAGONIST),
            ReceptorBinding("D2", ki_nm=2640, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("MOR", ki_nm=25000, activity=ReceptorActivity.AGONIST),
            ReceptorBinding("Sigma1", ki_nm=26000, activity=ReceptorActivity.AGONIST),
        ],
        primary_targets=["NMDA"],
        half_life_hours=2.5
    ),
    
    "Salvinorin_A": SubstancePharmacology(
        name="Salvinorin A",
        bindings=[
            ReceptorBinding("KOR", ki_nm=1.9, activity=ReceptorActivity.FULL_AGONIST, efficacy=1.0),
        ],
        primary_targets=["KOR"],
        half_life_hours=0.5  # Very short
    ),
    
    "Mescaline": SubstancePharmacology(
        name="Mescaline",
        bindings=[
            ReceptorBinding("5-HT2A", ki_nm=6200, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("5-HT2B", ki_nm=870, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("5-HT2C", ki_nm=1300, activity=ReceptorActivity.PARTIAL_AGONIST),
        ],
        primary_targets=["5-HT2A", "5-HT2B", "5-HT2C"],
        half_life_hours=6.0
    ),
    
    "Cannabis_THC": SubstancePharmacology(
        name="THC",
        bindings=[
            ReceptorBinding("CB1", ki_nm=40.7, activity=ReceptorActivity.PARTIAL_AGONIST),
            ReceptorBinding("CB2", ki_nm=36.4, activity=ReceptorActivity.PARTIAL_AGONIST),
        ],
        primary_targets=["CB1"],
        half_life_hours=25  # Lipophilic, long elimination
    )
}


# Effect-to-receptor-to-region mapping hypothesis
EFFECT_MECHANISM_MAPPING = {
    # Visual effects
    "geometry": {
        "primary_receptor": "5-HT2A",
        "brain_regions": ["visual_cortex_V1", "visual_association_cortex"],
        "mechanism": "Disruption of hierarchical predictive coding in visual cortex",
        "confidence": 0.9
    },
    "color_enhancement": {
        "primary_receptor": "5-HT2A",
        "brain_regions": ["visual_cortex_V4"],
        "mechanism": "Enhanced color processing pathway activation",
        "confidence": 0.8
    },
    "hallucination": {
        "primary_receptor": "5-HT2A",
        "brain_regions": ["visual_cortex", "temporal_cortex", "prefrontal_cortex"],
        "mechanism": "Disrupted reality monitoring, increased internal model dominance",
        "confidence": 0.85
    },
    
    # Cognitive effects
    "ego_dissolution": {
        "primary_receptor": "5-HT2A",
        "brain_regions": ["default_mode_network", "claustrum", "posterior_cingulate", "medial_prefrontal_cortex"],
        "mechanism": "Disruption of self-referential processing, DMN desynchronization",
        "confidence": 0.9
    },
    "time_dilation": {
        "primary_receptor": "CB1",
        "brain_regions": ["cerebellum", "basal_ganglia", "prefrontal_cortex"],
        "mechanism": "Disruption of internal clock and temporal processing",
        "secondary_receptors": ["5-HT2A"],
        "confidence": 0.7
    },
    "introspection": {
        "primary_receptor": "5-HT1A",
        "brain_regions": ["prefrontal_cortex", "anterior_cingulate"],
        "mechanism": "Enhanced metacognitive monitoring",
        "confidence": 0.6
    },
    "thought_loops": {
        "primary_receptor": "5-HT2A",
        "brain_regions": ["prefrontal_cortex", "cingulate_cortex"],
        "mechanism": "Disrupted executive function, working memory overload",
        "confidence": 0.7
    },
    
    # Emotional effects
    "euphoria_stimulant": {
        "primary_receptor": "D2",
        "brain_regions": ["nucleus_accumbens", "VTA", "prefrontal_cortex"],
        "mechanism": "Increased dopamine signaling in reward pathway",
        "confidence": 0.95
    },
    "euphoria_opioid": {
        "primary_receptor": "MOR",
        "brain_regions": ["nucleus_accumbens", "VTA", "periaqueductal_gray"],
        "mechanism": "Opioid receptor activation, dopamine disinhibition",
        "confidence": 0.95
    },
    "dysphoria": {
        "primary_receptor": "KOR",
        "brain_regions": ["nucleus_accumbens", "amygdala", "prefrontal_cortex"],
        "mechanism": "Kappa opioid-mediated aversion, dopamine suppression",
        "confidence": 0.9
    },
    
    # Physical effects
    "dissociation": {
        "primary_receptor": "NMDA",
        "brain_regions": ["hippocampus", "cortex_widespread", "thalamus"],
        "mechanism": "Disruption of glutamatergic signaling, cortical disconnection",
        "confidence": 0.95
    },
    "sedation": {
        "primary_receptor": "GABA_A",
        "brain_regions": ["thalamus", "cortex_widespread", "brainstem"],
        "mechanism": "Enhanced GABAergic inhibition, reduced arousal",
        "confidence": 0.95
    },
    "stimulation": {
        "primary_receptor": "D2",
        "brain_regions": ["striatum", "prefrontal_cortex", "motor_cortex"],
        "mechanism": "Increased catecholamine signaling",
        "secondary_receptors": ["NET", "DAT"],
        "confidence": 0.9
    }
}


class PharmacologyAnalyzer:
    """Analyze receptor binding profiles and predict effects using quantitative density maps"""
    
    def __init__(self):
        self.receptor_regions = RECEPTOR_BRAIN_REGIONS
        self.substance_data = SUBSTANCE_PHARMACOLOGY
        self.effect_mappings = EFFECT_MECHANISM_MAPPING
    
    def get_receptor_profile(self, substance: str) -> Optional[SubstancePharmacology]:
        """Get pharmacology data for a substance"""
        return self.substance_data.get(substance)
    
    def get_receptor_density_map(self, receptor: str) -> Dict[str, Dict]:
        """Get detailed density map for a receptor"""
        receptor_info = self.receptor_regions.get(receptor, {})
        return receptor_info.get("density_map", {})
    
    def get_high_density_regions(self, receptor: str, threshold_fmol_mg: float = 100) -> List[Tuple[str, float]]:
        """Get regions with density above threshold, sorted by density"""
        density_map = self.get_receptor_density_map(receptor)
        high_regions = [
            (region, data["density_fmol_mg"])
            for region, data in density_map.items()
            if data.get("density_fmol_mg", 0) >= threshold_fmol_mg
        ]
        return sorted(high_regions, key=lambda x: x[1], reverse=True)
    
    def get_brain_regions_for_substance(self, substance: str) -> Dict[str, List[Dict]]:
        """
        Get brain regions likely affected by a substance with quantitative estimates.
        Returns regions weighted by receptor density * binding affinity.
        """
        pharm = self.get_receptor_profile(substance)
        if not pharm:
            return {}
        
        region_activation = {}
        
        for binding in pharm.bindings:
            receptor_info = self.receptor_regions.get(binding.receptor, {})
            density_map = receptor_info.get("density_map", {})
            
            # Calculate occupancy estimate: higher affinity (lower Ki) = higher occupancy
            # Using simplified Hill equation approximation
            if binding.ki_nm and binding.ki_nm > 0:
                # Relative occupancy score (normalized, not actual %)
                affinity_score = 100 / (binding.ki_nm + 100)
            else:
                affinity_score = 0.5
            
            for region, data in density_map.items():
                density = data.get("density_fmol_mg", 0)
                cell_type = data.get("cell_type", "unknown")
                
                # Activation = density * affinity
                activation = density * affinity_score
                
                if region not in region_activation:
                    region_activation[region] = {
                        "total_activation": 0,
                        "receptors": [],
                        "cell_types": set()
                    }
                
                region_activation[region]["total_activation"] += activation
                region_activation[region]["receptors"].append({
                    "receptor": binding.receptor,
                    "density_fmol_mg": density,
                    "ki_nm": binding.ki_nm,
                    "activity": binding.activity.value,
                    "contribution": activation
                })
                region_activation[region]["cell_types"].add(cell_type)
        
        # Convert sets to lists and sort by activation
        for region in region_activation:
            region_activation[region]["cell_types"] = list(region_activation[region]["cell_types"])
        
        # Sort by total activation
        sorted_regions = dict(sorted(
            region_activation.items(),
            key=lambda x: x[1]["total_activation"],
            reverse=True
        ))
        
        return sorted_regions
    
    def predict_effects_from_pharmacology(self, substance: str) -> List[Dict]:
        """Predict likely effects based on receptor binding profile and density"""
        pharm = self.get_receptor_profile(substance)
        if not pharm:
            return []
        
        predicted_effects = []
        
        for binding in pharm.bindings:
            receptor_info = self.receptor_regions.get(binding.receptor, {})
            associated_effects = receptor_info.get("associated_effects", [])
            density_map = receptor_info.get("density_map", {})
            
            # Calculate max density for this receptor
            max_density = max(
                (d.get("density_fmol_mg", 0) for d in density_map.values()),
                default=0
            )
            
            for effect in associated_effects:
                # Intensity based on affinity and receptor density
                intensity = 1.0
                if binding.ki_nm:
                    if binding.ki_nm < 10:
                        intensity = 0.95
                    elif binding.ki_nm < 50:
                        intensity = 0.85
                    elif binding.ki_nm < 100:
                        intensity = 0.7
                    elif binding.ki_nm < 500:
                        intensity = 0.5
                    elif binding.ki_nm < 1000:
                        intensity = 0.35
                    else:
                        intensity = 0.2
                
                # Modulate by receptor density (high density = more effect)
                density_factor = min(max_density / 200, 1.5)  # Cap at 1.5x
                intensity *= density_factor
                
                # Get primary brain regions for this effect
                primary_regions = self.get_high_density_regions(binding.receptor, threshold_fmol_mg=50)[:5]
                
                predicted_effects.append({
                    "effect": effect,
                    "receptor": binding.receptor,
                    "predicted_intensity": min(intensity, 1.0),
                    "ki_nm": binding.ki_nm,
                    "activity": binding.activity.value,
                    "pet_tracer": receptor_info.get("pet_tracer"),
                    "primary_regions": [r[0] for r in primary_regions],
                    "max_density_fmol_mg": max_density
                })
        
        # Sort by predicted intensity
        return sorted(predicted_effects, key=lambda x: x["predicted_intensity"], reverse=True)
    
    def get_receptor_by_region(self, region_name: str) -> List[Dict]:
        """Get all receptors expressed in a specific brain region"""
        results = []
        
        for receptor, info in self.receptor_regions.items():
            density_map = info.get("density_map", {})
            
            # Search for region (partial match)
            for region, data in density_map.items():
                if region_name.lower() in region.lower():
                    results.append({
                        "receptor": receptor,
                        "region": region,
                        "density_fmol_mg": data.get("density_fmol_mg"),
                        "cell_type": data.get("cell_type"),
                        "function": info.get("function"),
                        "associated_effects": info.get("associated_effects", [])
                    })
        
        return sorted(results, key=lambda x: x.get("density_fmol_mg", 0), reverse=True)
    
    def compare_substances(self, sub1: str, sub2: str) -> Dict:
        """Compare pharmacological profiles of two substances"""
        pharm1 = self.get_receptor_profile(sub1)
        pharm2 = self.get_receptor_profile(sub2)
        
        if not pharm1 or not pharm2:
            return {"error": "Substance not found"}
        
        receptors1 = {b.receptor: b for b in pharm1.bindings}
        receptors2 = {b.receptor: b for b in pharm2.bindings}
        
        all_receptors = set(receptors1.keys()) | set(receptors2.keys())
        
        comparison = []
        for receptor in all_receptors:
            b1 = receptors1.get(receptor)
            b2 = receptors2.get(receptor)
            
            comparison.append({
                "receptor": receptor,
                f"{sub1}_ki": b1.ki_nm if b1 else None,
                f"{sub1}_activity": b1.activity.value if b1 else None,
                f"{sub2}_ki": b2.ki_nm if b2 else None,
                f"{sub2}_activity": b2.activity.value if b2 else None
            })
        
        return {
            "substances": [sub1, sub2],
            "comparison": comparison,
            "shared_targets": list(set(receptors1.keys()) & set(receptors2.keys())),
            "unique_to_1": list(set(receptors1.keys()) - set(receptors2.keys())),
            "unique_to_2": list(set(receptors2.keys()) - set(receptors1.keys()))
        }
    
    def calculate_regional_activation_profile(
        self, 
        substance: str,
        normalize: bool = True
    ) -> pd.DataFrame:
        """
        Calculate detailed regional activation profile as a DataFrame.
        Useful for creating brain maps.
        """
        regions = self.get_brain_regions_for_substance(substance)
        
        rows = []
        for region, data in regions.items():
            for receptor_contrib in data["receptors"]:
                rows.append({
                    "region": region,
                    "receptor": receptor_contrib["receptor"],
                    "density_fmol_mg": receptor_contrib["density_fmol_mg"],
                    "ki_nm": receptor_contrib["ki_nm"],
                    "activity": receptor_contrib["activity"],
                    "contribution": receptor_contrib["contribution"],
                    "cell_types": ", ".join(data["cell_types"])
                })
        
        df = pd.DataFrame(rows)
        
        if normalize and len(df) > 0:
            max_contrib = df["contribution"].max()
            if max_contrib > 0:
                df["normalized_contribution"] = df["contribution"] / max_contrib
        
        return df


def create_receptor_effect_dataframe() -> pd.DataFrame:
    """Create DataFrame mapping receptors to effects and regions with quantitative density"""
    rows = []
    
    for receptor, info in RECEPTOR_BRAIN_REGIONS.items():
        density_map = info.get("density_map", {})
        effects = info.get("associated_effects", [])
        
        for region, region_data in density_map.items():
            for effect in effects:
                rows.append({
                    "receptor": receptor,
                    "effect": effect,
                    "brain_region": region,
                    "density_fmol_mg": region_data.get("density_fmol_mg"),
                    "cell_type": region_data.get("cell_type"),
                    "function": info.get("function"),
                    "pet_tracer": info.get("pet_tracer"),
                    "subcellular": info.get("subcellular")
                })
    
    return pd.DataFrame(rows)


def create_density_heatmap_data() -> pd.DataFrame:
    """Create pivoted DataFrame suitable for heatmap visualization"""
    rows = []
    
    all_regions = set()
    for receptor, info in RECEPTOR_BRAIN_REGIONS.items():
        density_map = info.get("density_map", {})
        all_regions.update(density_map.keys())
    
    for receptor, info in RECEPTOR_BRAIN_REGIONS.items():
        density_map = info.get("density_map", {})
        row = {"receptor": receptor}
        for region in all_regions:
            if region in density_map:
                row[region] = density_map[region].get("density_fmol_mg", 0)
            else:
                row[region] = 0
        rows.append(row)
    
    return pd.DataFrame(rows).set_index("receptor")


if __name__ == "__main__":
    analyzer = PharmacologyAnalyzer()
    
    # Example: Analyze LSD with new detailed density data
    print("=" * 70)
    print("LSD PHARMACOLOGY ANALYSIS (Quantitative Density Maps)")
    print("=" * 70)
    
    lsd = analyzer.get_receptor_profile("LSD")
    if lsd:
        print(f"\nSubstance: {lsd.name}")
        print(f"Half-life: {lsd.half_life_hours}h")
        print(f"Primary targets: {lsd.primary_targets}")
        print("\nReceptor Bindings (sorted by affinity):")
        for b in sorted(lsd.bindings, key=lambda x: x.ki_nm or float('inf')):
            print(f"  {b.receptor}: Ki={b.ki_nm}nM ({b.activity.value})")
    
    print("\n" + "-" * 70)
    print("BRAIN REGIONS BY ACTIVATION (density × affinity)")
    print("-" * 70)
    regions = analyzer.get_brain_regions_for_substance("LSD")
    for i, (region, data) in enumerate(list(regions.items())[:15]):
        receptors = ", ".join([r["receptor"] for r in data["receptors"]])
        print(f"  {i+1:2}. {region}: {data['total_activation']:.1f}")
        print(f"      via: {receptors}")
    
    print("\n" + "-" * 70)
    print("PREDICTED EFFECTS (by intensity)")
    print("-" * 70)
    effects = analyzer.predict_effects_from_pharmacology("LSD")
    for e in effects[:12]:
        regions_str = ", ".join(e["primary_regions"][:3])
        print(f"  {e['effect']}: {e['predicted_intensity']:.2f}")
        print(f"    via {e['receptor']} (Ki={e['ki_nm']}nM), regions: {regions_str}")
    
    print("\n" + "-" * 70)
    print("CLAUSTRUM RECEPTOR EXPRESSION (high KOR density - salvia target)")
    print("-" * 70)
    claustrum_receptors = analyzer.get_receptor_by_region("claustrum")
    for r in claustrum_receptors[:8]:
        print(f"  {r['receptor']}: {r['density_fmol_mg']} fmol/mg ({r['cell_type']})")
    
    print("\n" + "-" * 70)
    print("5-HT2A HIGH-DENSITY REGIONS (psychedelic primary target)")
    print("-" * 70)
    ht2a_regions = analyzer.get_high_density_regions("5-HT2A", threshold_fmol_mg=25)
    for region, density in ht2a_regions:
        print(f"  {region}: {density} fmol/mg")
    
    print("\n" + "-" * 70)
    print("SUBSTANCE COMPARISON: LSD vs Ketamine")
    print("-" * 70)
    comparison = analyzer.compare_substances("LSD", "Ketamine")
    print(f"Shared targets: {comparison['shared_targets']}")
    print(f"Unique to LSD: {comparison['unique_to_1']}")
    print(f"Unique to Ketamine: {comparison['unique_to_2']}")
