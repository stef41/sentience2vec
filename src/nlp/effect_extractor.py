"""
LLM-based Effect Extraction
Uses local LLMs (via Ollama) to extract structured phenomenological data
"""

import json
import re
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict, field
from enum import Enum
import ollama


class EffectCategory(Enum):
    VISUAL = "visual"
    AUDITORY = "auditory"
    TACTILE = "tactile"
    COGNITIVE = "cognitive"
    PHYSICAL = "physical"
    EMOTIONAL = "emotional"
    TRANSPERSONAL = "transpersonal"


class Valence(Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"
    NEUTRAL = "neutral"
    MIXED = "mixed"


@dataclass
class ExtractedEffect:
    """Structured representation of an extracted subjective effect"""
    name: str                          # Standardized effect name
    category: EffectCategory           # Effect category
    intensity: float                   # 0-1 scale
    valence: Valence                   # Emotional quality
    onset_time: Optional[str] = None   # When it occurred
    duration: Optional[str] = None     # How long it lasted
    body_location: Optional[str] = None  # If localized
    description: Optional[str] = None  # Original text excerpt
    confidence: float = 1.0            # Extraction confidence


@dataclass
class TemporalPhase:
    """Phase of experience with associated effects"""
    phase: str  # onset, comeup, peak, plateau, offset, afterglow
    start_time: Optional[str] = None
    effects: List[ExtractedEffect] = field(default_factory=list)


@dataclass
class ExtractedReport:
    """Fully extracted experience report"""
    report_id: str
    substance: str
    dose: Optional[str]
    route: Optional[str]
    set_description: Optional[str]      # Mindset going in
    setting_description: Optional[str]  # Physical environment
    phases: List[TemporalPhase] = field(default_factory=list)
    overall_valence: Optional[Valence] = None
    key_insights: List[str] = field(default_factory=list)
    body_parts_mentioned: List[str] = field(default_factory=list)
    sensory_modalities: List[str] = field(default_factory=list)


# PsychonautWiki Subjective Effect Index - standardized effect taxonomy
EFFECT_TAXONOMY = {
    "visual": [
        "color_enhancement", "pattern_recognition_enhancement", "visual_acuity_enhancement",
        "drifting", "morphing", "breathing", "melting", "flowing",
        "geometry", "tracers", "after_images", "color_shifting",
        "symmetrical_texture_repetition", "recursion", "scenery_slicing",
        "internal_hallucination", "external_hallucination",
        "autonomous_entity", "shadow_people", "transformations"
    ],
    "auditory": [
        "auditory_enhancement", "auditory_distortion", "auditory_hallucination",
        "music_enhancement"
    ],
    "cognitive": [
        "analysis_enhancement", "anxiety", "conceptual_thinking", 
        "creativity_enhancement", "cognitive_euphoria", "cognitive_dysphoria",
        "delusion", "deja_vu", "ego_dissolution", "ego_death",
        "emotion_enhancement", "empathy_enhancement", "focus_enhancement",
        "immersion_enhancement", "introspection", "memory_suppression",
        "novelty_enhancement", "personal_meaning_enhancement",
        "thought_acceleration", "thought_connectivity", "thought_loops",
        "time_dilation", "time_compression", "unity_interconnectedness"
    ],
    "physical": [
        "stimulation", "sedation", "physical_euphoria", "body_high",
        "tactile_enhancement", "pain_relief", "muscle_relaxation",
        "nausea", "appetite_suppression", "temperature_dysregulation",
        "pupil_dilation", "vasoconstriction", "increased_heart_rate"
    ],
    "transpersonal": [
        "ego_dissolution", "unity_interconnectedness", 
        "existential_self_realization", "spirituality_enhancement",
        "perception_of_eternalism", "perception_of_predeterminism"
    ]
}


class EffectExtractor:
    """Extract structured effects from experience reports using local LLM"""
    
    def __init__(self, model: str = "llama3.2"):
        self.model = model
        self.effect_taxonomy = EFFECT_TAXONOMY
        
    def _build_extraction_prompt(self, text: str) -> str:
        """Build prompt for effect extraction"""
        effects_list = []
        for category, effects in self.effect_taxonomy.items():
            effects_list.extend([f"{category}:{e}" for e in effects])
        
        prompt = f"""You are an expert at analyzing psychedelic experience reports and extracting structured phenomenological data.

Analyze the following experience report and extract:
1. All subjective effects experienced
2. The temporal phase when each occurred (onset, comeup, peak, plateau, offset, afterglow)
3. Intensity of each effect (mild, moderate, strong, extreme)
4. Emotional valence (positive, negative, neutral, mixed)
5. Any body parts or locations mentioned
6. Key insights or realizations

Use ONLY effects from this standardized taxonomy:
{json.dumps(self.effect_taxonomy, indent=2)}

If an effect doesn't match the taxonomy exactly, map it to the closest match and note the original description.

Return your analysis as JSON with this structure:
{{
    "phases": [
        {{
            "phase": "peak",
            "effects": [
                {{
                    "name": "geometry",
                    "category": "visual",
                    "intensity": 0.8,
                    "valence": "positive",
                    "description": "original text excerpt"
                }}
            ]
        }}
    ],
    "overall_valence": "positive",
    "key_insights": ["insight 1", "insight 2"],
    "body_parts_mentioned": ["head", "chest"],
    "set_description": "description of mindset",
    "setting_description": "description of environment"
}}

EXPERIENCE REPORT:
{text[:8000]}

JSON OUTPUT:"""
        
        return prompt
    
    def extract_effects(self, report_text: str, report_id: str = "") -> ExtractedReport:
        """Extract structured effects from report text"""
        prompt = self._build_extraction_prompt(report_text)
        
        try:
            response = ollama.generate(
                model=self.model,
                prompt=prompt,
                options={
                    "temperature": 0.3,  # Lower for more consistent extraction
                    "num_predict": 2000
                }
            )
            
            # Parse JSON from response
            response_text = response["response"]
            
            # Find JSON in response
            json_match = re.search(r'\{[\s\S]*\}', response_text)
            if json_match:
                data = json.loads(json_match.group())
            else:
                data = {}
            
            # Build ExtractedReport
            phases = []
            for phase_data in data.get("phases", []):
                effects = []
                for eff in phase_data.get("effects", []):
                    try:
                        effects.append(ExtractedEffect(
                            name=eff.get("name", "unknown"),
                            category=EffectCategory(eff.get("category", "cognitive")),
                            intensity=float(eff.get("intensity", 0.5)),
                            valence=Valence(eff.get("valence", "neutral")),
                            description=eff.get("description")
                        ))
                    except (ValueError, KeyError):
                        continue
                
                phases.append(TemporalPhase(
                    phase=phase_data.get("phase", "unknown"),
                    effects=effects
                ))
            
            return ExtractedReport(
                report_id=report_id,
                substance="",  # Would be filled from report metadata
                dose=None,
                route=None,
                set_description=data.get("set_description"),
                setting_description=data.get("setting_description"),
                phases=phases,
                overall_valence=Valence(data.get("overall_valence", "neutral")) if data.get("overall_valence") else None,
                key_insights=data.get("key_insights", []),
                body_parts_mentioned=data.get("body_parts_mentioned", []),
                sensory_modalities=[]
            )
            
        except Exception as e:
            print(f"Extraction error: {e}")
            return ExtractedReport(
                report_id=report_id,
                substance="",
                dose=None,
                route=None,
                phases=[],
                key_insights=[]
            )
    
    def batch_extract(
        self, 
        reports: List[Dict],
        output_file: str = "data/processed/extracted_effects.json"
    ) -> List[ExtractedReport]:
        """Extract effects from multiple reports"""
        results = []
        
        for i, report in enumerate(reports):
            print(f"[{i+1}/{len(reports)}] Extracting: {report.get('id', 'unknown')}")
            
            extracted = self.extract_effects(
                report.get("text", ""),
                report_id=report.get("id", str(i))
            )
            extracted.substance = report.get("substance", "")
            
            # Get dose info if available
            substances_detail = report.get("substances_detail", [])
            if substances_detail:
                extracted.dose = substances_detail[0].get("amount")
                extracted.route = substances_detail[0].get("route")
            
            results.append(extracted)
        
        # Save results
        with open(output_file, "w") as f:
            json.dump([asdict(r) for r in results], f, indent=2, default=str)
        
        return results


def compute_effect_frequencies(extracted_reports: List[ExtractedReport]) -> Dict:
    """Compute frequency statistics for effects"""
    effect_counts = {}
    substance_effects = {}
    phase_effects = {}
    
    for report in extracted_reports:
        substance = report.substance
        
        for phase in report.phases:
            phase_name = phase.phase
            
            for effect in phase.effects:
                # Overall counts
                effect_counts[effect.name] = effect_counts.get(effect.name, 0) + 1
                
                # Per-substance counts
                if substance not in substance_effects:
                    substance_effects[substance] = {}
                substance_effects[substance][effect.name] = \
                    substance_effects[substance].get(effect.name, 0) + 1
                
                # Per-phase counts
                if phase_name not in phase_effects:
                    phase_effects[phase_name] = {}
                phase_effects[phase_name][effect.name] = \
                    phase_effects[phase_name].get(effect.name, 0) + 1
    
    return {
        "total_reports": len(extracted_reports),
        "effect_counts": effect_counts,
        "by_substance": substance_effects,
        "by_phase": phase_effects
    }


if __name__ == "__main__":
    # Test extraction on sample text
    sample_text = """
    T+0:00 - Dropped 150ug of LSD on an empty stomach. Feeling nervous but excited.
    
    T+0:30 - Starting to feel the first alerts. Slight body tingling and a shift in 
    perception. Colors seem slightly more vivid.
    
    T+1:00 - Definitely coming up now. Strong body high with waves of euphoria.
    The walls are starting to breathe gently. Geometric patterns appearing when 
    I close my eyes.
    
    T+2:00 - Peak hitting hard. Intense closed-eye visuals - intricate mandalas 
    and fractal patterns in every color. Time feels meaningless. Having deep 
    thoughts about consciousness and the nature of reality.
    
    T+3:00 - Still peaking. Experienced a moment of complete ego dissolution 
    where I felt connected to everything in the universe. Overwhelming but 
    beautiful. Tears of joy.
    
    T+5:00 - Coming down slowly. Visuals less intense but still present. 
    Feeling profoundly grateful and at peace. Deep insights about my life 
    and relationships.
    
    T+8:00 - Mostly baseline. Slight visual enhancement persists. 
    Afterglow of contentment and clarity.
    """
    
    extractor = EffectExtractor()
    
    print("Testing effect extraction...")
    print("(Requires Ollama running with llama3.2 model)")
    
    try:
        result = extractor.extract_effects(sample_text, "test_001")
        print("\n=== Extracted Report ===")
        print(json.dumps(asdict(result), indent=2, default=str))
    except Exception as e:
        print(f"Could not run extraction (is Ollama running?): {e}")
        print("\nTo use this module:")
        print("1. Install Ollama: https://ollama.ai")
        print("2. Pull a model: ollama pull llama3.2")
        print("3. Run Ollama in the background")
