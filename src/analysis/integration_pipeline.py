"""
Integration Pipeline
Connects all modules: data collection → NLP → pharmacology → knowledge graph → analysis
"""

import asyncio
import json
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from pathlib import Path
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class AnalysisResult:
    """Result from integrated analysis"""
    substance: str
    n_reports_analyzed: int
    
    # Effect frequencies from empirical data
    effect_frequencies: Dict[str, float]
    
    # Predicted effects from pharmacology
    predicted_effects: List[Dict]
    
    # Brain regions implicated
    brain_regions: Dict[str, float]
    
    # Receptor binding profile
    receptor_profile: Dict[str, Any]
    
    # Effect co-occurrence patterns
    effect_clusters: List[List[str]]
    
    # Model-data concordance
    prediction_accuracy: Optional[float] = None


class ConsciousnessResearchPipeline:
    """
    Main pipeline integrating all research components
    """
    
    def __init__(
        self,
        data_dir: str = "./data",
        cache_dir: str = "./cache",
        output_dir: str = "./results"
    ):
        self.data_dir = Path(data_dir)
        self.cache_dir = Path(cache_dir)
        self.output_dir = Path(output_dir)
        
        # Create directories
        for d in [self.data_dir, self.cache_dir, self.output_dir]:
            d.mkdir(parents=True, exist_ok=True)
        
        # Initialize components lazily
        self._psychonautwiki = None
        self._tripsit = None
        self._erowid = None
        self._effect_extractor = None
        self._knowledge_graph = None
        self._pharmacology = None
    
    @property
    def psychonautwiki(self):
        if self._psychonautwiki is None:
            from consciousness_study.src.data_collection.psychonautwiki_client import PsychonautWikiAPI
            self._psychonautwiki = PsychonautWikiAPI()
        return self._psychonautwiki
    
    @property
    def tripsit(self):
        if self._tripsit is None:
            from consciousness_study.src.data_collection.tripsit_client import TripSitAPI
            self._tripsit = TripSitAPI()
        return self._tripsit
    
    @property
    def erowid(self):
        if self._erowid is None:
            from consciousness_study.src.data_collection.erowid_scraper import ErowidScraper
            self._erowid = ErowidScraper(cache_dir=str(self.cache_dir / "erowid"))
        return self._erowid
    
    @property
    def effect_extractor(self):
        if self._effect_extractor is None:
            from consciousness_study.src.nlp.effect_extractor import EffectExtractor
            self._effect_extractor = EffectExtractor()
        return self._effect_extractor
    
    @property
    def knowledge_graph(self):
        if self._knowledge_graph is None:
            from consciousness_study.src.knowledge_graph.graph_builder import build_base_graph
            self._knowledge_graph = build_base_graph()
        return self._knowledge_graph
    
    @property
    def pharmacology(self):
        if self._pharmacology is None:
            from consciousness_study.src.pharmacology.receptor_mapping import PharmacologyAnalyzer
            self._pharmacology = PharmacologyAnalyzer()
        return self._pharmacology
    
    async def collect_substance_data(
        self,
        substance_name: str,
        max_reports: int = 100
    ) -> Dict[str, Any]:
        """
        Collect all available data for a substance from all sources
        """
        logger.info(f"Collecting data for {substance_name}...")
        
        data = {
            "substance": substance_name,
            "psychonautwiki": None,
            "tripsit": None,
            "erowid_reports": [],
            "pharmacology": None
        }
        
        # PsychonautWiki data
        try:
            pw_data = await asyncio.to_thread(
                self.psychonautwiki.get_substance, 
                substance_name
            )
            if pw_data:
                data["psychonautwiki"] = asdict(pw_data)
                logger.info(f"  PsychonautWiki: Found {substance_name}")
        except Exception as e:
            logger.warning(f"  PsychonautWiki error: {e}")
        
        # TripSit data
        try:
            ts_data = await asyncio.to_thread(
                self.tripsit.get_drug,
                substance_name
            )
            if ts_data:
                data["tripsit"] = ts_data
                logger.info(f"  TripSit: Found {substance_name}")
        except Exception as e:
            logger.warning(f"  TripSit error: {e}")
        
        # Erowid experience reports
        try:
            reports = await asyncio.to_thread(
                self.erowid.scrape_substance,
                substance_name,
                max_reports=max_reports
            )
            data["erowid_reports"] = [asdict(r) for r in reports]
            logger.info(f"  Erowid: Found {len(reports)} reports")
        except Exception as e:
            logger.warning(f"  Erowid error: {e}")
        
        # Pharmacology data
        pharm = self.pharmacology.get_receptor_profile(substance_name)
        if pharm:
            data["pharmacology"] = {
                "bindings": [
                    {
                        "receptor": b.receptor,
                        "ki_nm": b.ki_nm,
                        "activity": b.activity.value
                    }
                    for b in pharm.bindings
                ],
                "primary_targets": pharm.primary_targets,
                "half_life_hours": pharm.half_life_hours
            }
            logger.info(f"  Pharmacology: {len(pharm.bindings)} receptor bindings")
        
        # Cache collected data
        cache_path = self.cache_dir / f"{substance_name.lower()}_data.json"
        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2, default=str)
        
        return data
    
    async def extract_effects_from_reports(
        self,
        reports: List[Dict],
        batch_size: int = 10
    ) -> List[Dict]:
        """
        Extract structured effect data from experience reports using LLM
        """
        logger.info(f"Extracting effects from {len(reports)} reports...")
        
        extracted = []
        
        for i in range(0, len(reports), batch_size):
            batch = reports[i:i+batch_size]
            
            for report in batch:
                text = report.get("text", "")
                if not text or len(text) < 100:
                    continue
                
                try:
                    result = await asyncio.to_thread(
                        self.effect_extractor.extract_effects,
                        text
                    )
                    extracted.append({
                        "report_id": report.get("url", f"report_{i}"),
                        "substance": report.get("substance"),
                        "dose": report.get("dose"),
                        "extracted": asdict(result)
                    })
                except Exception as e:
                    logger.warning(f"  Extraction error: {e}")
            
            logger.info(f"  Processed {min(i+batch_size, len(reports))}/{len(reports)}")
        
        return extracted
    
    def compute_effect_statistics(
        self,
        extracted_reports: List[Dict]
    ) -> Dict[str, Any]:
        """
        Compute statistical summary of effects across reports
        """
        effect_counts = defaultdict(int)
        effect_intensities = defaultdict(list)
        total_reports = len(extracted_reports)
        
        for report in extracted_reports:
            extracted = report.get("extracted", {})
            effects = extracted.get("effects", [])
            
            for effect in effects:
                name = effect.get("effect_name", "")
                if name:
                    effect_counts[name] += 1
                    if effect.get("intensity"):
                        effect_intensities[name].append(effect["intensity"])
        
        # Compute frequencies and average intensities
        stats = {}
        for effect, count in effect_counts.items():
            freq = count / total_reports if total_reports > 0 else 0
            avg_intensity = (
                sum(effect_intensities[effect]) / len(effect_intensities[effect])
                if effect_intensities[effect] else None
            )
            stats[effect] = {
                "frequency": freq,
                "count": count,
                "avg_intensity": avg_intensity
            }
        
        # Sort by frequency
        sorted_stats = dict(sorted(
            stats.items(),
            key=lambda x: x[1]["frequency"],
            reverse=True
        ))
        
        return {
            "total_reports": total_reports,
            "unique_effects": len(effect_counts),
            "effects": sorted_stats
        }
    
    def map_effects_to_brain_regions(
        self,
        effect_stats: Dict[str, Any]
    ) -> Dict[str, float]:
        """
        Map observed effects to brain regions using pharmacology data
        """
        region_scores = defaultdict(float)
        
        for effect_name, data in effect_stats.get("effects", {}).items():
            freq = data.get("frequency", 0)
            
            # Get brain regions for this effect from knowledge graph
            regions = self.knowledge_graph.get_brain_regions_for_effect(effect_name)
            
            for region_info in regions:
                region = region_info["region"]
                score = region_info["score"]
                region_scores[region] += freq * score
        
        # Normalize scores
        max_score = max(region_scores.values()) if region_scores else 1
        normalized = {
            k: v / max_score 
            for k, v in region_scores.items()
        }
        
        return dict(sorted(
            normalized.items(),
            key=lambda x: x[1],
            reverse=True
        ))
    
    def compare_predicted_vs_observed(
        self,
        substance: str,
        observed_effects: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Compare pharmacology-predicted effects with empirically observed effects
        """
        # Get predicted effects from pharmacology
        predicted = self.pharmacology.predict_effects_from_pharmacology(substance)
        predicted_set = {e["effect"] for e in predicted}
        
        # Get observed effects
        observed_set = set(observed_effects.get("effects", {}).keys())
        
        # Calculate overlap
        true_positives = predicted_set & observed_set
        false_positives = predicted_set - observed_set
        false_negatives = observed_set - predicted_set
        
        precision = len(true_positives) / len(predicted_set) if predicted_set else 0
        recall = len(true_positives) / len(observed_set) if observed_set else 0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0
        
        return {
            "predicted_effects": list(predicted_set),
            "observed_effects": list(observed_set),
            "correctly_predicted": list(true_positives),
            "overpredicted": list(false_positives),
            "underpredicted": list(false_negatives),
            "precision": precision,
            "recall": recall,
            "f1_score": f1
        }
    
    async def run_full_analysis(
        self,
        substance: str,
        max_reports: int = 50
    ) -> AnalysisResult:
        """
        Run complete analysis pipeline for a substance
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"FULL ANALYSIS: {substance}")
        logger.info(f"{'='*60}\n")
        
        # 1. Collect data
        data = await self.collect_substance_data(substance, max_reports)
        
        # 2. Extract effects from reports
        extracted = []
        if data["erowid_reports"]:
            extracted = await self.extract_effects_from_reports(
                data["erowid_reports"],
                batch_size=5
            )
        
        # 3. Compute effect statistics
        effect_stats = self.compute_effect_statistics(extracted)
        
        # 4. Map to brain regions
        brain_regions = self.map_effects_to_brain_regions(effect_stats)
        
        # 5. Compare predicted vs observed
        comparison = self.compare_predicted_vs_observed(substance, effect_stats)
        
        # 6. Get receptor profile
        receptor_profile = {}
        if data["pharmacology"]:
            receptor_profile = {
                b["receptor"]: {"ki_nm": b["ki_nm"], "activity": b["activity"]}
                for b in data["pharmacology"]["bindings"]
            }
        
        # 7. Build result
        result = AnalysisResult(
            substance=substance,
            n_reports_analyzed=len(extracted),
            effect_frequencies={
                k: v["frequency"] 
                for k, v in effect_stats.get("effects", {}).items()
            },
            predicted_effects=self.pharmacology.predict_effects_from_pharmacology(substance),
            brain_regions=brain_regions,
            receptor_profile=receptor_profile,
            effect_clusters=[],  # TODO: clustering analysis
            prediction_accuracy=comparison.get("f1_score")
        )
        
        # Save results
        output_path = self.output_dir / f"{substance.lower()}_analysis.json"
        with open(output_path, "w") as f:
            json.dump(asdict(result), f, indent=2, default=str)
        
        logger.info(f"\nResults saved to {output_path}")
        
        return result
    
    def generate_report(self, result: AnalysisResult) -> str:
        """Generate human-readable report from analysis"""
        
        report = f"""
{'='*60}
CONSCIOUSNESS RESEARCH REPORT: {result.substance}
{'='*60}

SUMMARY
-------
Reports Analyzed: {result.n_reports_analyzed}
Prediction Accuracy (F1): {result.prediction_accuracy:.2f if result.prediction_accuracy else 'N/A'}

RECEPTOR BINDING PROFILE
------------------------
"""
        for receptor, data in result.receptor_profile.items():
            report += f"  {receptor}: Ki={data['ki_nm']}nM ({data['activity']})\n"
        
        report += """
TOP 10 OBSERVED EFFECTS (by frequency)
--------------------------------------
"""
        top_effects = sorted(
            result.effect_frequencies.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        for effect, freq in top_effects:
            report += f"  {effect}: {freq*100:.1f}%\n"
        
        report += """
BRAIN REGIONS IMPLICATED
------------------------
"""
        top_regions = sorted(
            result.brain_regions.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        for region, score in top_regions:
            report += f"  {region}: {score:.2f}\n"
        
        report += """
MECHANISTIC INTERPRETATION
--------------------------
"""
        if result.receptor_profile:
            primary = list(result.receptor_profile.keys())[:3]
            report += f"Primary targets: {', '.join(primary)}\n"
            report += "Proposed mechanism: "
            
            if "5-HT2A" in result.receptor_profile:
                report += "Serotonergic psychedelia via 5-HT2A agonism in visual cortex and DMN disruption. "
            if "NMDA" in result.receptor_profile:
                report += "Glutamatergic dissociation via NMDA antagonism in cortex and hippocampus. "
            if "KOR" in result.receptor_profile:
                report += "Kappa-opioid mediated dysphoria and reality distortion via claustrum. "
        
        return report


async def main():
    """Example usage of the pipeline"""
    
    pipeline = ConsciousnessResearchPipeline(
        data_dir="./data",
        cache_dir="./cache",
        output_dir="./results"
    )
    
    # Analyze a substance
    substance = "LSD"
    
    try:
        result = await pipeline.run_full_analysis(
            substance=substance,
            max_reports=20  # Start small for testing
        )
        
        # Generate and print report
        report = pipeline.generate_report(result)
        print(report)
        
    except Exception as e:
        logger.error(f"Pipeline error: {e}")
        raise


if __name__ == "__main__":
    asyncio.run(main())
