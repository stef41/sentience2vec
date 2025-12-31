"""
PsychonautWiki API Client
Fetches substance data via GraphQL API
"""

import json
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, asdict
import requests


@dataclass
class DoseRange:
    min: Optional[float] = None
    max: Optional[float] = None


@dataclass
class Dosage:
    units: Optional[str] = None
    threshold: Optional[float] = None
    light: Optional[DoseRange] = None
    common: Optional[DoseRange] = None
    strong: Optional[DoseRange] = None
    heavy: Optional[float] = None


@dataclass
class RouteOfAdministration:
    name: str
    dose: Optional[Dosage] = None
    duration: Optional[Dict] = None
    bioavailability: Optional[Dict] = None


@dataclass  
class SubjectiveEffect:
    name: str
    url: Optional[str] = None


@dataclass
class Substance:
    name: str
    common_names: List[str]
    chemical_class: Optional[str]
    psychoactive_class: Optional[str]
    roas: List[RouteOfAdministration]
    effects: List[SubjectiveEffect]
    tolerance: Optional[Dict]
    cross_tolerances: List[str]
    addictionPotential: Optional[str]
    toxicity: List[str]
    uncertainInteractions: List[str]
    unsafeInteractions: List[str]
    dangerousInteractions: List[str]


class PsychonautWikiAPI:
    """Client for PsychonautWiki GraphQL API"""
    
    ENDPOINT = "https://api.psychonautwiki.org/"
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json"
        })
    
    def _query(self, query: str, variables: Optional[Dict] = None) -> Dict:
        """Execute GraphQL query"""
        payload = {"query": query}
        if variables:
            payload["variables"] = variables
        
        response = self.session.post(self.ENDPOINT, json=payload)
        response.raise_for_status()
        
        result = response.json()
        if "errors" in result:
            raise Exception(f"GraphQL errors: {result['errors']}")
        
        return result.get("data", {})
    
    def get_substance(self, name: str) -> Optional[Dict]:
        """Fetch detailed substance information"""
        query = """
        query GetSubstance($name: String!) {
            substances(query: $name) {
                name
                commonNames
                class {
                    chemical
                    psychoactive
                }
                roas {
                    name
                    dose {
                        units
                        threshold
                        light { min max }
                        common { min max }
                        strong { min max }
                        heavy
                    }
                    duration {
                        onset { min max units }
                        comeup { min max units }
                        peak { min max units }
                        offset { min max units }
                        afterglow { min max units }
                        total { min max units }
                    }
                    bioavailability { min max }
                }
                effects {
                    name
                    url
                }
                tolerance {
                    full
                    half
                    zero
                }
                crossTolerances
                addictionPotential
                toxicity
                uncertainInteractions {
                    name
                }
                unsafeInteractions {
                    name
                }
                dangerousInteractions {
                    name
                }
            }
        }
        """
        
        data = self._query(query, {"name": name})
        substances = data.get("substances", [])
        return substances[0] if substances else None
    
    def get_all_substances(self) -> List[Dict]:
        """Fetch list of all substances"""
        query = """
        {
            substances(limit: 500) {
                name
                commonNames
                class {
                    chemical
                    psychoactive
                }
            }
        }
        """
        
        data = self._query(query)
        return data.get("substances", [])
    
    def get_effects_for_substance(self, name: str) -> List[Dict]:
        """Get effects associated with a substance"""
        query = """
        query GetEffects($name: String!) {
            substances(query: $name) {
                name
                effects {
                    name
                    url
                }
            }
        }
        """
        
        data = self._query(query, {"name": name})
        substances = data.get("substances", [])
        if substances:
            return substances[0].get("effects", [])
        return []
    
    def search_substances_by_effect(self, effect_name: str) -> List[str]:
        """Find substances that produce a specific effect"""
        # Note: This requires fetching all substances and filtering
        # Could be optimized with server-side filtering if API supports it
        all_substances = self.get_all_substances()
        matching = []
        
        for sub in all_substances:
            full_data = self.get_substance(sub["name"])
            if full_data:
                effects = [e["name"].lower() for e in full_data.get("effects", [])]
                if effect_name.lower() in " ".join(effects):
                    matching.append(sub["name"])
        
        return matching


def save_substance_data(api: PsychonautWikiAPI, output_dir: str = "data/raw/psychonautwiki"):
    """Download and save all substance data"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("Fetching substance list...")
    substances = api.get_all_substances()
    print(f"Found {len(substances)} substances")
    
    # Save summary list
    with open(f"{output_dir}/substances_list.json", "w") as f:
        json.dump(substances, f, indent=2)
    
    # Fetch detailed data for each
    detailed_data = []
    for i, sub in enumerate(substances):
        name = sub["name"]
        print(f"[{i+1}/{len(substances)}] Fetching: {name}")
        
        try:
            full_data = api.get_substance(name)
            if full_data:
                detailed_data.append(full_data)
                
                # Save individual file
                safe_name = name.replace("/", "_").replace(" ", "_")
                with open(f"{output_dir}/{safe_name}.json", "w") as f:
                    json.dump(full_data, f, indent=2)
        except Exception as e:
            print(f"  Error: {e}")
    
    # Save combined dataset
    with open(f"{output_dir}/all_substances_detailed.json", "w") as f:
        json.dump(detailed_data, f, indent=2)
    
    print(f"Saved {len(detailed_data)} substances to {output_dir}")
    return detailed_data


if __name__ == "__main__":
    api = PsychonautWikiAPI()
    
    # Example: Fetch LSD data
    lsd = api.get_substance("LSD")
    if lsd:
        print(json.dumps(lsd, indent=2))
    
    # Uncomment to download all data:
    # save_substance_data(api)
