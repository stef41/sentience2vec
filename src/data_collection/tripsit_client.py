"""
TripSit API Client
Fetches drug factsheets and combination data
"""

import json
from typing import Optional, Dict, List
import requests


class TripSitAPI:
    """Client for TripSit Factsheet API"""
    
    BASE_URL = "http://tripbot.tripsit.me/api/tripsit"
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json"
        })
    
    def get_drug(self, name: str) -> Optional[Dict]:
        """
        Fetch drug factsheet data
        
        Returns dict with:
        - name, pretty_name, aliases
        - categories
        - properties (dose, duration, onset, effects, etc.)
        - formatted_* fields
        - combos (interaction data)
        - links
        """
        # Use path-based endpoint (correct format)
        url = f"{self.BASE_URL}/getDrug/{name}"
        
        response = self.session.get(url)
        response.raise_for_status()
        
        data = response.json()
        
        if data.get("err"):
            return None
        
        return data.get("data", [{}])[0] if data.get("data") else None
    
    def get_all_drugs(self) -> List[Dict]:
        """Fetch all available drug factsheets"""
        url = f"{self.BASE_URL}/getAllDrugs"
        
        response = self.session.get(url)
        response.raise_for_status()
        
        data = response.json()
        return list(data.get("data", {}).values())
    
    def get_all_drug_names(self) -> List[str]:
        """Get list of all drug names in the database"""
        url = f"{self.BASE_URL}/getAllDrugNames"
        
        response = self.session.get(url)
        response.raise_for_status()
        
        data = response.json()
        return data.get("data", [])
    
    def get_all_drug_name_aliases(self) -> Dict[str, str]:
        """Get mapping of aliases to canonical drug names"""
        url = f"{self.BASE_URL}/getAllDrugNamesByAlias"
        
        response = self.session.get(url)
        response.raise_for_status()
        
        data = response.json()
        return data.get("data", {})
    
    def get_all_categories(self) -> List[str]:
        """Get all drug categories"""
        url = f"{self.BASE_URL}/getAllCategories"
        
        response = self.session.get(url)
        response.raise_for_status()
        
        data = response.json()
        return data.get("data", [])
    
    def get_interactions(self, drug_a: str, drug_b: str) -> Optional[Dict]:
        """
        Get interaction data between two drugs
        
        Status levels:
        - Low Risk & Synergy
        - Low Risk & No Synergy  
        - Low Risk & Decrease
        - Caution
        - Unsafe
        - Dangerous
        - Serotonin Syndrome
        """
        url = f"{self.BASE_URL}/getInteraction"
        params = {"drugA": drug_a, "drugB": drug_b}
        
        response = self.session.get(url, params=params)
        response.raise_for_status()
        
        data = response.json()
        
        if data.get("err"):
            return None
        
        return data.get("data", {})
    
    def parse_combo_data(self, drug_data: Dict) -> Dict[str, Dict]:
        """
        Parse combination/interaction data from drug factsheet
        
        Returns dict mapping drug names to interaction info:
        {
            "alcohol": {"status": "Dangerous", "note": "..."},
            "cannabis": {"status": "Low Risk & Synergy", "note": "..."},
            ...
        }
        """
        combos = {}
        combo_data = drug_data.get("combos", {})
        
        for drug_name, interaction in combo_data.items():
            status = interaction.get("status", "Unknown")
            note = interaction.get("note", "")
            combos[drug_name] = {
                "status": status,
                "note": note
            }
        
        return combos


def save_tripsit_data(api: TripSitAPI, output_dir: str = "data/raw/tripsit"):
    """Download and save all TripSit factsheet data"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("Fetching drug list...")
    drug_names = api.get_all_drug_names()
    print(f"Found {len(drug_names)} drugs")
    
    # Save name list
    with open(f"{output_dir}/drug_names.json", "w") as f:
        json.dump(drug_names, f, indent=2)
    
    # Save aliases
    print("Fetching aliases...")
    aliases = api.get_all_drug_name_aliases()
    with open(f"{output_dir}/drug_aliases.json", "w") as f:
        json.dump(aliases, f, indent=2)
    
    # Save categories
    print("Fetching categories...")
    categories = api.get_all_categories()
    with open(f"{output_dir}/categories.json", "w") as f:
        json.dump(categories, f, indent=2)
    
    # Fetch all detailed data
    print("Fetching all drug data...")
    all_drugs = api.get_all_drugs()
    with open(f"{output_dir}/all_drugs.json", "w") as f:
        json.dump(all_drugs, f, indent=2)
    
    # Save individual files
    for drug in all_drugs:
        name = drug.get("name", "unknown")
        safe_name = name.replace("/", "_").replace(" ", "_")
        with open(f"{output_dir}/drugs/{safe_name}.json", "w") as f:
            json.dump(drug, f, indent=2)
    
    print(f"Saved {len(all_drugs)} drugs to {output_dir}")
    return all_drugs


# Interaction matrix categories
INTERACTION_LEVELS = {
    "Low Risk & Synergy": 1,
    "Low Risk & No Synergy": 2,
    "Low Risk & Decrease": 3,
    "Caution": 4,
    "Unsafe": 5,
    "Dangerous": 6,
    "Serotonin Syndrome": 7
}


def build_interaction_matrix(api: TripSitAPI, drug_list: List[str]) -> Dict:
    """
    Build pairwise interaction matrix for a list of drugs
    
    Returns dict with:
    - matrix: 2D array of interaction levels
    - drugs: list of drug names (row/col labels)
    - interactions: detailed interaction info
    """
    import numpy as np
    
    n = len(drug_list)
    matrix = np.zeros((n, n), dtype=int)
    interactions = {}
    
    for i, drug_a in enumerate(drug_list):
        for j, drug_b in enumerate(drug_list):
            if i >= j:  # Skip diagonal and lower triangle
                continue
            
            interaction = api.get_interactions(drug_a, drug_b)
            if interaction:
                status = interaction.get("status", "Unknown")
                level = INTERACTION_LEVELS.get(status, 0)
                matrix[i, j] = level
                matrix[j, i] = level
                
                key = f"{drug_a}+{drug_b}"
                interactions[key] = interaction
    
    return {
        "matrix": matrix.tolist(),
        "drugs": drug_list,
        "interactions": interactions
    }


if __name__ == "__main__":
    api = TripSitAPI()
    
    # Example: Fetch LSD data
    lsd = api.get_drug("lsd")
    if lsd:
        print("=== LSD Factsheet ===")
        print(json.dumps(lsd, indent=2))
        
        print("\n=== Combinations ===")
        combos = api.parse_combo_data(lsd)
        for drug, info in combos.items():
            print(f"  {drug}: {info['status']}")
    
    # Example: Check specific interaction
    interaction = api.get_interactions("lsd", "mdma")
    if interaction:
        print(f"\n=== LSD + MDMA ===")
        print(json.dumps(interaction, indent=2))
    
    # Uncomment to download all data:
    # save_tripsit_data(api)
