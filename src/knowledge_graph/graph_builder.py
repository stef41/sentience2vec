"""
Knowledge Graph Builder
Integrates substance, effect, receptor, and brain region data
"""

import json
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict
import networkx as nx
from enum import Enum


class NodeType(Enum):
    SUBSTANCE = "substance"
    RECEPTOR = "receptor"
    BRAIN_REGION = "brain_region"
    EFFECT = "effect"
    EFFECT_CATEGORY = "effect_category"
    EXPERIENCE_REPORT = "experience_report"
    MECHANISM = "mechanism"


class EdgeType(Enum):
    # Substance relationships
    BINDS_TO = "binds_to"                 # substance -> receptor
    PRODUCES_EFFECT = "produces_effect"   # substance -> effect
    METABOLIZES_TO = "metabolizes_to"     # substance -> substance
    
    # Receptor relationships
    LOCATED_IN = "located_in"             # receptor -> brain_region
    MODULATES = "modulates"               # receptor -> effect
    
    # Effect relationships
    BELONGS_TO_CATEGORY = "belongs_to_category"  # effect -> effect_category
    CO_OCCURS_WITH = "co_occurs_with"     # effect -> effect
    MEDIATED_BY = "mediated_by"           # effect -> mechanism
    
    # Brain region relationships
    CONNECTED_TO = "connected_to"         # brain_region -> brain_region
    PROCESSES = "processes"               # brain_region -> effect
    
    # Report relationships
    DESCRIBES_SUBSTANCE = "describes_substance"  # report -> substance
    REPORTS_EFFECT = "reports_effect"     # report -> effect


@dataclass
class GraphNode:
    id: str
    node_type: NodeType
    properties: Dict[str, Any]


@dataclass  
class GraphEdge:
    source: str
    target: str
    edge_type: EdgeType
    properties: Dict[str, Any]


class ConsciousnessKnowledgeGraph:
    """
    Multi-relational knowledge graph connecting:
    - Substances (drugs, compounds)
    - Receptors (molecular targets)
    - Brain regions (neuroanatomy)
    - Effects (subjective experiences)
    - Experience reports (empirical data)
    """
    
    def __init__(self):
        self.graph = nx.MultiDiGraph()
        self._effect_cooccurrence = defaultdict(lambda: defaultdict(int))
    
    def add_node(self, node: GraphNode) -> None:
        """Add a node to the graph"""
        self.graph.add_node(
            node.id,
            node_type=node.node_type.value,
            **node.properties
        )
    
    def add_edge(self, edge: GraphEdge) -> None:
        """Add an edge to the graph"""
        self.graph.add_edge(
            edge.source,
            edge.target,
            edge_type=edge.edge_type.value,
            **edge.properties
        )
    
    def add_substance(
        self,
        name: str,
        class_: Optional[str] = None,
        common_names: List[str] = None,
        pharmacology: Dict = None
    ) -> str:
        """Add a substance node"""
        node_id = f"substance:{name.lower().replace(' ', '_')}"
        self.add_node(GraphNode(
            id=node_id,
            node_type=NodeType.SUBSTANCE,
            properties={
                "name": name,
                "class": class_,
                "common_names": common_names or [],
                "pharmacology": pharmacology or {}
            }
        ))
        return node_id
    
    def add_receptor(
        self,
        name: str,
        receptor_type: Optional[str] = None,
        neurotransmitter: Optional[str] = None
    ) -> str:
        """Add a receptor node"""
        node_id = f"receptor:{name.lower().replace('-', '_')}"
        self.add_node(GraphNode(
            id=node_id,
            node_type=NodeType.RECEPTOR,
            properties={
                "name": name,
                "receptor_type": receptor_type,
                "neurotransmitter": neurotransmitter
            }
        ))
        return node_id
    
    def add_brain_region(
        self,
        name: str,
        parent_structure: Optional[str] = None,
        brodmann_areas: List[int] = None
    ) -> str:
        """Add a brain region node"""
        node_id = f"brain_region:{name.lower().replace(' ', '_')}"
        self.add_node(GraphNode(
            id=node_id,
            node_type=NodeType.BRAIN_REGION,
            properties={
                "name": name,
                "parent_structure": parent_structure,
                "brodmann_areas": brodmann_areas or []
            }
        ))
        return node_id
    
    def add_effect(
        self,
        name: str,
        category: Optional[str] = None,
        description: Optional[str] = None
    ) -> str:
        """Add an effect node"""
        node_id = f"effect:{name.lower().replace(' ', '_')}"
        self.add_node(GraphNode(
            id=node_id,
            node_type=NodeType.EFFECT,
            properties={
                "name": name,
                "category": category,
                "description": description
            }
        ))
        
        # Link to category if provided
        if category:
            cat_id = f"effect_category:{category.lower().replace(' ', '_')}"
            if cat_id not in self.graph:
                self.add_node(GraphNode(
                    id=cat_id,
                    node_type=NodeType.EFFECT_CATEGORY,
                    properties={"name": category}
                ))
            self.add_edge(GraphEdge(
                source=node_id,
                target=cat_id,
                edge_type=EdgeType.BELONGS_TO_CATEGORY,
                properties={}
            ))
        
        return node_id
    
    def add_binding(
        self,
        substance_id: str,
        receptor_id: str,
        ki_nm: Optional[float] = None,
        activity: Optional[str] = None
    ) -> None:
        """Add substance-receptor binding edge"""
        self.add_edge(GraphEdge(
            source=substance_id,
            target=receptor_id,
            edge_type=EdgeType.BINDS_TO,
            properties={
                "ki_nm": ki_nm,
                "activity": activity,
                "affinity_strength": self._ki_to_strength(ki_nm)
            }
        ))
    
    def add_receptor_localization(
        self,
        receptor_id: str,
        brain_region_id: str,
        density: str = "moderate"
    ) -> None:
        """Add receptor-brain region localization edge"""
        self.add_edge(GraphEdge(
            source=receptor_id,
            target=brain_region_id,
            edge_type=EdgeType.LOCATED_IN,
            properties={"density": density}
        ))
    
    def add_effect_relationship(
        self,
        receptor_id: str,
        effect_id: str,
        mechanism: Optional[str] = None,
        confidence: float = 0.5
    ) -> None:
        """Add receptor-effect modulation edge"""
        self.add_edge(GraphEdge(
            source=receptor_id,
            target=effect_id,
            edge_type=EdgeType.MODULATES,
            properties={
                "mechanism": mechanism,
                "confidence": confidence
            }
        ))
    
    def record_effect_cooccurrence(self, effects: List[str]) -> None:
        """Record co-occurrence of effects in the same experience"""
        for i, effect1 in enumerate(effects):
            for effect2 in effects[i+1:]:
                self._effect_cooccurrence[effect1][effect2] += 1
                self._effect_cooccurrence[effect2][effect1] += 1
    
    def build_cooccurrence_edges(self, min_count: int = 3) -> None:
        """Build co-occurrence edges from recorded data"""
        for effect1, cooccurs in self._effect_cooccurrence.items():
            effect1_id = f"effect:{effect1.lower().replace(' ', '_')}"
            if effect1_id not in self.graph:
                continue
                
            for effect2, count in cooccurs.items():
                if count >= min_count:
                    effect2_id = f"effect:{effect2.lower().replace(' ', '_')}"
                    if effect2_id not in self.graph:
                        continue
                    
                    self.add_edge(GraphEdge(
                        source=effect1_id,
                        target=effect2_id,
                        edge_type=EdgeType.CO_OCCURS_WITH,
                        properties={"count": count}
                    ))
    
    def _ki_to_strength(self, ki_nm: Optional[float]) -> str:
        """Convert Ki value to qualitative strength"""
        if ki_nm is None:
            return "unknown"
        if ki_nm < 10:
            return "very_high"
        if ki_nm < 100:
            return "high"
        if ki_nm < 1000:
            return "moderate"
        if ki_nm < 10000:
            return "low"
        return "very_low"
    
    # Query methods
    
    def get_substance_effects(self, substance_name: str) -> List[Dict]:
        """Get all effects associated with a substance"""
        substance_id = f"substance:{substance_name.lower().replace(' ', '_')}"
        
        if substance_id not in self.graph:
            return []
        
        effects = []
        
        # Direct substance -> effect edges
        for _, target, data in self.graph.out_edges(substance_id, data=True):
            if data.get("edge_type") == EdgeType.PRODUCES_EFFECT.value:
                effects.append({
                    "effect": target,
                    "relationship": "direct",
                    **data
                })
        
        # Indirect: substance -> receptor -> effect
        for _, receptor, bind_data in self.graph.out_edges(substance_id, data=True):
            if bind_data.get("edge_type") == EdgeType.BINDS_TO.value:
                for _, effect, mod_data in self.graph.out_edges(receptor, data=True):
                    if mod_data.get("edge_type") == EdgeType.MODULATES.value:
                        effects.append({
                            "effect": effect,
                            "relationship": "via_receptor",
                            "receptor": receptor,
                            "binding_affinity": bind_data.get("ki_nm"),
                            "mechanism_confidence": mod_data.get("confidence")
                        })
        
        return effects
    
    def get_effect_mechanisms(self, effect_name: str) -> List[Dict]:
        """Get mechanistic pathways for an effect"""
        effect_id = f"effect:{effect_name.lower().replace(' ', '_')}"
        
        if effect_id not in self.graph:
            return []
        
        mechanisms = []
        
        # Find receptors that modulate this effect
        for source, _, data in self.graph.in_edges(effect_id, data=True):
            if data.get("edge_type") == EdgeType.MODULATES.value:
                receptor_data = self.graph.nodes.get(source, {})
                
                # Find brain regions for this receptor
                regions = []
                for _, region, loc_data in self.graph.out_edges(source, data=True):
                    if loc_data.get("edge_type") == EdgeType.LOCATED_IN.value:
                        regions.append({
                            "region": region,
                            "density": loc_data.get("density")
                        })
                
                mechanisms.append({
                    "receptor": source,
                    "receptor_name": receptor_data.get("name"),
                    "mechanism": data.get("mechanism"),
                    "confidence": data.get("confidence"),
                    "brain_regions": regions
                })
        
        return mechanisms
    
    def find_similar_substances(self, substance_name: str, top_k: int = 5) -> List[Tuple[str, float]]:
        """Find substances with similar receptor binding profiles"""
        substance_id = f"substance:{substance_name.lower().replace(' ', '_')}"
        
        if substance_id not in self.graph:
            return []
        
        # Get receptors this substance binds to
        source_receptors = set()
        for _, target, data in self.graph.out_edges(substance_id, data=True):
            if data.get("edge_type") == EdgeType.BINDS_TO.value:
                source_receptors.add(target)
        
        if not source_receptors:
            return []
        
        # Find other substances and compute Jaccard similarity
        similarities = []
        
        for node in self.graph.nodes():
            if node.startswith("substance:") and node != substance_id:
                other_receptors = set()
                for _, target, data in self.graph.out_edges(node, data=True):
                    if data.get("edge_type") == EdgeType.BINDS_TO.value:
                        other_receptors.add(target)
                
                if other_receptors:
                    jaccard = len(source_receptors & other_receptors) / len(source_receptors | other_receptors)
                    similarities.append((node, jaccard))
        
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_k]
    
    def get_brain_regions_for_effect(self, effect_name: str) -> List[Dict]:
        """Get brain regions implicated in producing an effect"""
        mechanisms = self.get_effect_mechanisms(effect_name)
        
        region_scores = defaultdict(float)
        
        for mech in mechanisms:
            confidence = mech.get("confidence", 0.5)
            for region_info in mech.get("brain_regions", []):
                region = region_info["region"]
                density = region_info["density"]
                
                # Score based on confidence and density
                density_weight = {"high": 1.0, "moderate": 0.5, "low": 0.25}.get(density, 0.5)
                region_scores[region] += confidence * density_weight
        
        return sorted(
            [{"region": r, "score": s} for r, s in region_scores.items()],
            key=lambda x: x["score"],
            reverse=True
        )
    
    # Export methods
    
    def to_networkx(self) -> nx.MultiDiGraph:
        """Return the underlying NetworkX graph"""
        return self.graph
    
    def to_json(self) -> str:
        """Export graph to JSON"""
        data = {
            "nodes": [
                {"id": n, **self.graph.nodes[n]}
                for n in self.graph.nodes()
            ],
            "edges": [
                {"source": u, "target": v, **d}
                for u, v, d in self.graph.edges(data=True)
            ]
        }
        return json.dumps(data, indent=2)
    
    def to_cypher(self) -> List[str]:
        """Generate Cypher statements for Neo4j import"""
        statements = []
        
        # Create nodes
        for node_id in self.graph.nodes():
            node = self.graph.nodes[node_id]
            node_type = node.get("node_type", "Node")
            props = {k: v for k, v in node.items() if k != "node_type" and v is not None}
            props_str = json.dumps(props)
            statements.append(
                f"CREATE (n:{node_type} {{id: '{node_id}', {props_str[1:-1]}}})"
            )
        
        # Create edges
        for u, v, data in self.graph.edges(data=True):
            edge_type = data.get("edge_type", "RELATED_TO")
            props = {k: v for k, v in data.items() if k != "edge_type" and v is not None}
            props_str = json.dumps(props) if props else ""
            if props_str:
                statements.append(
                    f"MATCH (a {{id: '{u}'}}), (b {{id: '{v}'}}) CREATE (a)-[:{edge_type} {props_str}]->(b)"
                )
            else:
                statements.append(
                    f"MATCH (a {{id: '{u}'}}), (b {{id: '{v}'}}) CREATE (a)-[:{edge_type}]->(b)"
                )
        
        return statements
    
    def summary(self) -> Dict:
        """Get summary statistics of the graph"""
        node_types = defaultdict(int)
        edge_types = defaultdict(int)
        
        for node in self.graph.nodes():
            ntype = self.graph.nodes[node].get("node_type", "unknown")
            node_types[ntype] += 1
        
        for _, _, data in self.graph.edges(data=True):
            etype = data.get("edge_type", "unknown")
            edge_types[etype] += 1
        
        return {
            "total_nodes": self.graph.number_of_nodes(),
            "total_edges": self.graph.number_of_edges(),
            "node_types": dict(node_types),
            "edge_types": dict(edge_types)
        }


def build_base_graph() -> ConsciousnessKnowledgeGraph:
    """Build initial graph with core pharmacology data"""
    from consciousness_study.src.pharmacology.receptor_mapping import (
        RECEPTOR_BRAIN_REGIONS,
        SUBSTANCE_PHARMACOLOGY,
        EFFECT_MECHANISM_MAPPING
    )
    
    kg = ConsciousnessKnowledgeGraph()
    
    # Add receptors and their brain region localizations
    for receptor_name, info in RECEPTOR_BRAIN_REGIONS.items():
        receptor_id = kg.add_receptor(
            name=receptor_name,
            neurotransmitter=receptor_name.split("_")[0] if "_" in receptor_name else receptor_name[:2]
        )
        
        # Add brain regions and localization edges
        for region in info.get("high_density", []):
            region_id = kg.add_brain_region(name=region)
            kg.add_receptor_localization(receptor_id, region_id, density="high")
        
        for region in info.get("moderate_density", []):
            region_id = kg.add_brain_region(name=region)
            kg.add_receptor_localization(receptor_id, region_id, density="moderate")
        
        # Add effects associated with receptor
        for effect in info.get("associated_effects", []):
            effect_id = kg.add_effect(name=effect)
            kg.add_effect_relationship(receptor_id, effect_id, confidence=0.7)
    
    # Add substances and their bindings
    for name, pharm in SUBSTANCE_PHARMACOLOGY.items():
        substance_id = kg.add_substance(
            name=pharm.name,
            pharmacology={
                "half_life_hours": pharm.half_life_hours,
                "primary_targets": pharm.primary_targets
            }
        )
        
        for binding in pharm.bindings:
            receptor_id = f"receptor:{binding.receptor.lower().replace('-', '_')}"
            if receptor_id in kg.graph:
                kg.add_binding(
                    substance_id,
                    receptor_id,
                    ki_nm=binding.ki_nm,
                    activity=binding.activity.value
                )
    
    # Add effect mechanism mappings
    for effect_name, mech in EFFECT_MECHANISM_MAPPING.items():
        effect_id = kg.add_effect(name=effect_name)
        
        receptor_id = f"receptor:{mech['primary_receptor'].lower().replace('-', '_')}"
        if receptor_id in kg.graph:
            kg.add_effect_relationship(
                receptor_id,
                effect_id,
                mechanism=mech.get("mechanism"),
                confidence=mech.get("confidence", 0.5)
            )
        
        for region in mech.get("brain_regions", []):
            region_id = kg.add_brain_region(name=region)
            kg.add_edge(GraphEdge(
                source=region_id,
                target=effect_id,
                edge_type=EdgeType.PROCESSES,
                properties={"confidence": mech.get("confidence", 0.5)}
            ))
    
    return kg


if __name__ == "__main__":
    print("Building consciousness knowledge graph...")
    
    # Build base graph from pharmacology data
    kg = build_base_graph()
    
    print(f"\n=== Graph Summary ===")
    summary = kg.summary()
    print(f"Total nodes: {summary['total_nodes']}")
    print(f"Total edges: {summary['total_edges']}")
    print(f"\nNode types: {summary['node_types']}")
    print(f"Edge types: {summary['edge_types']}")
    
    print("\n=== Example Queries ===")
    
    print("\n--- Effects of LSD ---")
    effects = kg.get_substance_effects("LSD")
    for e in effects[:10]:
        print(f"  {e['effect']}: {e['relationship']}")
    
    print("\n--- Mechanisms of Ego Dissolution ---")
    mechanisms = kg.get_effect_mechanisms("ego_dissolution")
    for m in mechanisms:
        print(f"  {m['receptor_name']}: {m.get('mechanism', 'unknown')}")
        for r in m.get("brain_regions", []):
            print(f"    - {r['region']} ({r['density']})")
    
    print("\n--- Similar Substances to LSD ---")
    similar = kg.find_similar_substances("LSD")
    for name, score in similar:
        print(f"  {name}: {score:.2f}")
    
    print("\n--- Brain Regions for Hallucination ---")
    regions = kg.get_brain_regions_for_effect("hallucination")
    for r in regions:
        print(f"  {r['region']}: {r['score']:.2f}")
