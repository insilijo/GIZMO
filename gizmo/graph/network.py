"""
Core NetworkX-based metabolite graph.

Topology: directed graph with four node types:
  metabolite, reaction, disease, gene

Core bipartite edges:
  metabolite → reaction  (substrate/modifier)
  reaction   → metabolite (product)

Clinical overlay:
  disease → gene → reaction → metabolite
"""

from __future__ import annotations

from typing import Iterable

import networkx as nx

from gizmo.schema import (
    DiseaseEdge,
    DiseaseNode,
    EdgeRole,
    GeneNode,
    MetaboliteNode,
    ReactionEdge,
    ReactionNode,
)

_NODE_TYPES = {"metabolite", "reaction", "disease", "gene"}


class MetaboliteGraph:
    """Directed graph integrating metabolite reactions with clinical annotations."""

    def __init__(self) -> None:
        self._g: nx.DiGraph = nx.DiGraph()

    # ------------------------------------------------------------------
    # Node management
    # ------------------------------------------------------------------

    def add_metabolite(self, node: MetaboliteNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_reaction(self, node: ReactionNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_disease(self, node: DiseaseNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_gene(self, node: GeneNode) -> None:
        self._g.add_node(node.node_id, **node.model_dump())

    def add_metabolites(self, nodes: Iterable[MetaboliteNode]) -> None:
        for n in nodes:
            self.add_metabolite(n)

    def add_reactions(self, nodes: Iterable[ReactionNode]) -> None:
        for n in nodes:
            self.add_reaction(n)

    def add_diseases(self, nodes: Iterable[DiseaseNode]) -> None:
        for n in nodes:
            self.add_disease(n)

    def add_genes(self, nodes: Iterable[GeneNode]) -> None:
        for n in nodes:
            self.add_gene(n)

    # ------------------------------------------------------------------
    # Edge management
    # ------------------------------------------------------------------

    def add_edge(self, edge: ReactionEdge) -> None:
        self._g.add_edge(
            edge.source,
            edge.target,
            role=edge.role.value,
            stoichiometry=edge.stoichiometry,
            compartment=edge.compartment,
        )

    def add_disease_edge(self, edge: DiseaseEdge) -> None:
        self._g.add_edge(
            edge.source,
            edge.target,
            edge_type=edge.edge_type.value,
            score=edge.score,
            evidence_count=edge.evidence_count,
            source_db=edge.source_db,
        )

    def add_edges(self, edges: Iterable[ReactionEdge]) -> None:
        for e in edges:
            self.add_edge(e)

    def add_disease_edges(self, edges: Iterable[DiseaseEdge]) -> None:
        for e in edges:
            self.add_disease_edge(e)

    # ------------------------------------------------------------------
    # Typed node accessors
    # ------------------------------------------------------------------

    @property
    def graph(self) -> nx.DiGraph:
        return self._g

    def _nodes_of_type(self, node_type: str) -> list[str]:
        return [n for n, d in self._g.nodes(data=True) if d.get("node_type") == node_type]

    def metabolite_nodes(self) -> list[str]:
        return self._nodes_of_type("metabolite")

    def reaction_nodes(self) -> list[str]:
        return self._nodes_of_type("reaction")

    def disease_nodes(self) -> list[str]:
        return self._nodes_of_type("disease")

    def gene_nodes(self) -> list[str]:
        return self._nodes_of_type("gene")

    def currency_nodes(self) -> list[str]:
        """Metabolite nodes flagged as currency metabolites."""
        return [
            n
            for n, d in self._g.nodes(data=True)
            if d.get("node_type") == "metabolite" and d.get("is_currency", False)
        ]

    # ------------------------------------------------------------------
    # Traversal helpers
    # ------------------------------------------------------------------

    def neighbors_of_metabolite(self, node_id: str) -> list[str]:
        """Reaction nodes that involve this metabolite."""
        pred = [n for n in self._g.predecessors(node_id) if self._g.nodes[n].get("node_type") == "reaction"]
        succ = [n for n in self._g.successors(node_id) if self._g.nodes[n].get("node_type") == "reaction"]
        return list(set(pred + succ))

    def diseases_for_metabolite(self, node_id: str) -> list[str]:
        """Disease nodes directly linked to a metabolite (biomarker edges)."""
        return [
            n
            for n in self._g.predecessors(node_id)
            if self._g.nodes[n].get("node_type") == "disease"
        ]

    def diseases_for_reaction(self, node_id: str) -> list[str]:
        """Disease nodes linked to a reaction via pathway association."""
        return [
            n
            for n in self._g.predecessors(node_id)
            if self._g.nodes[n].get("node_type") == "disease"
        ]

    def metabolite_subgraph(self, node_ids: Iterable[str]) -> MetaboliteGraph:
        """
        Induced subgraph over a set of metabolite IDs, including reaction nodes
        that connect only those metabolites. Disease/gene nodes are excluded.
        """
        node_ids = set(node_ids)
        reaction_ids: set[str] = set()
        for rxn in self.reaction_nodes():
            rxn_mets = {
                n
                for n in list(self._g.predecessors(rxn)) + list(self._g.successors(rxn))
                if self._g.nodes[n].get("node_type") == "metabolite"
            }
            if rxn_mets and rxn_mets.issubset(node_ids):
                reaction_ids.add(rxn)

        sub_nx = self._g.subgraph(node_ids | reaction_ids).copy()
        mg = MetaboliteGraph()
        mg._g = sub_nx
        return mg

    # ------------------------------------------------------------------
    # In-place mutation helpers
    # ------------------------------------------------------------------

    def flag_currency(self, node_ids: Iterable[str]) -> int:
        """
        Mark nodes as currency metabolites in-place.
        Returns count of nodes flagged.
        """
        count = 0
        for nid in node_ids:
            if nid in self._g and self._g.nodes[nid].get("node_type") == "metabolite":
                self._g.nodes[nid]["is_currency"] = True
                count += 1
        return count

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def summary(self) -> dict:
        return {
            "metabolites": len(self.metabolite_nodes()),
            "reactions": len(self.reaction_nodes()),
            "diseases": len(self.disease_nodes()),
            "genes": len(self.gene_nodes()),
            "currency_flagged": len(self.currency_nodes()),
            "edges": self._g.number_of_edges(),
        }

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"MetaboliteGraph("
            f"{s['metabolites']} metabolites [{s['currency_flagged']} currency], "
            f"{s['reactions']} reactions, "
            f"{s['diseases']} diseases, "
            f"{s['genes']} genes, "
            f"{s['edges']} edges)"
        )
