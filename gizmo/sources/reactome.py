"""
Reactome REST API client.

License: Reactome data is CC BY 4.0 — https://reactome.org/license
API docs: https://reactome.org/ContentService/

Key endpoints used:
  /data/pathways/top/{species}         — top-level pathways
  /data/pathway/{id}/containedEvents   — reactions in a pathway
  /data/reaction/{id}/participants     — metabolites/participants
  /data/query/{id}                     — generic entity lookup
"""

from __future__ import annotations

import logging
from typing import Any, Optional
from urllib.parse import urljoin

import requests

from gizmo.schema import EdgeRole, MetaboliteNode, ReactionEdge, ReactionNode

log = logging.getLogger(__name__)

_BASE = "https://reactome.org/ContentService/"
_SPECIES_DEFAULT = "Homo sapiens"


class ReactomeClient:
    """Thin wrapper around the Reactome Content Service REST API."""

    def __init__(self, base_url: str = _BASE, timeout: int = 30) -> None:
        self.base_url = base_url
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})
        self.timeout = timeout

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------

    def _get(self, path: str, **params: Any) -> Any:
        url = urljoin(self.base_url, path)
        resp = self.session.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    # ------------------------------------------------------------------
    # Pathways
    # ------------------------------------------------------------------

    def top_pathways(self, species: str = _SPECIES_DEFAULT) -> list[dict]:
        """Return top-level pathway stubs for a species."""
        return self._get(f"data/pathways/top/{species}")

    def pathway_reactions(self, pathway_stid: str) -> list[dict]:
        """Return events (reactions) contained in a pathway."""
        return self._get(f"data/pathway/{pathway_stid}/containedEvents")

    # ------------------------------------------------------------------
    # Reactions → graph nodes/edges
    # ------------------------------------------------------------------

    def reaction_detail(self, reaction_stid: str) -> dict:
        return self._get(f"data/query/{reaction_stid}")

    def parse_reaction(
        self,
        detail: dict,
        pathway_stids: Optional[list[str]] = None,
    ) -> tuple[ReactionNode, list[MetaboliteNode], list[ReactionEdge]]:
        """
        Convert a Reactome reaction detail dict into graph primitives.

        Returns (ReactionNode, [MetaboliteNode, ...], [ReactionEdge, ...]).
        Metabolite nodes are compartment-aware; node IDs are
        "CHEBI:XXXXX@compartment" when both are available, else "reactome:{stid}".
        """
        rxn_stid = detail.get("stId", detail.get("dbId", "unknown"))
        rxn_id = f"reactome:{rxn_stid}"

        ec_list: list[str] = []
        for cat in detail.get("catalystActivity", []):
            act = cat.get("activity", {})
            if acc := act.get("accession"):
                ec_list.append(acc)

        rxn_node = ReactionNode(
            node_id=rxn_id,
            reactome_id=rxn_stid,
            name=detail.get("displayName", rxn_stid),
            reversible=detail.get("isReversible", False),
            ec_numbers=ec_list,
            pathways=pathway_stids or [],
            species=detail.get("species", [{}])[0].get("displayName") if detail.get("species") else None,
        )

        metabolites: list[MetaboliteNode] = []
        edges: list[ReactionEdge] = []

        def _parse_participants(entries: list[dict], role: EdgeRole) -> None:
            for entry in entries:
                met_node, edge = _participant_to_node_edge(entry, rxn_id, role)
                if met_node:
                    metabolites.append(met_node)
                    edges.append(edge)

        _parse_participants(detail.get("input", []), EdgeRole.SUBSTRATE)
        _parse_participants(detail.get("output", []), EdgeRole.PRODUCT)
        _parse_participants(detail.get("catalystActivity", []), EdgeRole.MODIFIER)

        return rxn_node, metabolites, edges


def _participant_to_node_edge(
    entry: dict, rxn_id: str, role: EdgeRole
) -> tuple[Optional[MetaboliteNode], Optional[ReactionEdge]]:
    """
    Map a Reactome participant/input/output dict to a MetaboliteNode + ReactionEdge.
    Returns (None, None) for non-small-molecule participants (complexes, proteins).
    """
    # Only process SimpleEntity (small molecules) and ChemicalDrug
    schema_class = entry.get("schemaClass", "")
    if schema_class not in {"SimpleEntity", "ChemicalDrug", "OtherEntity"}:
        return None, None

    stid = entry.get("stId") or str(entry.get("dbId", ""))
    chebi_xrefs = [
        ref.get("identifier")
        for ref in entry.get("crossReference", [])
        if ref.get("databaseName") == "ChEBI"
    ]
    chebi_id = f"CHEBI:{chebi_xrefs[0]}" if chebi_xrefs else None

    compartment_name: Optional[str] = None
    if comp := entry.get("compartment"):
        compartment_name = comp.get("displayName") or comp.get("name")

    node_id = chebi_id or f"reactome:{stid}"
    if compartment_name:
        node_id = f"{node_id}@{compartment_name}"

    met_node = MetaboliteNode(
        node_id=node_id,
        chebi_id=chebi_id,
        reactome_id=stid,
        name=entry.get("displayName", stid),
        compartment=compartment_name,
    )

    if role == EdgeRole.SUBSTRATE:
        src, tgt = node_id, rxn_id
    else:  # PRODUCT or MODIFIER — reaction → metabolite
        src, tgt = rxn_id, node_id

    edge = ReactionEdge(
        source=src,
        target=tgt,
        role=role,
        stoichiometry=entry.get("stoichiometry", 1.0),
        compartment=compartment_name,
    )

    return met_node, edge
