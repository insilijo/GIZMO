"""
Node and edge schemas for the GIZMO metabolite graph.

Graph topology: directed multi-typed graph.
  - MetaboliteNode  (node_type="metabolite")
  - ReactionNode    (node_type="reaction")
  - DiseaseNode     (node_type="disease")
  - GeneNode        (node_type="gene")

Core bipartite edges:
  metabolite → reaction  (substrate/modifier)
  reaction   → metabolite (product)

Clinical overlay edges (DiseaseEdge):
  disease → gene        (genetic association, Open Targets)
  disease → reaction    (pathway association, Reactome disease)
  disease → metabolite  (biomarker / causal metabolite, Orphanet)
  gene    → reaction    (enzyme-reaction, via EC/Reactome)

Stoichiometry on ReactionEdge; association scores on DiseaseEdge.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


# ---------------------------------------------------------------------------
# Edge role enums
# ---------------------------------------------------------------------------

class EdgeRole(str, Enum):
    SUBSTRATE = "substrate"
    PRODUCT = "product"
    MODIFIER = "modifier"    # catalytic / inhibitory


class DiseaseEdgeType(str, Enum):
    GENE_ASSOCIATED = "gene_associated"          # disease ↔ gene (Open Targets)
    PATHWAY_ASSOCIATED = "pathway_associated"    # disease → reaction/pathway
    BIOMARKER = "biomarker"                      # disease → metabolite (Orphanet)
    CAUSAL = "causal"                            # inborn error: enzyme defect → substrate accumulation
    GENE_REACTION = "gene_reaction"              # gene → reaction (enzyme catalysis)


# ---------------------------------------------------------------------------
# Metabolite node
# ---------------------------------------------------------------------------

class MetaboliteNode(BaseModel):
    """
    Represents a chemical species (compartment-aware).

    Canonical node_id: "CHEBI:XXXXX" or "CHEBI:XXXXX@compartment".
    Falls back to "reactome:{stId}" or "pubchem:{cid}" when ChEBI is unknown.
    """

    node_id: str
    chebi_id: Optional[str] = None
    metanetx_id: Optional[str] = None      # MNXM_XXXXXX
    reactome_id: Optional[str] = None      # R-ALL-XXXXXX PhysicalEntity stId
    pubchem_cid: Optional[str] = None
    metabolon_name: Optional[str] = None   # Metabolon BIOCHEMICAL field
    name: str
    formula: Optional[str] = None
    charge: Optional[int] = None
    mass: Optional[float] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    smiles: Optional[str] = None
    compartment: Optional[str] = None
    # Metabolon platform / chromatography context
    platform: Optional[str] = None        # e.g. "lc/ms neg"
    retention_index: Optional[float] = None
    # Graph analysis flags
    is_currency: bool = False              # ATP, NAD, H2O etc. — flag but don't delete
    node_type: str = "metabolite"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Reaction node
# ---------------------------------------------------------------------------

class ReactionNode(BaseModel):
    """Represents a biochemical reaction or transport event."""

    node_id: str                           # "reactome:R-HSA-XXXXXX" or "mnxr:MNXR_XXXXX"
    reactome_id: Optional[str] = None
    metanetx_id: Optional[str] = None     # MNXR_XXXXXX
    name: str
    reversible: bool = False
    direction: Optional[str] = None       # "left-to-right" | "right-to-left" | "bidirectional"
    ec_numbers: list[str] = Field(default_factory=list)
    gene_symbols: list[str] = Field(default_factory=list)  # HGNC symbols of catalysing genes
    pathways: list[str] = Field(default_factory=list)      # Reactome pathway stIDs
    species: Optional[str] = None
    node_type: str = "reaction"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Disease node
# ---------------------------------------------------------------------------

class DiseaseNode(BaseModel):
    """
    Represents a disease entity.

    Canonical node_id: "MONDO:XXXXXXX"
    Cross-references are stored as strings; we do not embed OMIM data.
    """

    node_id: str                           # "MONDO:XXXXXXX"
    mondo_id: Optional[str] = None
    name: str
    synonyms: list[str] = Field(default_factory=list)
    definition: Optional[str] = None
    # Cross-references (ID strings only — no embedded restricted data)
    xref_omim: list[str] = Field(default_factory=list)     # "OMIM:XXXXXX"
    xref_orphanet: list[str] = Field(default_factory=list) # "Orphanet:XXXXX"
    xref_doid: list[str] = Field(default_factory=list)     # "DOID:XXXXX"
    xref_icd10: list[str] = Field(default_factory=list)    # "ICD10:XXXXX"
    xref_mesh: list[str] = Field(default_factory=list)     # "MeSH:DXXXXXX"
    # Classification flags
    is_rare: bool = False
    is_inborn_error_of_metabolism: bool = False
    node_type: str = "disease"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Gene node
# ---------------------------------------------------------------------------

class GeneNode(BaseModel):
    """
    Represents a human gene.

    Canonical node_id: "ENSG:ENSGXXXXXXXXXXXX" (Ensembl)
    Falls back to "HGNC:XXXXX".
    """

    node_id: str
    ensembl_id: Optional[str] = None
    hgnc_id: Optional[str] = None
    symbol: str
    name: Optional[str] = None
    node_type: str = "gene"

    model_config = ConfigDict(frozen=True)


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------

class ReactionEdge(BaseModel):
    """Directed edge between a metabolite and a reaction node."""

    source: str
    target: str
    role: EdgeRole
    stoichiometry: float = 1.0
    compartment: Optional[str] = None

    model_config = ConfigDict(frozen=True)


class DiseaseEdge(BaseModel):
    """
    Directed edge connecting disease/gene nodes to metabolite/reaction nodes.
    score: Open Targets association score [0, 1] where available.
    """

    source: str
    target: str
    edge_type: DiseaseEdgeType
    score: Optional[float] = None
    evidence_count: Optional[int] = None
    source_db: Optional[str] = None       # "open_targets" | "orphanet" | "reactome"

    model_config = ConfigDict(frozen=True)
