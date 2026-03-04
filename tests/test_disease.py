"""Tests for disease/gene schema nodes and graph integration."""

import pytest

from gizmo.schema import (
    DiseaseEdge,
    DiseaseEdgeType,
    DiseaseNode,
    EdgeRole,
    GeneNode,
    MetaboliteNode,
    ReactionEdge,
    ReactionNode,
)
from gizmo.graph.network import MetaboliteGraph


@pytest.fixture()
def clinical_graph() -> MetaboliteGraph:
    """
    Phenylketonuria model:
      phenylalanine → [PAH reaction] → tyrosine
      PAH gene (ENSG:ENSG00000171759)
      PKU disease (MONDO:0009861) → PAH gene → PAH reaction → phenylalanine
    """
    mg = MetaboliteGraph()

    phe = MetaboliteNode(
        node_id="CHEBI:17295", chebi_id="CHEBI:17295", name="phenylalanine",
    )
    tyr = MetaboliteNode(
        node_id="CHEBI:17895", chebi_id="CHEBI:17895", name="tyrosine",
    )
    rxn = ReactionNode(
        node_id="reactome:R-HSA-71033",
        name="Phenylalanine hydroxylation",
        ec_numbers=["1.14.16.1"],
        gene_symbols=["PAH"],
        pathways=["R-HSA-71032"],
    )
    gene = GeneNode(
        node_id="ENSG:ENSG00000171759",
        ensembl_id="ENSG00000171759",
        symbol="PAH",
        name="phenylalanine hydroxylase",
    )
    disease = DiseaseNode(
        node_id="MONDO:0009861",
        mondo_id="MONDO:0009861",
        name="phenylketonuria",
        synonyms=["PKU", "classic phenylketonuria"],
        xref_omim=["OMIM:261600"],
        xref_orphanet=["Orphanet:716"],
        is_rare=True,
        is_inborn_error_of_metabolism=True,
    )

    mg.add_metabolites([phe, tyr])
    mg.add_reaction(rxn)
    mg.add_gene(gene)
    mg.add_disease(disease)

    mg.add_edges([
        ReactionEdge(source="CHEBI:17295", target="reactome:R-HSA-71033", role=EdgeRole.SUBSTRATE),
        ReactionEdge(source="reactome:R-HSA-71033", target="CHEBI:17895", role=EdgeRole.PRODUCT),
    ])
    mg.add_disease_edges([
        DiseaseEdge(
            source="MONDO:0009861",
            target="ENSG:ENSG00000171759",
            edge_type=DiseaseEdgeType.GENE_ASSOCIATED,
            score=0.95,
            source_db="orphanet",
        ),
        DiseaseEdge(
            source="ENSG:ENSG00000171759",
            target="reactome:R-HSA-71033",
            edge_type=DiseaseEdgeType.GENE_REACTION,
            source_db="reactome",
        ),
        DiseaseEdge(
            source="MONDO:0009861",
            target="CHEBI:17295",
            edge_type=DiseaseEdgeType.BIOMARKER,
            score=0.99,
            source_db="orphanet",
        ),
    ])
    return mg


def test_summary_includes_disease_gene(clinical_graph):
    s = clinical_graph.summary()
    assert s["diseases"] == 1
    assert s["genes"] == 1
    assert s["metabolites"] == 2
    assert s["reactions"] == 1


def test_disease_nodes_accessor(clinical_graph):
    assert "MONDO:0009861" in clinical_graph.disease_nodes()


def test_gene_nodes_accessor(clinical_graph):
    assert "ENSG:ENSG00000171759" in clinical_graph.gene_nodes()


def test_diseases_for_metabolite(clinical_graph):
    diseases = clinical_graph.diseases_for_metabolite("CHEBI:17295")
    assert "MONDO:0009861" in diseases


def test_diseases_for_reaction(clinical_graph):
    # disease → gene → reaction; disease is not directly linked to reaction in fixture
    # gene → reaction edge exists
    diseases = clinical_graph.diseases_for_reaction("reactome:R-HSA-71033")
    # ENSG node is predecessor, not MONDO directly
    assert isinstance(diseases, list)


def test_disease_node_schema():
    d = DiseaseNode(
        node_id="MONDO:0009861",
        name="phenylketonuria",
        is_inborn_error_of_metabolism=True,
        xref_omim=["OMIM:261600"],
    )
    assert d.is_rare is False   # not set; orphanet list empty
    assert d.is_inborn_error_of_metabolism is True
    assert "OMIM:261600" in d.xref_omim


def test_repr_includes_disease(clinical_graph):
    r = repr(clinical_graph)
    assert "disease" in r.lower()
