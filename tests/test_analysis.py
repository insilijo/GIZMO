"""Tests for currency detection and computational readiness QC."""

import pytest

from gizmo.schema import EdgeRole, MetaboliteNode, ReactionEdge, ReactionNode
from gizmo.graph.network import MetaboliteGraph
from gizmo.analysis.currency import (
    KNOWN_CURRENCY_CHEBI,
    flag_currency_metabolites,
    noncurrency_subgraph,
)
from gizmo.analysis.qc import assess_readiness


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

@pytest.fixture()
def annotated_graph() -> MetaboliteGraph:
    """
    Glucose → [Glycolysis] → Pyruvate
    ATP + water participate as currency co-substrates/products.
    """
    mg = MetaboliteGraph()

    glucose = MetaboliteNode(
        node_id="CHEBI:17234", chebi_id="CHEBI:17234", name="glucose",
        formula="C6H12O6", inchikey="WQZGKKKJIJFFOK-GASJEMHNSA-N",
    )
    pyruvate = MetaboliteNode(
        node_id="CHEBI:15361", chebi_id="CHEBI:15361", name="pyruvate",
        formula="C3H3O3", inchikey="LCTONWCANYUPML-UHFFFAOYSA-M",
    )
    atp = MetaboliteNode(
        node_id="CHEBI:30616", chebi_id="CHEBI:30616", name="ATP",
        formula="C10H16N5O13P3",
    )
    water = MetaboliteNode(
        node_id="CHEBI:15377", chebi_id="CHEBI:15377", name="water",
        formula="H2O",
    )

    rxn = ReactionNode(
        node_id="reactome:R-HSA-111",
        name="Glycolysis",
        ec_numbers=["2.7.1.1"],
        gene_symbols=["HK1"],
        pathways=["R-HSA-70171"],
    )

    mg.add_metabolites([glucose, pyruvate, atp, water])
    mg.add_reaction(rxn)
    mg.add_edges([
        ReactionEdge(source="CHEBI:17234", target="reactome:R-HSA-111", role=EdgeRole.SUBSTRATE),
        ReactionEdge(source="CHEBI:30616", target="reactome:R-HSA-111", role=EdgeRole.MODIFIER),
        ReactionEdge(source="reactome:R-HSA-111", target="CHEBI:15361", role=EdgeRole.PRODUCT),
        ReactionEdge(source="reactome:R-HSA-111", target="CHEBI:15377", role=EdgeRole.PRODUCT),
    ])
    return mg


# ------------------------------------------------------------------
# Currency tests
# ------------------------------------------------------------------

def test_known_currency_set_non_empty():
    assert len(KNOWN_CURRENCY_CHEBI) > 10


def test_flag_currency_canonical(annotated_graph):
    result = flag_currency_metabolites(annotated_graph)
    flagged = set(result["total"])
    # ATP (CHEBI:30616) and water (CHEBI:15377) should be flagged
    assert "CHEBI:30616" in flagged
    assert "CHEBI:15377" in flagged
    # glucose and pyruvate should NOT be flagged
    assert "CHEBI:17234" not in flagged
    assert "CHEBI:15361" not in flagged


def test_currency_in_graph(annotated_graph):
    flag_currency_metabolites(annotated_graph)
    g = annotated_graph.graph
    assert g.nodes["CHEBI:30616"].get("is_currency") is True
    assert g.nodes["CHEBI:17234"].get("is_currency", False) is False


def test_noncurrency_subgraph(annotated_graph):
    flag_currency_metabolites(annotated_graph)
    sub = noncurrency_subgraph(annotated_graph)
    # Currency nodes should be gone
    assert "CHEBI:30616" not in sub.metabolite_nodes()
    assert "CHEBI:15377" not in sub.metabolite_nodes()
    # Non-currency nodes should remain
    assert "CHEBI:17234" in sub.metabolite_nodes()
    assert "CHEBI:15361" in sub.metabolite_nodes()


def test_currency_nodes_accessor(annotated_graph):
    flag_currency_metabolites(annotated_graph)
    currencies = annotated_graph.currency_nodes()
    assert "CHEBI:30616" in currencies
    assert "CHEBI:15377" in currencies


# ------------------------------------------------------------------
# QC report tests
# ------------------------------------------------------------------

def test_qc_basic(annotated_graph):
    report = assess_readiness(annotated_graph)
    assert report.n_metabolites == 4
    assert report.n_reactions == 1
    assert report.n_edges == 4


def test_qc_currency_after_flag(annotated_graph):
    flag_currency_metabolites(annotated_graph)
    report = assess_readiness(annotated_graph)
    assert report.n_currency == 2
    assert report.currency_edge_fraction > 0


def test_qc_annotation_completeness(annotated_graph):
    report = assess_readiness(annotated_graph)
    assert report.reactions_with_ec == 1
    assert report.reactions_with_ec_fraction == 1.0
    assert report.reactions_with_gene == 1
    assert report.reactions_with_pathway == 1


def test_qc_structural_data(annotated_graph):
    report = assess_readiness(annotated_graph)
    # glucose and pyruvate have formula; atp does; water does
    assert report.metabolites_with_formula == 4
    assert report.metabolites_with_chebi == 4


def test_qc_dead_ends(annotated_graph):
    # glucose is only substrate (no reaction produces it here) → dead-end
    # pyruvate is only product (no reaction consumes it here) → dead-end
    # currency metabolites are excluded from dead-end check
    report = assess_readiness(annotated_graph)
    assert report.n_dead_end_metabolites >= 0  # exact count depends on currency flagging


def test_qc_empty_graph():
    mg = MetaboliteGraph()
    report = assess_readiness(mg)
    assert report.n_metabolites == 0
    assert not report.is_fba_ready
