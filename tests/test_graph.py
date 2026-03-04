"""Tests for core MetaboliteGraph construction and traversal."""

import pytest

from gizmo.schema import EdgeRole, MetaboliteNode, ReactionEdge, ReactionNode
from gizmo.graph.network import MetaboliteGraph


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

@pytest.fixture()
def simple_graph() -> MetaboliteGraph:
    """
    Minimal graph:
      glucose (substrate) → R_glycolysis → pyruvate (product)
                          → R_glycolysis → ATP (product)
      ATP (modifier)      → R_glycolysis
    """
    mg = MetaboliteGraph()

    glucose = MetaboliteNode(node_id="CHEBI:17234", chebi_id="CHEBI:17234", name="glucose")
    pyruvate = MetaboliteNode(node_id="CHEBI:15361", chebi_id="CHEBI:15361", name="pyruvate")
    atp = MetaboliteNode(
        node_id="CHEBI:30616@cytosol", chebi_id="CHEBI:30616", name="ATP", compartment="cytosol"
    )
    rxn = ReactionNode(
        node_id="reactome:R-HSA-111",
        reactome_id="R-HSA-111",
        name="Glycolysis step",
        reversible=False,
    )

    mg.add_metabolites([glucose, pyruvate, atp])
    mg.add_reaction(rxn)

    mg.add_edges([
        ReactionEdge(source="CHEBI:17234", target="reactome:R-HSA-111", role=EdgeRole.SUBSTRATE),
        ReactionEdge(source="reactome:R-HSA-111", target="CHEBI:15361", role=EdgeRole.PRODUCT),
        ReactionEdge(source="reactome:R-HSA-111", target="CHEBI:30616@cytosol", role=EdgeRole.PRODUCT),
    ])

    return mg


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------

def test_summary(simple_graph):
    s = simple_graph.summary()
    assert s["metabolites"] == 3
    assert s["reactions"] == 1
    assert s["edges"] == 3


def test_metabolite_nodes(simple_graph):
    mets = simple_graph.metabolite_nodes()
    assert "CHEBI:17234" in mets
    assert "reactome:R-HSA-111" not in mets


def test_reaction_nodes(simple_graph):
    rxns = simple_graph.reaction_nodes()
    assert "reactome:R-HSA-111" in rxns
    assert len(rxns) == 1


def test_neighbors_of_metabolite(simple_graph):
    neighbors = simple_graph.neighbors_of_metabolite("CHEBI:15361")
    assert "reactome:R-HSA-111" in neighbors


def test_repr(simple_graph):
    r = repr(simple_graph)
    assert "3 metabolites" in r
    assert "1 reactions" in r


def test_empty_graph():
    mg = MetaboliteGraph()
    s = mg.summary()
    assert s["metabolites"] == 0
    assert s["reactions"] == 0
    assert s["edges"] == 0
