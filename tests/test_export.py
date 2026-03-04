"""Tests for GraphML and JSON round-trip export/import."""

import tempfile
from pathlib import Path

import pytest

from gizmo.schema import EdgeRole, MetaboliteNode, ReactionEdge, ReactionNode
from gizmo.graph.network import MetaboliteGraph
from gizmo.export.graphml import write_graphml, read_graphml
from gizmo.export.json_export import write_json, read_json


@pytest.fixture()
def mg() -> MetaboliteGraph:
    g = MetaboliteGraph()
    g.add_metabolite(MetaboliteNode(node_id="CHEBI:1", name="A", ec_numbers=[], pathways=[]))
    g.add_reaction(ReactionNode(
        node_id="reactome:R1", name="R1", ec_numbers=["1.1.1.1"], pathways=["R-HSA-1"]
    ))
    g.add_edge(ReactionEdge(source="CHEBI:1", target="reactome:R1", role=EdgeRole.SUBSTRATE))
    return g


def test_json_roundtrip(mg):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = Path(f.name)
    try:
        write_json(mg, path)
        mg2 = read_json(path)
        assert mg2.summary() == mg.summary()
        assert "CHEBI:1" in mg2.metabolite_nodes()
    finally:
        path.unlink(missing_ok=True)


def test_graphml_roundtrip(mg):
    with tempfile.NamedTemporaryFile(suffix=".graphml", delete=False) as f:
        path = Path(f.name)
    try:
        write_graphml(mg, path)
        mg2 = read_graphml(path)
        assert mg2.summary() == mg.summary()
        assert "CHEBI:1" in mg2.metabolite_nodes()
    finally:
        path.unlink(missing_ok=True)
