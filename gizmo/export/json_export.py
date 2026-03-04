"""
JSON export/import using NetworkX node-link format.

This is the lossless round-trip format — all Python types are preserved.
The output is compatible with D3.js force-directed graphs and
can be loaded directly into PyTorch Geometric via custom loaders.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Union

import networkx as nx
from networkx.readwrite import json_graph

from gizmo.graph.network import MetaboliteGraph


def write_json(mg: MetaboliteGraph, path: Union[str, Path], indent: int = 2) -> None:
    """Serialise a MetaboliteGraph to JSON node-link format."""
    data = json_graph.node_link_data(mg.graph)
    Path(path).write_text(json.dumps(data, indent=indent, default=str))


def read_json(path: Union[str, Path]) -> MetaboliteGraph:
    """Load a MetaboliteGraph from JSON node-link format."""
    data = json.loads(Path(path).read_text())
    g = json_graph.node_link_graph(data, directed=True, multigraph=False)
    mg = MetaboliteGraph()
    mg._g = g
    return mg
