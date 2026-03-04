"""
Currency metabolite detection and flagging.

Currency metabolites are ubiquitous cofactors/co-substrates that participate
in hundreds of reactions. Including them in metabolite–metabolite co-occurrence
graphs creates spurious edges (everything is connected via ATP/NAD+).

Strategy:
  1. Flag a canonical known set by ChEBI ID.
  2. Optionally flag by degree heuristic: metabolites with reaction-degree
     > mean + k*std (default k=3) are likely unlisted currencies.
  3. Flagged nodes are NOT deleted — stoichiometry data is preserved for FBA.
     Graph kernels / co-occurrence projections should filter on is_currency=False.

ChEBI IDs sourced from Recon3D / published currency lists:
  Brunk et al. 2018, Thiele et al. 2013, Hucka et al. 2019.
"""

from __future__ import annotations

import statistics
from typing import Optional

from gizmo.graph.network import MetaboliteGraph

# ---------------------------------------------------------------------------
# Canonical currency metabolite ChEBI IDs
# (compartment-agnostic — match by chebi_id attribute, not node_id)
# ---------------------------------------------------------------------------

KNOWN_CURRENCY_CHEBI: dict[str, str] = {
    # Adenine nucleotides
    "CHEBI:30616": "ATP",
    "CHEBI:456216": "ADP",
    "CHEBI:16027": "AMP",
    "CHEBI:17877": "adenosine",
    # Guanine nucleotides
    "CHEBI:37565": "GTP",
    "CHEBI:58189": "GDP",
    "CHEBI:5978":  "GMP",
    # Pyridine nucleotides — oxidised
    "CHEBI:57540": "NAD+",
    "CHEBI:58349": "NADP+",
    # Pyridine nucleotides — reduced
    "CHEBI:57945": "NADH",
    "CHEBI:57783": "NADPH",
    # Flavin nucleotides
    "CHEBI:57692": "FAD",
    "CHEBI:58307": "FADH2",
    # Coenzyme A and acyl carriers
    "CHEBI:57287": "CoA",
    "CHEBI:57288": "acetyl-CoA",
    # Phosphates
    "CHEBI:43474": "orthophosphate (Pi)",
    "CHEBI:33019": "pyrophosphate (PPi)",
    # Small inorganic / proton carriers
    "CHEBI:15377": "water",
    "CHEBI:15378": "H+ (proton)",
    "CHEBI:16526": "CO2",
    "CHEBI:15379": "O2",
    "CHEBI:29033": "H2",
    # Ubiquinones
    "CHEBI:16389": "ubiquinone (CoQ)",
    "CHEBI:17976": "ubiquinol (CoQH2)",
    # Glutamate / 2-OG (transamination hub)
    "CHEBI:16015": "L-glutamate",
    "CHEBI:16810": "2-oxoglutarate",
    # Thioredoxin system (small, but extremely high degree)
    "CHEBI:15422": "glucose-6-phosphate (sometimes)",   # not always currency; flag for review
}

# ChEBI IDs that are borderline — flag but mark as "review"
BORDERLINE_CURRENCY_CHEBI: set[str] = {
    "CHEBI:15422",   # glucose-6-phosphate (central metabolite, not always excluded)
    "CHEBI:16015",   # L-glutamate
    "CHEBI:16810",   # 2-oxoglutarate
    "CHEBI:17877",   # adenosine
}


def flag_currency_metabolites(
    mg: MetaboliteGraph,
    *,
    degree_threshold_k: Optional[float] = 3.0,
    include_borderline: bool = False,
) -> dict[str, list[str]]:
    """
    Flag currency metabolites in-place on mg.

    Parameters
    ----------
    mg : MetaboliteGraph
        Graph to annotate.
    degree_threshold_k : float | None
        Flag metabolites whose reaction-degree > mean + k*std as statistical
        currency. Set to None to skip heuristic detection.
    include_borderline : bool
        Whether to also flag borderline currency metabolites.

    Returns
    -------
    dict with keys:
      "canonical"   — node_ids flagged by known ChEBI list
      "statistical" — node_ids flagged by degree heuristic
      "total"       — all flagged node_ids
    """
    g = mg.graph
    canonical: list[str] = []
    statistical: list[str] = []

    currency_chebi = set(KNOWN_CURRENCY_CHEBI.keys())
    if not include_borderline:
        currency_chebi -= BORDERLINE_CURRENCY_CHEBI

    # --- Pass 1: canonical ChEBI matching ---
    for nid, data in g.nodes(data=True):
        if data.get("node_type") != "metabolite":
            continue
        chebi = data.get("chebi_id")
        if chebi and chebi in currency_chebi:
            g.nodes[nid]["is_currency"] = True
            canonical.append(nid)

    # --- Pass 2: degree heuristic ---
    if degree_threshold_k is not None:
        met_nodes = mg.metabolite_nodes()
        # Degree = number of distinct reaction neighbours
        degrees = []
        for nid in met_nodes:
            rxn_neighbors = {
                n for n in list(g.predecessors(nid)) + list(g.successors(nid))
                if g.nodes[n].get("node_type") == "reaction"
            }
            degrees.append((nid, len(rxn_neighbors)))

        if len(degrees) > 2:
            vals = [d for _, d in degrees]
            mean = statistics.mean(vals)
            stdev = statistics.stdev(vals) if len(vals) > 1 else 0.0
            cutoff = mean + degree_threshold_k * stdev

            for nid, deg in degrees:
                if deg > cutoff and not g.nodes[nid].get("is_currency", False):
                    g.nodes[nid]["is_currency"] = True
                    statistical.append(nid)

    all_flagged = list(set(canonical + statistical))
    return {
        "canonical": canonical,
        "statistical": statistical,
        "total": all_flagged,
    }


def noncurrency_subgraph(mg: MetaboliteGraph) -> MetaboliteGraph:
    """
    Return a MetaboliteGraph with all currency-flagged metabolites removed.
    Reactions that become substrate-less or product-less are also removed.
    Use for metabolite co-occurrence / graph kernel analysis.
    """
    g = mg.graph
    keep_mets = {
        n for n in mg.metabolite_nodes()
        if not g.nodes[n].get("is_currency", False)
    }
    keep_rxns: set[str] = set()
    for rxn in mg.reaction_nodes():
        substrates = {
            n for n in g.predecessors(rxn)
            if g.nodes[n].get("node_type") == "metabolite"
        }
        products = {
            n for n in g.successors(rxn)
            if g.nodes[n].get("node_type") == "metabolite"
        }
        # Keep reaction only if it still has at least one non-currency participant on each side
        if substrates & keep_mets and products & keep_mets:
            keep_rxns.add(rxn)

    from gizmo.graph.network import MetaboliteGraph as MG
    sub = MG()
    sub._g = g.subgraph(keep_mets | keep_rxns).copy()
    return sub
