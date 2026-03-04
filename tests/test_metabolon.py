"""Tests for Metabolon data dictionary loading and ChEBI mapping."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from gizmo.sources.metabolon import MetabolonLoader, MatchReport, _clean, _is_numeric


@pytest.fixture()
def sample_csv(tmp_path: Path) -> Path:
    """Create a minimal Metabolon-format CSV for testing."""
    rows = [
        {
            "BIOCHEMICAL": "glucose",
            "PATHWAY": "carbohydrate",
            "ChemRICHClass": "Saccharides",
            "PUBCHEM": "5793",
            "INCHIKEY": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
            "PLATFORM": "lc/ms pos early",
            "MASS": 180.0634,
            "RI": 1234,
        },
        {
            "BIOCHEMICAL": "ATP",
            "PATHWAY": "nucleotide",
            "ChemRICHClass": "Nucleotides",
            "PUBCHEM": "5957",
            "INCHIKEY": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N",
            "PLATFORM": "lc/ms neg",
            "MASS": 507.0,
            "RI": 987,
        },
        {
            "BIOCHEMICAL": "unknown_lipid_xyz",
            "PATHWAY": "lipid",
            "ChemRICHClass": "Fatty Acids",
            "PUBCHEM": "",
            "INCHIKEY": "",
            "PLATFORM": "lc/ms neg",
            "MASS": "",
            "RI": "",
        },
    ]
    df = pd.DataFrame(rows)
    out = tmp_path / "metabolon_test.csv"
    df.to_csv(out, index=False)
    return out


def test_load_raw(sample_csv):
    loader = MetabolonLoader(sample_csv)
    df = loader.raw_df()
    assert len(df) == 3


def test_to_nodes_no_index(sample_csv):
    """Without MetaNetX index, all matches should be unmatched."""
    loader = MetabolonLoader(sample_csv)
    nodes, report = loader.to_metabolite_nodes()
    assert report.total == 3
    assert report.exact_inchikey == 0
    assert report.unmatched == 3
    assert len(nodes) == 3


def test_node_fields(sample_csv):
    loader = MetabolonLoader(sample_csv)
    nodes, _ = loader.to_metabolite_nodes()
    glucose_node = next(n for n in nodes if n.metabolon_name == "glucose")
    assert glucose_node.inchikey == "WQZGKKKJIJFFOK-GASJEMHNSA-N"
    assert glucose_node.pubchem_cid == "5793"
    assert glucose_node.mass == 180.0634


def test_unknown_fallback_id(sample_csv):
    loader = MetabolonLoader(sample_csv)
    nodes, _ = loader.to_metabolite_nodes()
    unknown = next(n for n in nodes if n.metabolon_name == "unknown_lipid_xyz")
    assert unknown.node_id.startswith("metabolon:")


def test_metanetx_index_enrichment(sample_csv):
    """Simulate MetaNetX index with known InChIKey → ChEBI mapping."""
    loader = MetabolonLoader(sample_csv)
    # Manually inject the index to simulate a loaded MetaNetX file
    loader._inchikey_to_chebi["WQZGKKKJIJFFOK-GASJEMHNSA-N"] = "CHEBI:17234"  # glucose
    loader._inchikey_to_chebi["ZKHQWZAMYRWXGA-KQYNXXCUSA-N"] = "CHEBI:30616"  # ATP

    nodes, report = loader.to_metabolite_nodes()

    assert report.exact_inchikey == 2
    assert report.unmatched == 1
    assert report.chebi_coverage == pytest.approx(2 / 3)

    glucose_node = next(n for n in nodes if n.metabolon_name == "glucose")
    assert glucose_node.chebi_id == "CHEBI:17234"
    assert glucose_node.node_id == "CHEBI:17234"


def test_match_report_str():
    r = MatchReport(total=100, exact_inchikey=70, pubchem_inchikey=10, unmatched=20)
    s = str(r)
    assert "80.0%" in s   # chebi_coverage = 80/100


def test_clean_helper():
    assert _clean("") is None
    assert _clean(float("nan")) is None
    assert _clean("ABCDEF") == "ABCDEF"


def test_is_numeric():
    assert _is_numeric("123.4") is True
    assert _is_numeric("") is False
    assert _is_numeric(None) is False
