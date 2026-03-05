"""
MetaNetX FTP/REST client for cross-referencing and ID mapping.

License: MetaNetX data is CC BY 4.0 — https://www.metanetx.org/mnxdoc/mnxref.html

MetaNetX provides flat TSV files via FTP and a REST API.
Key files used:
  chem_xref.tsv   — maps external IDs (ChEBI, KEGG, HMDB, …) → MNX IDs
  reac_xref.tsv   — maps external reaction IDs → MNXR IDs
  chem_prop.tsv   — name, formula, charge, InChIKey for each MNX compound
  reac_prop.tsv   — stoichiometry, direction, EC for each MNXR reaction

We download and cache these files locally; they are released under CC BY 4.0
and do NOT carry KEGG or HMDB downstream restrictions.
"""

from __future__ import annotations

import gzip
import logging
import os
from pathlib import Path
from typing import Iterator, Optional
from urllib.request import urlretrieve

import pandas as pd

log = logging.getLogger(__name__)

_FTP_BASE = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/"
_FILES = {
    "chem_xref": "chem_xref.tsv",
    "reac_xref": "reac_xref.tsv",
    "chem_prop": "chem_prop.tsv",
    "reac_prop": "reac_prop.tsv",
}


def _mnx_header_info(path: Path) -> tuple[list[str], int]:
    """
    Scan a MetaNetX TSV and return (column_names, data_start_line).

    MetaNetX files prefix their column-header row with '#' (e.g. '#ID\\t...')
    alongside hundreds of other '#'-comment lines.  This function finds the
    last '#'-prefixed line with ≥3 non-empty tab-separated fields whose first
    field is a bare word — that is the column header.
    """
    header_cols: list[str] = []
    data_start: int = 0
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh):
            if line.startswith("#"):
                parts = line[1:].rstrip("\n").split("\t")
                non_empty = [p for p in parts if p.strip()]
                if len(non_empty) >= 3 and parts[0] and " " not in parts[0]:
                    header_cols = parts
                    data_start = lineno + 1
            else:
                break
    return header_cols, data_start


def _read_mnx_tsv(path: Path) -> pd.DataFrame:
    """
    Read a MetaNetX flat-file TSV correctly.

    MetaNetX files have hundreds of '#'-prefixed comment lines followed by a
    column header that is *also* '#'-prefixed (e.g. '#ID\\tname\\t...' or
    '#source\\tID\\t...').  Using pd.read_csv(comment='#') drops that header,
    so the first data row becomes the column names — producing garbage.
    """
    header_cols, data_start = _mnx_header_info(path)
    if not header_cols:
        return pd.read_csv(path, sep="\t", comment="#", header=0, low_memory=False)
    return pd.read_csv(
        path, sep="\t", skiprows=data_start, header=None, names=header_cols, low_memory=False
    )


class MetaNetXClient:
    """
    Downloads and parses MetaNetX reference files for ID mapping and
    stoichiometry enrichment.
    """

    def __init__(self, cache_dir: str | Path = "data/raw/metanetx") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._chem_xref: Optional[pd.DataFrame] = None
        self._chem_prop: Optional[pd.DataFrame] = None
        self._reac_xref: Optional[pd.DataFrame] = None
        self._reac_prop: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # Download helpers
    # ------------------------------------------------------------------

    def download(self, force: bool = False) -> None:
        """Download all reference TSVs if not already cached."""
        for key, fname in _FILES.items():
            dest = self.cache_dir / fname
            if dest.exists() and not force:
                log.info("Cache hit: %s", dest)
                continue
            url = _FTP_BASE + fname
            log.info("Downloading %s → %s", url, dest)
            urlretrieve(url, dest)

    def _load(self, key: str) -> pd.DataFrame:
        path = self.cache_dir / _FILES[key]
        if not path.exists():
            raise FileNotFoundError(f"MetaNetX file not found: {path}. Run .download() first.")
        return _read_mnx_tsv(path)

    # ------------------------------------------------------------------
    # Compound cross-references
    # ------------------------------------------------------------------

    @property
    def chem_xref(self) -> pd.DataFrame:
        if self._chem_xref is None:
            df = self._load("chem_xref")
            # columns: source, ID, mnx_id, description
            self._chem_xref = df
        return self._chem_xref

    @property
    def chem_prop(self) -> pd.DataFrame:
        if self._chem_prop is None:
            self._chem_prop = self._load("chem_prop")
        return self._chem_prop

    def chebi_to_mnx(self, chebi_ids: list[str]) -> dict[str, str]:
        """Map ChEBI IDs to MetaNetX MNX compound IDs."""
        xref = self.chem_xref
        # source column looks like "chebi:CHEBI:15422" or "chebi:15422"
        chebi_col = xref[xref["source"].str.startswith("chebi:", na=False)].copy()
        chebi_col["chebi_id"] = chebi_col["source"].str.replace(
            r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
        )
        mapping = chebi_col.set_index("chebi_id")["ID"].to_dict()
        return {cid: mapping[cid] for cid in chebi_ids if cid in mapping}

    def mnx_compound_properties(self, mnx_ids: list[str]) -> pd.DataFrame:
        """Return chem_prop rows for a list of MNX compound IDs."""
        prop = self.chem_prop
        id_col = prop.columns[0]  # "ID" or "mnx_id" depending on version
        return prop[prop[id_col].isin(mnx_ids)].copy()

    # ------------------------------------------------------------------
    # Reaction cross-references
    # ------------------------------------------------------------------

    @property
    def reac_xref(self) -> pd.DataFrame:
        if self._reac_xref is None:
            self._reac_xref = self._load("reac_xref")
        return self._reac_xref

    @property
    def reac_prop(self) -> pd.DataFrame:
        if self._reac_prop is None:
            self._reac_prop = self._load("reac_prop")
        return self._reac_prop

    def reactome_to_mnxr(self, reactome_stids: list[str]) -> dict[str, str]:
        """Map Reactome stIDs to MetaNetX MNXR reaction IDs."""
        xref = self.reac_xref
        reactome_rows = xref[xref["source"].str.startswith("reactome:", na=False)].copy()
        reactome_rows["reactome_stid"] = reactome_rows["source"].str.replace(
            "reactome:", "", regex=False
        )
        mapping = reactome_rows.set_index("reactome_stid")["ID"].to_dict()
        return {sid: mapping[sid] for sid in reactome_stids if sid in mapping}
