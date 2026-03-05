"""
Metabolon data dictionary loader and ChEBI mapper.

Input: Metabolon global data dictionary CSV (PMC Open Access subset, 2024-04-14)
Columns: BIOCHEMICAL, PATHWAY, ChemRICHClass, PUBCHEM, INCHIKEY, PLATFORM, MASS, RI

ChEBI mapping strategy (in order of preference):
  1. InChIKey → MetaNetX chem_prop.tsv local lookup (fast, no rate limiting)
  2. InChIKey → ChEBI OLS4 REST API
  3. PubChem CID → PubChem PUG REST → InChIKey → step 1/2
  4. Unmatched: node_id = "metabolon:{BIOCHEMICAL}" (name-based fallback)

Match confidence levels:
  "exact_inchikey"   — InChIKey matched in MetaNetX / OLS4
  "pubchem_inchikey" — resolved via PubChem CID → InChIKey → ChEBI
  "unmatched"        — no ChEBI found

Usage::

    from gizmo.sources.metabolon import MetabolonLoader
    loader = MetabolonLoader("gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv")
    # Fast path: use MetaNetX local file for bulk InChIKey matching
    loader.load_metanetx_index("data/raw/metanetx/chem_prop.tsv", "data/raw/metanetx/chem_xref.tsv")
    nodes, report = loader.to_metabolite_nodes()
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from gizmo.schema import MetaboliteNode
from gizmo.sources.metanetx import _mnx_header_info

log = logging.getLogger(__name__)

_OLS4_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"
_PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"


@dataclass
class MatchReport:
    total: int = 0
    exact_inchikey: int = 0        # InChIKey exact match via MetaNetX (local)
    connectivity_inchikey: int = 0 # InChIKey connectivity-layer match via MetaNetX (local)
    pubchem_mnx: int = 0           # PubChem CID match via MetaNetX chem_xref (local)
    api_inchikey: int = 0          # InChIKey match via ChEBI OLS4 API
    pubchem_inchikey: int = 0      # PubChem CID → InChIKey → ChEBI via API
    unmatched: int = 0

    @property
    def chebi_coverage(self) -> float:
        if self.total == 0:
            return 0.0
        matched = (self.exact_inchikey + self.connectivity_inchikey
                   + self.pubchem_mnx + self.api_inchikey + self.pubchem_inchikey)
        return matched / self.total

    def __str__(self) -> str:
        return (
            f"Metabolon mapping: {self.total} compounds | "
            f"InChIKey(exact): {self.exact_inchikey} | "
            f"InChIKey(connectivity): {self.connectivity_inchikey} | "
            f"PubChem(local): {self.pubchem_mnx} | "
            f"InChIKey(API): {self.api_inchikey} | "
            f"PubChem(API): {self.pubchem_inchikey} | "
            f"Unmatched: {self.unmatched} | "
            f"ChEBI coverage: {self.chebi_coverage*100:.1f}%"
        )


class MetabolonLoader:
    """
    Load Metabolon data dictionary and map compounds to MetaboliteNode objects.

    The core matching is done via InChIKey against the MetaNetX chem_prop table
    (which is cached locally). This avoids any API rate-limiting for bulk loads.
    """

    def __init__(self, csv_path: str | Path) -> None:
        self.csv_path = Path(csv_path)
        self._df: Optional[pd.DataFrame] = None
        # InChIKey → ChEBI ID (from MetaNetX chem_prop + chem_xref, exact match)
        self._inchikey_to_chebi: dict[str, str] = {}
        # InChIKey connectivity prefix (14 chars) → ChEBI ID (fallback for stereoisomers)
        self._connectivity_to_chebi: dict[str, str] = {}
        # PubChem CID string → ChEBI ID (from MetaNetX chem_xref pubchem: entries)
        self._pubchem_to_chebi: dict[str, str] = {}
        # MNX compound ID → ChEBI ID (from chem_xref)
        self._mnx_to_chebi: dict[str, str] = {}

    # ------------------------------------------------------------------
    # Index loading
    # ------------------------------------------------------------------

    def load_metanetx_index(
        self,
        chem_prop_path: str | Path,
        chem_xref_path: str | Path,
    ) -> int:
        """
        Build InChIKey → ChEBI lookup from MetaNetX flat files.

        Streams both files line-by-line to avoid loading 700MB+ DataFrames.
        Returns number of indexed InChIKey → ChEBI entries.

        chem_prop.tsv columns (v4.x): ID, name, reference, formula, charge, mass, InChI, InChIKey, SMILES
        chem_xref.tsv columns: source, ID, description
          (source looks like "chebi:CHEBI:15422" or "chebi:15422")
        """
        import re

        chem_prop_path = Path(chem_prop_path)
        chem_xref_path = Path(chem_xref_path)

        if not chem_prop_path.exists() or not chem_xref_path.exists():
            log.warning(
                "MetaNetX files not found (%s, %s). Run MetaNetXClient().download() first.",
                chem_prop_path, chem_xref_path,
            )
            return 0

        # --- Pass 1: stream chem_xref ---
        #   chebi: rows  → {MNX_ID: ChEBI_ID}
        #   pubchem: rows → {pubchem_cid: MNX_ID}  (joined with chebi map after pass)
        xref_cols, xref_start = _mnx_header_info(chem_xref_path)
        # Typical columns: source (0), ID (1), description (2)
        mnx_idx = xref_cols.index("ID") if "ID" in xref_cols else 1

        mnx_to_chebi: dict[str, str] = {}
        mnx_from_pubchem: dict[str, str] = {}   # pubchem_cid → mnx_id (temporary)
        with open(chem_xref_path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh):
                if lineno < xref_start:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= mnx_idx:
                    continue
                source = parts[0]
                mnx_id = parts[mnx_idx]
                if source.startswith("chebi:"):
                    # "chebi:CHEBI:15422" → "CHEBI:15422";  "chebi:15422" → "CHEBI:15422"
                    chebi_id = "CHEBI:" + re.sub(r"^chebi:(?:CHEBI:)?", "", source)
                    mnx_to_chebi[mnx_id] = chebi_id
                elif source.startswith("pubchem:"):
                    cid = source[8:]   # strip "pubchem:"
                    mnx_from_pubchem[cid] = mnx_id

        if not mnx_to_chebi:
            log.warning("No ChEBI entries found in chem_xref.tsv")
            return 0

        self._mnx_to_chebi = mnx_to_chebi
        # Join pubchem → mnx with mnx → chebi to get pubchem → chebi
        self._pubchem_to_chebi = {
            cid: mnx_to_chebi[mid]
            for cid, mid in mnx_from_pubchem.items()
            if mid in mnx_to_chebi
        }
        log.info(
            "chem_xref: %d MNX→ChEBI, %d PubChem→ChEBI entries",
            len(mnx_to_chebi), len(self._pubchem_to_chebi),
        )

        # --- Pass 2: stream chem_prop → build {InChIKey: ChEBI_ID} ---
        prop_cols, prop_start = _mnx_header_info(chem_prop_path)

        try:
            id_idx = prop_cols.index("ID")
        except ValueError:
            id_idx = 0
        try:
            ik_idx = next(i for i, c in enumerate(prop_cols) if c.lower() == "inchikey")
        except StopIteration:
            log.warning("No InChIKey column found in chem_prop.tsv")
            return 0

        self._inchikey_to_chebi = {}
        self._connectivity_to_chebi = {}
        count = 0
        with open(chem_prop_path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh):
                if lineno < prop_start:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(id_idx, ik_idx):
                    continue
                mnx_id = parts[id_idx]
                inchikey = parts[ik_idx].strip()
                if inchikey and mnx_id in mnx_to_chebi:
                    chebi = mnx_to_chebi[mnx_id]
                    self._inchikey_to_chebi[inchikey] = chebi
                    # Connectivity prefix index (first 14 chars = connectivity layer).
                    # Only store the first entry seen — gives a canonical ChEBI for
                    # any stereoisomer / charge variant with the same skeleton.
                    pfx = inchikey[:14]
                    if pfx not in self._connectivity_to_chebi:
                        self._connectivity_to_chebi[pfx] = chebi
                    count += 1

        log.info(
            "MetaNetX index: %d exact InChIKey, %d connectivity prefix entries",
            count, len(self._connectivity_to_chebi),
        )
        return count

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self.csv_path, encoding="utf-8-sig")
        return self._df

    def raw_df(self) -> pd.DataFrame:
        """Return the raw Metabolon CSV as a DataFrame."""
        return self._load_df().copy()

    # ------------------------------------------------------------------
    # Node generation
    # ------------------------------------------------------------------

    def to_metabolite_nodes(
        self,
        *,
        api_fallback: bool = False,
        api_rate_limit_s: float = 0.2,
    ) -> tuple[list[MetaboliteNode], MatchReport]:
        """
        Convert Metabolon CSV rows to MetaboliteNode objects with ChEBI IDs where possible.

        Parameters
        ----------
        api_fallback : bool
            If True, use ChEBI OLS4 / PubChem REST for rows not matched via MetaNetX.
        api_rate_limit_s : float
            Seconds to sleep between API calls (only when api_fallback=True).

        Returns
        -------
        (nodes, report)
        """
        import time

        df = self._load_df()
        report = MatchReport(total=len(df))
        nodes: list[MetaboliteNode] = []

        for _, row in df.iterrows():
            biochemical = str(row.get("BIOCHEMICAL", "")).strip()
            inchikey = _clean(row.get("INCHIKEY"))
            pubchem_cid = _clean(row.get("PUBCHEM"))
            pathway = _clean(row.get("PATHWAY"))
            platform = _clean(row.get("PLATFORM"))
            mass_raw = row.get("MASS")
            mass = float(mass_raw) if _is_numeric(mass_raw) else None
            ri_raw = row.get("RI")
            ri = float(ri_raw) if _is_numeric(ri_raw) else None

            chebi_id: Optional[str] = None
            matched_inchikey: Optional[str] = inchikey
            confidence = "unmatched"

            # --- Match 1: MetaNetX InChIKey exact (local) ---
            if inchikey and inchikey in self._inchikey_to_chebi:
                chebi_id = self._inchikey_to_chebi[inchikey]
                confidence = "exact_inchikey"

            # --- Match 1b: MetaNetX InChIKey connectivity layer (local) ---
            # Catches stereoisomers / charge variants with the same carbon skeleton.
            elif inchikey and len(inchikey) >= 14 and inchikey[:14] in self._connectivity_to_chebi:
                chebi_id = self._connectivity_to_chebi[inchikey[:14]]
                confidence = "connectivity_inchikey"

            # --- Match 2: MetaNetX PubChem CID index (local) ---
            elif pubchem_cid and pubchem_cid in self._pubchem_to_chebi:
                chebi_id = self._pubchem_to_chebi[pubchem_cid]
                confidence = "pubchem_mnx"

            # --- Match 3: API fallback ---
            elif api_fallback:
                if inchikey:
                    chebi_id = _chebi_from_inchikey_api(inchikey)
                    if chebi_id:
                        confidence = "api_inchikey"

                if not chebi_id and pubchem_cid:
                    try:
                        resolved_ik = _inchikey_from_pubchem(pubchem_cid)
                        if resolved_ik:
                            matched_inchikey = resolved_ik
                            chebi_id = _chebi_from_inchikey_api(resolved_ik)
                            if chebi_id:
                                confidence = "pubchem_inchikey"
                    except Exception as exc:
                        log.debug("PubChem lookup failed for %s: %s", pubchem_cid, exc)

                if api_rate_limit_s > 0:
                    time.sleep(api_rate_limit_s)

            # --- Build node_id ---
            if chebi_id:
                node_id = chebi_id
            elif pubchem_cid:
                node_id = f"pubchem:{pubchem_cid}"
            else:
                # Sanitise biochemical name for use as ID
                safe_name = biochemical.replace(" ", "_").replace("/", "-")[:80]
                node_id = f"metabolon:{safe_name}"

            node = MetaboliteNode(
                node_id=node_id,
                chebi_id=chebi_id,
                pubchem_cid=pubchem_cid,
                metabolon_name=biochemical,
                name=biochemical,
                inchikey=matched_inchikey,
                mass=mass,
                platform=platform,
                retention_index=ri,
            )
            nodes.append(node)

            if confidence == "exact_inchikey":
                report.exact_inchikey += 1
            elif confidence == "connectivity_inchikey":
                report.connectivity_inchikey += 1
            elif confidence == "pubchem_mnx":
                report.pubchem_mnx += 1
            elif confidence == "api_inchikey":
                report.api_inchikey += 1
            elif confidence == "pubchem_inchikey":
                report.pubchem_inchikey += 1
            else:
                report.unmatched += 1

        log.info("%s", report)
        return nodes, report


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _chebi_from_inchikey_api(inchikey: str) -> Optional[str]:
    """Query ChEBI OLS4 REST for a ChEBI ID by InChIKey. Returns first match or None."""
    try:
        resp = requests.get(
            _OLS4_SEARCH,
            params={
                "q": inchikey,
                "ontology": "chebi",
                "fieldList": "id,obo_id",
                "exact": "true",
                "rows": 1,
            },
            timeout=15,
        )
        resp.raise_for_status()
        docs = resp.json().get("response", {}).get("docs", [])
        if docs:
            obo_id = docs[0].get("obo_id", "")
            if obo_id.startswith("CHEBI:"):
                return obo_id
    except Exception as exc:
        log.debug("ChEBI OLS4 lookup failed for %s: %s", inchikey, exc)
    return None


def _inchikey_from_pubchem(cid: str) -> Optional[str]:
    """Query PubChem PUG REST to get InChIKey for a CID."""
    try:
        url = _PUBCHEM_REST.format(cid=cid.strip())
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        props = resp.json().get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("InChIKey")
    except Exception as exc:
        log.debug("PubChem lookup failed for CID %s: %s", cid, exc)
    return None


def _clean(val: object) -> Optional[str]:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    # Convert float-encoded integers (e.g. 5793.0) to int strings
    if isinstance(val, float) and val == int(val):
        val = int(val)
    s = str(val).strip()
    return s if s else None


def _is_numeric(val: object) -> bool:
    try:
        float(val)  # type: ignore[arg-type]
        return True
    except (TypeError, ValueError):
        return False
