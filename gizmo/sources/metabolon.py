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

log = logging.getLogger(__name__)

_OLS4_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"
_PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"


@dataclass
class MatchReport:
    total: int = 0
    exact_inchikey: int = 0
    pubchem_inchikey: int = 0
    unmatched: int = 0

    @property
    def chebi_coverage(self) -> float:
        if self.total == 0:
            return 0.0
        return (self.exact_inchikey + self.pubchem_inchikey) / self.total

    def __str__(self) -> str:
        return (
            f"Metabolon mapping: {self.total} compounds | "
            f"InChIKey hit: {self.exact_inchikey} | "
            f"PubChem hit: {self.pubchem_inchikey} | "
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
        # InChIKey → ChEBI ID index (populated by load_metanetx_index)
        self._inchikey_to_chebi: dict[str, str] = {}
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
        Returns number of indexed entries.

        chem_prop.tsv columns (v4.x): ID, Description, Formula, Charge, Mass, InChI, InChIKey, SMILES
        chem_xref.tsv columns: source, ID, description
          (source looks like "chebi:CHEBI:15422")
        """
        chem_prop_path = Path(chem_prop_path)
        chem_xref_path = Path(chem_xref_path)

        if not chem_prop_path.exists() or not chem_xref_path.exists():
            log.warning(
                "MetaNetX files not found (%s, %s). Run MetaNetXClient().download() first.",
                chem_prop_path, chem_xref_path,
            )
            return 0

        # Load chem_xref to build MNX → ChEBI
        xref_df = pd.read_csv(chem_xref_path, sep="\t", comment="#", header=0, low_memory=False)
        # Find the ChEBI source column
        chebi_rows = xref_df[xref_df.iloc[:, 0].str.startswith("chebi:", na=False)].copy()
        if chebi_rows.empty:
            log.warning("No ChEBI entries found in chem_xref.tsv")
        else:
            # columns[0] = source (chebi:CHEBI:XXXXX), columns[1] = MNX ID
            chebi_rows = chebi_rows.copy()
            chebi_rows["chebi_id"] = chebi_rows.iloc[:, 0].str.replace(
                r"^chebi:(?:CHEBI:)?", "CHEBI:", regex=True
            )
            mnx_col = chebi_rows.columns[1]
            self._mnx_to_chebi = dict(zip(chebi_rows[mnx_col], chebi_rows["chebi_id"]))

        # Load chem_prop to build InChIKey → MNX → ChEBI
        prop_df = pd.read_csv(chem_prop_path, sep="\t", comment="#", header=0, low_memory=False)
        # Find InChIKey column (varies by version: "InChIKey", "inchikey")
        inchikey_col = next(
            (c for c in prop_df.columns if c.lower() == "inchikey"), None
        )
        id_col = prop_df.columns[0]  # MNX ID column

        if inchikey_col is None:
            log.warning("No InChIKey column found in chem_prop.tsv")
            return 0

        count = 0
        for _, row in prop_df.iterrows():
            ik = row.get(inchikey_col)
            mnx_id = row.get(id_col)
            if pd.isna(ik) or pd.isna(mnx_id):
                continue
            chebi_id = self._mnx_to_chebi.get(mnx_id)
            if chebi_id:
                self._inchikey_to_chebi[str(ik)] = chebi_id
                count += 1

        log.info("MetaNetX InChIKey index: %d entries mapped to ChEBI", count)
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

            # --- Match 1: MetaNetX InChIKey index ---
            if inchikey and inchikey in self._inchikey_to_chebi:
                chebi_id = self._inchikey_to_chebi[inchikey]
                confidence = "exact_inchikey"

            # --- Match 2: API fallback ---
            elif api_fallback:
                if inchikey:
                    chebi_id = _chebi_from_inchikey_api(inchikey)
                    if chebi_id:
                        confidence = "exact_inchikey"

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
                if "pubchem" in confidence:
                    report.pubchem_inchikey += 1
                else:
                    report.exact_inchikey += 1
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
