# GIZMO

**Graph-Integrated Zone of Metabolite Operations**

A metabolite-centered reaction graph for systems biology, cheminformatics, and clinical translation.
Built entirely on open-licensed data sources — no HMDB or KEGG dependency.

---

## What it does

GIZMO builds a directed graph where **metabolites** and **reactions** are first-class nodes, enriched with **disease** and **gene** nodes from clinical databases. The graph is designed to plug directly into:

- **GATOM / shiny-gam** — GraphML export, bipartite topology
- **Graph kernels** — Weisfeiler-Lehman, random walk (currency-filtered subgraph)
- **PyTorch Geometric** — node-link JSON → `Data` object
- **COBRApy / FBA** — stoichiometry and reversibility preserved on every reaction node
- **Cytoscape / Gephi** — GraphML with typed nodes and scored edges

---

## Data sources

All sources are **CC BY 4.0** — safe for academic and commercial use downstream.

| Source | What it provides | License |
|---|---|---|
| [Reactome](https://reactome.org) | Reactions, pathways, EC numbers, gene associations | CC BY 4.0 |
| [ChEBI](https://www.ebi.ac.uk/chebi/) | Structural metadata (formula, InChI, SMILES, charge) | CC BY 4.0 |
| [MetaNetX](https://www.metanetx.org) | Cross-reference hub: ChEBI ↔ Reactome ↔ MNX IDs, stoichiometry | CC BY 4.0 |
| [MONDO](https://mondo.monarchinitiative.org) | Unified disease ontology, IEM subset, ICD-10/OMIM xrefs | CC BY 4.0 |
| [Open Targets](https://platform.opentargets.org) | Scored gene–disease associations (GWAS, rare variant, etc.) | CC BY 4.0 |
| [Orphanet / Orphadata](https://www.orphadata.com) | Rare disease classifications, gene–disease associations | CC BY 4.0 |
| Metabolon global dictionary | 5 395 measured metabolites with InChIKey, PubChem CID, platform | Internal / PMC OA |

---

## Installation

**Requirements:** Python ≥ 3.10, [Poetry](https://python-poetry.org/docs/#installation) ≥ 1.8

```bash
git clone <repo>
cd GIZMO

# Core install (graph, sources, analysis, export)
poetry install

# With Jupyter notebooks
poetry install --with notebooks

# With RDKit for structural fingerprints / cheminformatics
poetry install --with cheminformatics

# With PyTorch Geometric for GNN training
poetry install --with torch

# Everything
poetry install --with notebooks,cheminformatics,torch
```

Activate the environment:

```bash
poetry shell
# or prefix commands with: poetry run
```

---

## Quick start

```python
from gizmo.sources.reactome import ReactomeClient
from gizmo.graph.network import MetaboliteGraph
from gizmo.analysis.currency import flag_currency_metabolites
from gizmo.analysis.qc import assess_readiness

# 1. Build graph from Reactome
client = ReactomeClient()
mg = MetaboliteGraph()

reactions = client.pathway_reactions("R-HSA-70171")   # glycolysis
for stub in reactions:
    detail = client.reaction_detail(stub["stId"])
    rxn, mets, edges = client.parse_reaction(detail, pathway_stids=["R-HSA-70171"])
    mg.add_reaction(rxn)
    mg.add_metabolites(mets)
    mg.add_edges(edges)

# 2. Flag currency metabolites (ATP, NAD+, water, etc.)
flag_currency_metabolites(mg)

# 3. Run QC
report = assess_readiness(mg)
report.print_summary()

print(mg)
# MetaboliteGraph(42 metabolites [6 currency], 13 reactions, 0 diseases, 0 genes, 67 edges)
```

---

## One-time data downloads

These are large files (~750 MB total); run once and they are cached locally.

```python
# MetaNetX cross-reference tables (~500 MB)
from gizmo.sources.metanetx import MetaNetXClient
MetaNetXClient().download()

# MONDO disease ontology (~200 MB)
from gizmo.sources.mondo import MondoClient
MondoClient().download()

# Orphanet rare disease XML (~50 MB)
from gizmo.sources.orphanet import OrphanetClient
OrphanetClient().download()
```

---

## Project layout

```
GIZMO/
├── gizmo/
│   ├── schema.py              # Node and edge schemas (Pydantic v2)
│   ├── graph/
│   │   └── network.py         # MetaboliteGraph — core NetworkX DiGraph
│   ├── sources/
│   │   ├── reactome.py        # Reactome REST API client
│   │   ├── chebi.py           # ChEBI via EBI OLS4 REST
│   │   ├── metanetx.py        # MetaNetX flat-file download + ID mapping
│   │   ├── mondo.py           # MONDO OBO parser → DiseaseNode
│   │   ├── open_targets.py    # Open Targets GraphQL client
│   │   ├── orphanet.py        # Orphadata XML parser
│   │   └── metabolon.py       # Metabolon data dictionary loader + ChEBI mapper
│   ├── analysis/
│   │   ├── currency.py        # Currency metabolite detection and filtering
│   │   └── qc.py              # Computational readiness report
│   ├── export/
│   │   ├── graphml.py         # GraphML read/write (GATOM / Cytoscape)
│   │   └── json_export.py     # JSON node-link read/write (D3 / PyG)
│   └── cli.py                 # `gizmo` command-line interface
├── tests/                     # pytest suite (34 tests)
├── notebooks/
│   └── 01_gizmo_walkthrough.ipynb
├── data/
│   ├── raw/                   # Downloaded source files (gitignored)
│   └── processed/             # Exported graphs
├── pyproject.toml
└── poetry.lock
```

---

## Graph schema

### Node types

| `node_type` | Canonical `node_id` | Key attributes |
|---|---|---|
| `metabolite` | `CHEBI:XXXXX` or `CHEBI:XXXXX@compartment` | `chebi_id`, `metanetx_id`, `pubchem_cid`, `metabolon_name`, `formula`, `inchikey`, `mass`, `is_currency` |
| `reaction` | `reactome:R-HSA-XXXXX` | `ec_numbers`, `gene_symbols`, `pathways`, `reversible` |
| `disease` | `MONDO:XXXXXXX` | `xref_omim`, `xref_orphanet`, `xref_icd10`, `is_inborn_error_of_metabolism` |
| `gene` | `ENSG:ENSGXXXXXXXXXXXXX` | `ensembl_id`, `hgnc_id`, `symbol` |

### Edge types

| Edge class | Direction | `role` / `edge_type` |
|---|---|---|
| `ReactionEdge` | metabolite → reaction | `substrate`, `modifier` |
| `ReactionEdge` | reaction → metabolite | `product` |
| `DiseaseEdge` | disease → gene | `gene_associated` (Open Targets `score`) |
| `DiseaseEdge` | disease → metabolite | `biomarker`, `causal` |
| `DiseaseEdge` | disease → reaction | `pathway_associated` |
| `DiseaseEdge` | gene → reaction | `gene_reaction` |

---

## Currency metabolites

Currency metabolites (ATP, NAD⁺, water, etc.) appear in hundreds of reactions and create
spurious edges in metabolite co-occurrence graphs. GIZMO flags them rather than deleting
them — stoichiometry data is preserved for FBA.

```python
from gizmo.analysis.currency import flag_currency_metabolites, noncurrency_subgraph

# Flag in-place (canonical list + degree-outlier heuristic)
result = flag_currency_metabolites(mg)
print(result["canonical"])    # flagged by known ChEBI IDs
print(result["statistical"])  # flagged as degree outliers

# For graph kernels / co-occurrence: remove currency nodes
mg_clean = noncurrency_subgraph(mg)
```

Canonical list covers 24 metabolites including ATP/ADP/AMP, NAD⁺/NADH, NADP⁺/NADPH,
FAD/FADH₂, CoA, water, H⁺, Pi, PPi, O₂, CO₂. Borderline cases (L-glutamate,
2-oxoglutarate, glucose-6-phosphate) are excluded by default; pass `include_borderline=True`
to flag them.

---

## Metabolon compound mapping

```python
from gizmo.sources.metabolon import MetabolonLoader

loader = MetabolonLoader("gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv")

# Load MetaNetX InChIKey index for fast bulk matching (requires MetaNetX download)
loader.load_metanetx_index(
    "data/raw/metanetx/chem_prop.tsv",
    "data/raw/metanetx/chem_xref.tsv",
)

nodes, report = loader.to_metabolite_nodes(
    api_fallback=True,   # also hit ChEBI OLS4 / PubChem for unmatched rows
)
print(report)
# Metabolon mapping: 5395 compounds | InChIKey hit: 3241 | PubChem hit: 412 | Unmatched: 1742 | ChEBI coverage: 67.7%
```

Match confidence order:
1. MetaNetX InChIKey local index — instant, no network
2. ChEBI OLS4 REST — for InChIKeys absent from MetaNetX
3. PubChem CID → InChIKey → ChEBI — last resort
4. Fallback node ID: `metabolon:{BIOCHEMICAL_NAME}`

---

## Computational readiness QC

```python
from gizmo.analysis.qc import assess_readiness

report = assess_readiness(mg)
report.print_summary()   # rich table in terminal / notebook

report.is_fba_ready      # heuristic: >100 metabolites, >50 reactions,
                         # >30% EC coverage, >50% formula coverage, <10 components
```

The report checks:

- Currency edge fraction (should be < 50%)
- Dead-end metabolites (no upstream or downstream reaction — structural gaps)
- Orphan reactions (no substrates or no products found in graph)
- Weakly connected components and isolated nodes
- Annotation completeness: EC numbers, gene symbols, pathway membership
- Structural data coverage: molecular formula, InChIKey, ChEBI ID
- Compartment inventory
- Clinical overlay: disease–metabolite and disease–gene edge counts
- Metabolon → ChEBI resolution coverage

---

## CLI

```bash
# Computational readiness report on a saved graph
gizmo qc --graph data/processed/gizmo_graph.json

# Metabolon ChEBI mapping report
gizmo metabolon \
  --csv gizmo/sources/metabolon_data_dictionary_PMC_OA_subset_4.14.2024.csv \
  --metanetx-prop data/raw/metanetx/chem_prop.tsv \
  --metanetx-xref data/raw/metanetx/chem_xref.tsv
```

---

## Export formats

| Format | Use case | Function |
|---|---|---|
| JSON node-link | Lossless round-trip, D3.js, PyTorch Geometric | `write_json` / `read_json` |
| GraphML | GATOM, igraph (R), Cytoscape, Gephi | `write_graphml` / `read_graphml` |

List-valued node attributes (`ec_numbers`, `pathways`, `synonyms`) are
serialised as pipe-delimited strings in GraphML and restored on read.

```python
from gizmo.export.json_export import write_json, read_json
from gizmo.export.graphml import write_graphml, read_graphml

write_json(mg, "graph.json")
write_graphml(mg, "graph.graphml")

# Load in R
# library(igraph)
# g <- read_graph("graph.graphml", format = "graphml")
```

---

## Tests

```bash
poetry run pytest          # 34 tests
poetry run pytest -v       # verbose
poetry run pytest --cov=gizmo --cov-report=term-missing
```

---

## Disease and clinical integration

```python
from gizmo.sources.mondo import MondoClient
from gizmo.sources.orphanet import OrphanetClient
from gizmo.sources.open_targets import OpenTargetsClient

# Inborn errors of metabolism from MONDO
mondo = MondoClient()
iem = mondo.load_iem_subset()
mg.add_diseases(iem)

# Gene–disease associations from Orphanet
orphanet = OrphanetClient()
genes, edges = orphanet.load_gene_associations()
mg.add_genes(genes)
mg.add_disease_edges(edges)

# Scored associations from Open Targets (live API)
ot = OpenTargetsClient()
genes, edges = ot.gene_associations_for_disease("MONDO:0009861", min_score=0.1)
mg.add_genes(genes)
mg.add_disease_edges(edges)
```

---

## Notebook

```bash
poetry install --with notebooks
poetry run jupyter notebook notebooks/01_gizmo_walkthrough.ipynb
```

The walkthrough notebook covers the full pipeline end-to-end with visualisations.

---

## License

MIT. All embedded data is CC BY 4.0 (see Data sources above).
OMIM IDs appear as cross-reference strings only; no OMIM content is redistributed.
