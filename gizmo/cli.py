"""
GIZMO command-line interface.

Usage:
  gizmo --help
  gizmo qc --graph graph.json
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gizmo",
        description="Graph-Integrated Zone of Metabolite Operations",
    )
    sub = parser.add_subparsers(dest="command")

    # qc subcommand
    qc_p = sub.add_parser("qc", help="Run computational readiness report on a saved graph")
    qc_p.add_argument("--graph", required=True, help="Path to graph JSON file")

    # metabolon subcommand
    met_p = sub.add_parser("metabolon", help="Load Metabolon data dictionary and report ChEBI coverage")
    met_p.add_argument("--csv", required=True, help="Path to Metabolon CSV")
    met_p.add_argument(
        "--metanetx-prop", default="data/raw/metanetx/chem_prop.tsv",
        help="MetaNetX chem_prop.tsv path"
    )
    met_p.add_argument(
        "--metanetx-xref", default="data/raw/metanetx/chem_xref.tsv",
        help="MetaNetX chem_xref.tsv path"
    )

    args = parser.parse_args(argv)

    if args.command == "qc":
        from gizmo.export.json_export import read_json
        from gizmo.analysis.qc import assess_readiness

        mg = read_json(args.graph)
        report = assess_readiness(mg)
        report.print_summary()

    elif args.command == "metabolon":
        from gizmo.sources.metabolon import MetabolonLoader

        loader = MetabolonLoader(args.csv)
        loader.load_metanetx_index(args.metanetx_prop, args.metanetx_xref)
        _, report = loader.to_metabolite_nodes()
        print(report)

    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
