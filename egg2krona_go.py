#!/usr/bin/env python3
"""
Turn EggNOG-mapper annotations into a Krona plot of GO-term abundance.
"""

import argparse
import collections
import pathlib
import pandas as pd
from goatools.obo_parser import GODag
import subprocess
import sys
import tempfile


# ────────────────────────────────────────────────────────────────────
def find_header_row(xls_path: pathlib.Path) -> int:
    """Return the 0-based row index whose first cell is exactly 'query'."""
    col0 = pd.read_excel(xls_path, usecols=[0], header=None)[0]
    try:
        return col0[col0 == "query"].index[0]
    except IndexError:
        sys.exit(
            f"[ERROR] Could not find a row whose first cell is 'query' in {xls_path}"
        )


def load_gos(xls_path: pathlib.Path) -> pd.Series:
    """Read the Excel file and return the 'GOs' column as a Series."""
    header_row = find_header_row(xls_path)
    df = pd.read_excel(xls_path, header=header_row)
    if "GOs" not in df.columns:
        sys.exit("[ERROR] A column named 'GOs' was not found.")
    return df["GOs"].fillna("-")  # replace NaNs with '-'


def count_go_terms(go_series: pd.Series) -> collections.Counter:
    """Count GO IDs across all rows."""
    counts = collections.Counter()
    for gos in go_series:
        if gos == "-" or not gos.strip():
            continue
        for go_id in map(str.strip, gos.split(",")):
            if go_id.startswith("GO:"):
                counts[go_id] += 1
    return counts


def lineage_from_root(term, namespace_roots):
    """
    Build a lineage list (root ➞ … ➞ term).
    If a term has multiple parents we pick the first (deterministically).
    """
    lineage = [term.name]
    node = term
    while node.id not in namespace_roots:
        parents = sorted(node.parents, key=lambda t: t.id)
        if not parents:
            break
        node = parents[0]
        lineage.append(node.name)
    namespace = term.namespace.replace("_", " ").title()
    lineage[-1] = namespace  # replace root name
    return list(reversed(lineage))


def write_krona_input(counter, godag, out_tsv):
    """Write Krona TSV with COUNT first."""
    namespace_roots = {
        "GO:0008150",  # biological_process
        "GO:0003674",  # molecular_function
        "GO:0005575",  # cellular_component
    }
    with open(out_tsv, "w") as fh:
        for go_id, n in counter.items():
            if go_id not in godag:
                continue
            lineage = lineage_from_root(godag[go_id], namespace_roots)
            # *** FIX: count first! ***
            fh.write(f"{n}\t" + "\t".join(lineage) + "\n")


def run_krona(tsv_path, html_out):
    """Call ktImportText to generate the interactive plot."""
    try:
        subprocess.run(
            ["ktImportText", str(tsv_path), "-o", str(html_out)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        sys.exit(
            "KronaTools not found. Install it (`conda install -c bioconda krona`) "
            "or run ktImportText manually later."
        )


# ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Create a Krona plot of GO-term abundance from EggNOG output."
    )
    ap.add_argument("eggnog_excel", type=pathlib.Path)
    ap.add_argument(
        "--obo",
        required=True,
        type=pathlib.Path,
        help="Path to go-basic.obo (download from geneontology.org)",
    )
    ap.add_argument(
        "--out",
        default="krona_go.html",
        type=pathlib.Path,
        help="Name of output HTML (default: krona_go.html)",
    )
    args = ap.parse_args()

    print("[1/4] Reading EggNOG file…")
    go_series = load_gos(args.eggnog_excel)

    print("[2/4] Counting GO terms…")
    counts = count_go_terms(go_series)
    print(f"    → {len(counts):,} distinct GO terms, {sum(counts.values()):,} hits")

    print("[3/4] Loading GO DAG…")
    godag = GODag(str(args.obo), optional_attrs={"namespace"})

    print("[4/4] Building Krona plot…")
    with tempfile.TemporaryDirectory() as tmpdir:
        krona_tsv = pathlib.Path(tmpdir) / "krona_go.tsv"
        write_krona_input(counts, godag, krona_tsv)
        run_krona(krona_tsv, args.out)

    print(f"✓ Done! Open {args.out} in your browser.")


if __name__ == "__main__":
    main()
