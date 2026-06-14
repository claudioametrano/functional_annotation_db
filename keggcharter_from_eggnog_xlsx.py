#!/usr/bin/env python3
"""
eggNOG-mapper XLSX -> KEGGCharter, simplified and avoiding KEGG organism-list download.

This version calls KEGGCharter with -it TAXON instead of -tc Taxon.
That avoids KEGGCharter's problematic attempt to download:
    http://rest.kegg.jp/list/organism

Use this for one custom Prokka/eggNOG genome where you only want to color
KEGG reference maps with your detected KOs.
"""

import argparse
import re
import subprocess
from pathlib import Path
import pandas as pd


def read_eggnog_xlsx(xlsx_file):
    """Read eggNOG Excel output and locate the real header row containing KEGG_ko."""
    raw = pd.read_excel(xlsx_file, header=None, dtype=str)

    header_row_index = None
    for i, row in raw.iterrows():
        row_values = [str(x).strip() for x in row.tolist()]
        if "KEGG_ko" in row_values:
            header_row_index = i
            break

    if header_row_index is None:
        raise ValueError("Could not find a column named 'KEGG_ko'.")

    header = raw.iloc[header_row_index].tolist()
    data = raw.iloc[header_row_index + 1:].copy()
    data.columns = header
    data = data.dropna(how="all")
    return data


def extract_unique_kos(data):
    """
    Extract unique KEGG KO identifiers from KEGG_ko.

    Handles values like:
        ko:K00844
        ko:K00844,ko:K12407
        K00844,K12407
        -
    """
    if "KEGG_ko" not in data.columns:
        raise ValueError("The table does not contain a 'KEGG_ko' column.")

    all_kos = set()
    for value in data["KEGG_ko"].dropna():
        for ko in re.findall(r"K\d{5}", str(value)):
            all_kos.add(ko)

    return sorted(all_kos)


def write_keggcharter_input(kos, output_file):
    """
    Write minimal KEGGCharter input.

    Only one column is needed:
        KO
        K00001
        K00002

    The sample/genome label is passed separately with KEGGCharter -it TAXON.
    """
    output_file = Path(output_file)
    with output_file.open("w") as out:
        out.write("KO\n")
        for ko in kos:
            out.write(f"{ko}\n")


def clean_map_argument(map_argument):
    """
    Accept:
        --map 00195
        --map ko00195
        --map map00195
        --map 00195,00710,00910
    Return:
        ['00195'] or ['00195', '00710', '00910']
    """
    cleaned_maps = []

    for item in map_argument.split(","):
        item = item.strip()
        item = item.removeprefix("ko")
        item = item.removeprefix("map")

        if not re.fullmatch(r"\d{5}", item):
            raise ValueError(
                f"Invalid KEGG map ID: {item!r}. "
                "Use a five-digit KEGG pathway map ID, for example 00195."
            )

        if item not in cleaned_maps:
            cleaned_maps.append(item)

    return cleaned_maps


def run_keggcharter(input_tsv, output_dir, resource_dir, maps, taxon, threads=None, step=None):
    """
    Run KEGGCharter.

    Key point:
        use -it TAXON instead of -tc Taxon.

    This avoids the KEGG organism-list lookup that caused your crash.
    """
    command = [
        "keggcharter",
        "-f", str(input_tsv),
        "-koc", "KO",
        "-iq",
        "-it", str(taxon),
        "--map-all",
        "-mm", ",".join(maps),
        "-rd", str(resource_dir),
        "-o", str(output_dir),
    ]

    if threads is not None:
        command.extend(["-t", str(threads)])

    if step is not None:
        command.extend(["--step", str(step)])

    print("\nRunning KEGGCharter command:")
    print(" ".join(command))
    print()

    subprocess.run(command, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Convert eggNOG-mapper XLSX output to KEGGCharter maps."
    )

    parser.add_argument(
        "xlsx",
        help="eggNOG-mapper Excel output file, e.g. nostoc.emapper.annotations.xlsx",
    )

    parser.add_argument(
        "--map",
        required=True,
        help="KEGG pathway map ID. Example: --map 00195. Multiple: --map 00195,00710",
    )

    parser.add_argument(
        "--taxon",
        default="Genome",
        help="Custom genome/sample name shown by KEGGCharter. Default: Genome",
    )

    parser.add_argument(
        "--outdir",
        default="keggcharter_from_eggnog",
        help="Output directory. Default: keggcharter_from_eggnog",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Optional KEGGCharter threads, passed to -t.",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=None,
        help="Optional IDs per KEGG API request, passed to --step. Try 10 if API is unstable.",
    )

    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually run KEGGCharter after creating the input table.",
    )

    args = parser.parse_args()

    xlsx_file = Path(args.xlsx)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    maps = clean_map_argument(args.map)

    keggcharter_input = outdir / "keggcharter_unique_KO.tsv"
    resource_dir = outdir / "keggcharter_resources"
    maps_outdir = outdir / "keggcharter_maps"

    print(f"Reading eggNOG Excel file: {xlsx_file}")
    data = read_eggnog_xlsx(xlsx_file)

    print("Extracting unique KEGG KOs from column: KEGG_ko")
    kos = extract_unique_kos(data)

    if len(kos) == 0:
        raise ValueError("No KEGG KO identifiers were found in KEGG_ko.")

    print(f"Number of unique KOs found: {len(kos)}")

    print(f"Writing KEGGCharter input table: {keggcharter_input}")
    write_keggcharter_input(kos, keggcharter_input)

    print("\nCreated KEGGCharter input file:")
    print(f"  {keggcharter_input}")

    print("\nRequested KEGG pathway map(s):")
    for m in maps:
        print(f"  {m}")

    print("\nTaxon/sample label passed to KEGGCharter with -it:")
    print(f"  {args.taxon}")

    if args.run:
        run_keggcharter(
            input_tsv=keggcharter_input,
            output_dir=maps_outdir,
            resource_dir=resource_dir,
            maps=maps,
            taxon=args.taxon,
            threads=args.threads,
            step=args.step,
        )
        print("\nKEGGCharter finished.")
        print("Check PNG maps in:")
        print(f"  {maps_outdir}/maps")
    else:
        print("\nKEGGCharter was NOT run because --run was not used.")
        print("To run KEGGCharter, use the same command with --run.")


if __name__ == "__main__":
    main()
