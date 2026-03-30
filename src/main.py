#!/usr/bin/env python3
"""
main.py

Purpose:
    Main driver for the ORF analysis pipeline.

Role in Project:
    Calls other modules to handle sequence loading, ORF detection,
    and CSV reporting.

Input:
    DNA sequence fetched from NCBI via accession number.

Output:
    - output/orfs.csv          : flat table of every ORF found (one row per ORF)
    - output/orf_stats.csv     : per-ORF statistics (GC content, protein length, etc.)
    - output/stats_summary.txt : human-readable summary report

Example usage
-------------
# All defaults (ATG only, min length 30, nested ORFs included)
python main.py --accession NM_001301717 --email you@example.com

# All three start codons, minimum 60 nt, ignore nested ORFs
python main.py --accession NM_001301717 --email you@example.com \
    --start-codons ATG GTG TTG --min-length 60 --ignore-nested
"""

import argparse
import csv
import os
import sys

from src.input_lib.input_validate import run as validate_run
from src.orf_finder_lib.orf_finder import find_orfs, CSV_FIELDNAMES
from src.statistics_lib.statistics_summary import (
    calculate_orf_stats,
    write_stats_to_file,
    ORF_STATS_FIELDNAMES,
)

# Valid start codons the user is allowed to request
VALID_START_CODONS = {"ATG", "GTG", "TTG"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_nested(flat_list: list) -> list:
    """
    Return the subset of ORFs that are nested inside another ORF.
    An ORF is nested if its start position falls within another ORF's
    start-end range on the same strand and reading frame.
    """
    nested_orfs = []
    for i, orf in enumerate(flat_list):
        for j, other in enumerate(flat_list):
            if i == j:
                continue
            if orf["strand"] != other["strand"]:
                continue
            if orf["frame"] != other["frame"]:
                continue
            if other["end"] is None:
                continue
            if other["start"] < orf["start"] < other["end"]:
                nested_orfs.append(orf)
                break
    return nested_orfs


def _print_summary(nested: dict, flat_list: list) -> None:
    """Print a short summary of ORF counts to stdout."""
    complete   = nested["complete"]
    incomplete = nested["incomplete"]

    n_complete_canonical      = len(complete["canonical"])
    n_incomplete_canonical    = len(incomplete["canonical"])

    n_complete_noncanonical   = sum(
        len(v) for v in complete["noncanonical"].values()
    )
    n_incomplete_noncanonical = sum(
        len(v) for v in incomplete["noncanonical"].values()
    )

    total = (
        n_complete_canonical
        + n_incomplete_canonical
        + n_complete_noncanonical
        + n_incomplete_noncanonical
    )

    plus_strand  = sum(1 for o in flat_list if o.get("strand") == "+")
    minus_strand = sum(1 for o in flat_list if o.get("strand") == "-")

    print("\n========== ORF Summary ==========")
    print(f"  Total ORFs found            : {total}")
    print(f"  Forward strand (+)          : {plus_strand}")
    print(f"  Reverse strand (-)          : {minus_strand}")
    print(f"  Complete   (ATG)            : {n_complete_canonical}")
    print(f"  Incomplete (ATG)            : {n_incomplete_canonical}")
    print(f"  Complete   (non-canonical)  : {n_complete_noncanonical}")
    print(f"  Incomplete (non-canonical)  : {n_incomplete_noncanonical}")

    for sc in ("GTG", "TTG"):
        nc = len(complete["noncanonical"].get(sc, {}))
        ni = len(incomplete["noncanonical"].get(sc, {}))
        if nc + ni > 0:
            print(f"    {sc} — complete: {nc}, incomplete: {ni}")

    nested_found = _find_nested(flat_list)
    print(f"  Nested ORFs detected        : {len(nested_found)}")
    print("=================================\n")


def _write_csv(flat_list: list, output_path: str) -> None:
    """Write the flat ORF list to a CSV file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat_list)
    print(f"[INFO] ORF table written to: {output_path}")


def _write_csv_with_fields(data: list, output_path: str, fieldnames: list) -> None:
    """Write a list of dicts to a CSV file with specified fieldnames."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True) if os.path.dirname(output_path) else None
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)
    print(f"[INFO] Table written to: {output_path}")


def _validate_start_codons(requested: list) -> list:
    """
    Upper-case and validate the user-supplied start codons.
    Exits with a helpful message if any unrecognised codon is given.
    """
    upper   = [c.upper() for c in requested]
    unknown = [c for c in upper if c not in VALID_START_CODONS]
    if unknown:
        print(
            f"[ERROR] Unrecognised start codon(s): {', '.join(unknown)}\n"
            f"        Allowed values are: {', '.join(sorted(VALID_START_CODONS))}"
        )
        sys.exit(1)
    return upper


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ORF analysis pipeline — fetches a sequence from NCBI "
                    "and reports all ORFs as a CSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Sequence retrieval ---
    parser.add_argument(
        "--accession",
        type=str,
        help="NCBI accession number (e.g. NM_001301717)",
    )
    parser.add_argument(
        "--email",
        type=str,
        help="Email address required by NCBI Entrez",
    )

    # --- ORF finder options ---
    parser.add_argument(
        "--min-length",
        type=int,
        default=30,
        metavar="NT",
        help="Minimum ORF length in nucleotides",
    )
    parser.add_argument(
        "--start-codons",
        nargs="+",
        default=["ATG"],
        metavar="CODON",
        help=(
            "One or more start codons to search for. "
            "Supply as a space-separated list, e.g.: --start-codons ATG GTG TTG"
        ),
    )
    parser.add_argument(
        "--ignore-nested",
        action="store_true",
        default=False,
        help="Exclude ORFs whose start falls inside another ORF in the same frame",
    )

    # --- Output options ---
    parser.add_argument(
        "--output",
        type=str,
        default="output/orfs.csv",
        help="Path for the output CSV file",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Get accession and email (prompt if not supplied as flags)
    # ------------------------------------------------------------------
    accession = args.accession or input("Enter NCBI accession number: ").strip()
    email     = args.email     or input("Enter your email (required by NCBI): ").strip()

    # ------------------------------------------------------------------
    # 2. Validate user-supplied options
    # ------------------------------------------------------------------
    start_codons = _validate_start_codons(args.start_codons)

    if args.min_length < 3:
        print("[ERROR] --min-length must be at least 3 (one codon).")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 3. Fetch and validate the sequence
    # ------------------------------------------------------------------
    acc, clean_seq = validate_run(accession, email)

    if clean_seq is None:
        print("[ERROR] Pipeline failed: could not retrieve a valid sequence.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 4. Run the ORF finder
    # ------------------------------------------------------------------
    nested, flat_list = find_orfs(
        clean_seq,
        start_codons=start_codons,
        min_length=args.min_length,
        ignore_nested=args.ignore_nested,
    )

    # ------------------------------------------------------------------
    # 5. Print terminal summary
    # ------------------------------------------------------------------
    _print_summary(nested, flat_list)

    if not flat_list:
        print("[WARNING] No ORFs found. No output files written.")
        return

    # ------------------------------------------------------------------
    # 6. Write output files
    # ------------------------------------------------------------------
    # Main ORF table
    _write_csv(flat_list, args.output)

    # Per-ORF stats CSV
    per_orf_stats = calculate_orf_stats(flat_list, clean_seq)
    _write_csv_with_fields(
        per_orf_stats,
        "output/orf_stats.csv",
        ORF_STATS_FIELDNAMES,
    )

    # Human-readable summary report
    write_stats_to_file(flat_list, clean_seq, acc, outfile="output/stats_summary.txt")


if __name__ == "__main__":
    main()

