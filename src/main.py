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
    - orfs.csv  : flat table of every ORF found (one row per ORF)
    - Summary statistics printed to the terminal
"""

import argparse
import csv
import os
import sys
from pprint import pprint

from src.input_lib.input_validate import run as validate_run
from src.orf_finder_lib.orf_finder import find_orfs, CSV_FIELDNAMES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print_summary(nested: dict) -> None:
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

    print("\n========== ORF Summary ==========")
    print(f"  Total ORFs found            : {total}")
    print(f"  Complete   (ATG)            : {n_complete_canonical}")
    print(f"  Incomplete (ATG)            : {n_incomplete_canonical}")
    print(f"  Complete   (non-canonical)  : {n_complete_noncanonical}")
    print(f"  Incomplete (non-canonical)  : {n_incomplete_noncanonical}")

    # Per non-canonical start codon breakdown
    for sc in ("GTG", "TTG"):
        nc = len(complete["noncanonical"].get(sc, {}))
        ni = len(incomplete["noncanonical"].get(sc, {}))
        if nc + ni > 0:
            print(f"    {sc} — complete: {nc}, incomplete: {ni}")

    print("=================================\n")


def _write_csv(flat_list: list, output_path: str) -> None:
    """Write the flat ORF list to a CSV file."""
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(flat_list)
    print(f"[INFO] ORF table written to: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ORF analysis pipeline — fetches a sequence from NCBI "
                    "and reports all ORFs as a CSV file."
    )
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
    parser.add_argument(
        "--output",
        type=str,
        default="orfs.csv",
        help="Path for the output CSV file (default: orfs.csv)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print the full nested ORF dictionary to stdout",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Get accession and email (prompt if not supplied as flags)
    # ------------------------------------------------------------------
    accession = args.accession or input("Enter NCBI accession number: ").strip()
    email     = args.email     or input("Enter your email (required by NCBI): ").strip()

    # ------------------------------------------------------------------
    # 2. Fetch and validate the sequence
    # ------------------------------------------------------------------
    print(f"\n[INFO] Fetching sequence for accession: {accession}")
    acc, clean_seq = validate_run(accession, email)

    if clean_seq is None:
        print("[ERROR] Pipeline failed: could not retrieve a valid sequence.")
        sys.exit(1)

    print(f"[INFO] Sequence retrieved — {len(clean_seq):,} nucleotides")

    # ------------------------------------------------------------------
    # 3. Run the ORF finder
    # ------------------------------------------------------------------
    print("[INFO] Running ORF finder across all 3 forward reading frames...")
    nested, flat_list = find_orfs(clean_seq)

    # ------------------------------------------------------------------
    # 4. Report results
    # ------------------------------------------------------------------
    _print_summary(nested)

    if args.verbose:
        print("[DEBUG] Full nested ORF dictionary:")
        pprint(nested)

    if not flat_list:
        print("[WARNING] No ORFs found. No CSV written.")
        return

    _write_csv(flat_list, args.output)


if __name__ == "__main__":
    main()

