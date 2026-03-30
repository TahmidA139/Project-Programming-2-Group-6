#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py

Purpose:
    Main driver for the ORCA (ORF Recognition and Comparative Analysis)
    pipeline.

Role in Project:
    Calls other modules to handle sequence loading, ORF detection,
    statistics, and CSV/text reporting.  When --accession2 is supplied the
    pipeline runs on both sequences and produces comparative output.

Input:
    DNA sequence(s) fetched from NCBI via accession number(s).

Output (single-sequence mode):
    output/orfs.csv            : flat table of every ORF found
    output/orf_stats.csv       : per-ORF statistics
    output/stats_summary.txt   : human-readable summary report

Output (comparative mode, additional files):
    output/orfs_seq2.csv           : ORF table for the second sequence
    output/orf_stats_seq2.csv      : per-ORF statistics for the second sequence
    output/stats_summary_seq2.txt  : human-readable summary for the second sequence
    output/comparative_report.txt  : side-by-side comparative report
    output/comparative_stats.csv   : codon-usage delta table

Example usage
-------------
# Single sequence — all defaults (ATG only, min length 30, nested included)
python main.py --accession NM_001301717 --email you@example.com

# Comparative — two accessions
python main.py --accession NM_001301717 --accession2 NM_001256799 \\
    --email you@example.com

# All three start codons, minimum 60 nt, ignore nested ORFs
python main.py --accession NM_001301717 --email you@example.com \\
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
    write_comparative_report,
    write_comparative_csv,
    ORF_STATS_FIELDNAMES,
)
from src.analysis_lib.orf_analysis import compare_orf_sets

# Valid start codons the user is allowed to request
VALID_START_CODONS = {"ATG", "GTG", "TTG"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_nested(flat_list: list) -> list:
    """Return the subset of ORFs that are nested inside another ORF."""
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


def _print_summary(nested: dict, flat_list: list, label: str = "") -> None:
    """Print a short summary of ORF counts to stdout."""
    complete   = nested["complete"]
    incomplete = nested["incomplete"]

    n_complete_canonical      = len(complete["canonical"])
    n_incomplete_canonical    = len(incomplete["canonical"])
    n_complete_noncanonical   = sum(len(v) for v in complete["noncanonical"].values())
    n_incomplete_noncanonical = sum(len(v) for v in incomplete["noncanonical"].values())

    total        = (n_complete_canonical + n_incomplete_canonical
                    + n_complete_noncanonical + n_incomplete_noncanonical)
    plus_strand  = sum(1 for o in flat_list if o.get("strand") == "+")
    minus_strand = sum(1 for o in flat_list if o.get("strand") == "-")

    header = f" ORF Summary{' — ' + label if label else ''} "
    print(f"\n{'=' * 10}{header}{'=' * 10}")
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
    print("=" * (20 + len(header)))


def _write_csv(flat_list: list, output_path: str) -> None:
    """Write the flat ORF list to a CSV file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True) \
        if os.path.dirname(output_path) else None
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat_list)
    print(f"[INFO] ORF table written to: {output_path}")


def _write_csv_with_fields(data: list, output_path: str, fieldnames: list) -> None:
    """Write a list of dicts to a CSV file with specified fieldnames."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True) \
        if os.path.dirname(output_path) else None
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


def _run_single_sequence(
    accession:     str,
    email:         str,
    start_codons:  list,
    min_length:    int,
    ignore_nested: bool,
    orf_csv:       str,
    stats_csv:     str,
    summary_txt:   str,
    label:         str = "",
) -> tuple[str, str, list, list] | tuple[None, None, None, None]:
    """
    Fetch, validate, analyse, and write outputs for a single accession.

    Returns
    -------
    (accession, clean_seq, nested_dict, flat_list) on success,
    (None, None, None, None) on failure.
    """
    # 1. Fetch and validate
    acc, clean_seq = validate_run(accession, email)
    if clean_seq is None:
        print(f"[ERROR] Pipeline failed for accession '{accession}'.")
        return None, None, None, None

    # 2. Find ORFs
    nested, flat_list = find_orfs(
        clean_seq,
        start_codons=start_codons,
        min_length=min_length,
        ignore_nested=ignore_nested,
    )

    # 3. Print terminal summary
    _print_summary(nested, flat_list, label=label or accession)

    if not flat_list:
        print(f"[WARNING] No ORFs found for '{accession}'. No output files written.")
        return acc, clean_seq, nested, flat_list

    # 4. Write output files
    _write_csv(flat_list, orf_csv)

    per_orf_stats = calculate_orf_stats(flat_list, clean_seq)
    _write_csv_with_fields(per_orf_stats, stats_csv, ORF_STATS_FIELDNAMES)

    write_stats_to_file(flat_list, clean_seq, acc, outfile=summary_txt)

    return acc, clean_seq, nested, flat_list


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ORCA — ORF Recognition and Comparative Analysis.\n"
            "Fetches sequence(s) from NCBI and reports all ORFs as CSV files.\n"
            "Supply --accession2 to enable side-by-side comparative analysis."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Sequence retrieval ────────────────────────────────────────────────
    parser.add_argument(
        "--accession",
        type=str,
        help="NCBI accession number for sequence 1 (e.g. NM_001301717)",
    )
    parser.add_argument(
        "--accession2",
        type=str,
        default=None,
        help=(
            "NCBI accession number for sequence 2. "
            "When provided, comparative analysis is performed between "
            "the two sequences."
        ),
    )
    parser.add_argument(
        "--email",
        type=str,
        help="Email address required by NCBI Entrez",
    )

    # ── ORF finder options ────────────────────────────────────────────────
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

    # ── Output options ────────────────────────────────────────────────────
    parser.add_argument(
        "--output",
        type=str,
        default="output/orfs.csv",
        help="Path for the primary ORF output CSV file",
    )

    args = parser.parse_args()

    # ── 1. Collect accession(s) and email ─────────────────────────────────
    accession  = args.accession  or input("Enter NCBI accession number: ").strip()
    accession2 = args.accession2  # None if not supplied
    email      = args.email      or input("Enter your email (required by NCBI): ").strip()

    comparative = accession2 is not None

    # ── 2. Validate shared options ────────────────────────────────────────
    start_codons = _validate_start_codons(args.start_codons)

    if args.min_length < 3:
        print("[ERROR] --min-length must be at least 3 (one codon).")
        sys.exit(1)

    # ── 3. Run pipeline for sequence 1 ───────────────────────────────────
    print(f"\n[ORCA] Processing sequence 1: {accession}")

    acc1, seq1, nested1, flat1 = _run_single_sequence(
        accession     = accession,
        email         = email,
        start_codons  = start_codons,
        min_length    = args.min_length,
        ignore_nested = args.ignore_nested,
        orf_csv       = args.output,
        stats_csv     = "output/orf_stats.csv",
        summary_txt   = "output/stats_summary.txt",
        label         = "Sequence 1",
    )

    if acc1 is None:
        print("[ERROR] Pipeline failed: could not retrieve a valid sequence.")
        sys.exit(1)

    # ── 4. Run pipeline for sequence 2 (if requested) ────────────────────
    if comparative:
        print(f"\n[ORCA] Processing sequence 2: {accession2}")

        acc2, seq2, nested2, flat2 = _run_single_sequence(
            accession     = accession2,
            email         = email,
            start_codons  = start_codons,
            min_length    = args.min_length,
            ignore_nested = args.ignore_nested,
            orf_csv       = "output/orfs_seq2.csv",
            stats_csv     = "output/orf_stats_seq2.csv",
            summary_txt   = "output/stats_summary_seq2.txt",
            label         = "Sequence 2",
        )

        if acc2 is None:
            print("[ERROR] Pipeline failed for sequence 2.")
            sys.exit(1)

        # ── 5. Comparative analysis ───────────────────────────────────────
        print("\n[ORCA] Running comparative analysis…")

        # Enrich flat lists with codon_usage + gc_content_pct before comparing
        enriched1 = calculate_orf_stats(flat1, seq1) if flat1 else []
        enriched2 = calculate_orf_stats(flat2, seq2) if flat2 else []

        comparison = compare_orf_sets(
            flat_list_1 = enriched1,
            flat_list_2 = enriched2,
            acc1        = acc1,
            acc2        = acc2,
            seq1        = seq1,
            seq2        = seq2,
        )

        write_comparative_report(comparison, outfile="output/comparative_report.txt")
        write_comparative_csv(comparison,    outfile="output/comparative_stats.csv")

        # ── 6. Console summary of comparative findings ────────────────────
        print("\n========== Comparative Summary ==========")
        print(f"  {'Metric':<30} {acc1:>14}  {acc2:>14}")
        print("  " + "-" * 62)
        for label, key in [
            ("Total ORFs",         "total_orfs"),
            ("Complete ORFs",      "complete_orfs"),
            ("Incomplete ORFs",    "incomplete_orfs"),
            ("Nested ORFs",        "nested_orfs"),
            ("Forward strand (+)", "plus_strand_orfs"),
            ("Reverse strand (-)", "minus_strand_orfs"),
        ]:
            v1, v2 = comparison[key]
            print(f"  {label:<30} {str(v1):>14}  {str(v2):>14}")

        gc1, gc2 = comparison["genomic_gc"]
        print(f"  {'Genomic GC%':<30} {str(gc1):>14}  {str(gc2):>14}")

        shared_n = len(comparison["shared_start_sites"])
        u1_n     = len(comparison["unique_to_seq1"])
        u2_n     = len(comparison["unique_to_seq2"])
        print(f"\n  Shared ORF start positions  : {shared_n}")
        print(f"  Unique to {acc1:<14}: {u1_n}")
        print(f"  Unique to {acc2:<14}: {u2_n}")
        print("=========================================\n")

    else:
        # Single-sequence mode — nothing extra to do; files already written.
        if not flat1:
            print("[WARNING] No ORFs found. No output files written.")


if __name__ == "__main__":
    main()

