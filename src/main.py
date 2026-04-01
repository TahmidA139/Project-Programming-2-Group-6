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
    output/cleaned_sequence.fasta            : cleaned sequence
    output/orf_stats.csv                     : per-ORF statistics
    output/stats_summary.txt                 : human-readable summary report
"""

import argparse
import csv
import os
import sys

from src.input_lib.input_validate import run as validate_run, validate_start_codons
from src.orf_finder_lib.orf_finder import find_orfs, CSV_FIELDNAMES
from src.orf_finder_lib.output_writer import write_combined_csv, print_summary
from src.graphics_lib.graphics import plot_orf_map, plot_comparative_orf_map

VALID_START_CODONS = {"ATG", "GTG", "TTG"}

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
    print_summary(nested, flat_list, label=label or accession)

    if not flat_list:
        print(f"[WARNING] No ORFs found for '{accession}'. No output files written.")
        return acc, clean_seq, nested, flat_list


    return acc, clean_seq, nested, flat_list


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
    start_codons = validate_start_codons(args.start_codons)

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
            
    write_combined_csv(
            acc1=acc1, flat1=flat1, seq1=seq1,
            output_path=args.output,
            acc2=acc2 if comparative else None,
            flat2=flat2 if comparative else None,
            seq2=seq2 if comparative else None,
        )
   
    # ── 5. ORF map ────────────────────────────────────────────────────────
    if comparative:
        plot_comparative_orf_map(
            flat1=flat1, seq_len1=len(seq1), acc1=acc1,
            flat2=flat2, seq_len2=len(seq2), acc2=acc2,
            output_path="output/orf_map.png",
        )
    else:
        plot_orf_map(
            flat_list=flat1, seq_len=len(seq1),
            accession=acc1, output_path="output/orf_map.png",
        )
        
if __name__ == "__main__":
    main()

