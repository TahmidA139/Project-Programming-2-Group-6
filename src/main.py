#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py

Everyone has contributed to this file

Purpose:
    Main driver for the ORCA (ORF Recognition and Comparative Analysis)
    pipeline.

Role in Project:
    Calls other modules to handle sequence loading, ORF detection,
    statistics, and CSV/text reporting.  When --accession2 / --fasta2 is
    supplied the pipeline runs on both sequences and produces comparative
    output.

Input:
    DNA sequence(s) fetched from NCBI via accession number(s),
    OR loaded from local single-sequence FASTA files via --fasta / --fasta2.

Output (single-sequence mode):
    output/cleaned_sequence_1.fasta          : cleaned sequence
    output/orf_stats.csv                     : per-ORF statistics
    output/stats_summary.txt                 : human-readable summary report

Output (comparative mode, additional files):
    output/comp_cleaned_sequence_1.fasta     : cleaned sequence 1
    output/comp_cleaned_sequence_2.fasta     : cleaned sequence 2
"""

import argparse
import csv
import os
import sys

from src.input_validate import (run as validate_run, validate_start_codons,)
from src.graphics import plot_orf_map, plot_comparative_orf_map, plot_codon_usage_comparison
from src.orf_finder_lib.orf_finder import find_orfs, CSV_FIELDNAMES
from src.analysis_lib.orf_analysis import calculate_orf_stats, find_repeated_orfs
from src.analysis_lib.statistics_summary import (write_stats_to_file, write_orf_comparison_report, write_comparative_csv, write_combined_csv, print_summary,)


def _run_single_sequence(
    accession:     str,
    email:         str,
    start_codons:  list,
    min_length:    int,
    orf_csv:       str,
    stats_csv:     str,
    summary_txt:   str,
    outdir:        str = "output",
    label:         str = "",
    fasta_file:    str | None = None,
    comparative:   bool = False,
    seq_num:       int  = 1,
) -> tuple[str, str, list, list] | tuple[None, None, None, None]:
    """
    Fetch/load, validate, analyse, and write outputs for a single sequence.

    Parameters
    ----------
    accession   : NCBI accession number. Ignored when fasta_file is given.
    email       : User email for NCBI Entrez.
    start_codons: List of start codons to search for.
    min_length  : Minimum ORF length in nucleotides.
    orf_csv     : Output path for the ORF CSV file.
    stats_csv   : Output path for the stats CSV file.
    summary_txt : Output path for the human-readable summary.
    outdir      : Output directory passed through to validate_run so cleaned
                  FASTA files land in the correct location.
    label       : Human-friendly label printed in terminal output.
    fasta_file  : Path to a local FASTA file (mutually exclusive with NCBI).
                  When provided, the sequence is loaded from disk instead of
                  being fetched from NCBI.  The file must contain exactly one
                  sequence — if it contains more, the pipeline exits with an
                  error (handled inside load_fasta_from_file).
    comparative : When True, the cleaned FASTA is written with the comp_
                  prefix (e.g. comp_cleaned_sequence_1.fasta).
    seq_num     : Sequence number (1 or 2) used in the cleaned FASTA filename.

    Returns
    -------
    (accession, clean_seq, nested_dict, flat_list) on success,
    (None, None, None, None) on failure.
    """
    acc, clean_seq, _, _ = validate_run(
        accession,
        email,
        outdir=outdir,
        fasta_file=fasta_file,
        comparative=comparative,
        seq_num=seq_num,
    )
    if clean_seq is None:
        src = f"file '{fasta_file}'" if fasta_file else f"accession '{accession}'"
        print(f"[ERROR] Pipeline failed for {src}.")
        return None, None, None, None

    # 2. Find ORFs
    nested, flat_list = find_orfs(
        clean_seq,
        start_codons=start_codons,
        min_length=min_length,
    )

    # 3. Print terminal summary
    print_summary(nested, flat_list, label=label or acc)

    if not flat_list:
        src = f"file '{fasta_file}'" if fasta_file else f"accession '{accession}'"
        print(f"[WARNING] No ORFs found for {src}. No output files written.")
        return acc, clean_seq, nested, flat_list

    return acc, clean_seq, nested, flat_list


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ORCA — ORF Recognition and Comparative Analysis.\n"
            "Fetches sequence(s) from NCBI or loads them from local FASTA files,\n"
            "then reports all ORFs as CSV files.\n"
            "Supply --accession2 or --fasta2 to enable side-by-side comparative analysis.\n"
            "\n"
            "Input source rules:\n"
            "  Sequence 1: provide EITHER --accession OR --fasta, not both.\n"
            "  Sequence 2: provide EITHER --accession2 OR --fasta2, not both.\n"
            "  Local FASTA files must contain exactly ONE sequence.\n"
            "  Files with multiple sequences will cause the pipeline to exit\n"
            "  with a clear error message."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Sequence retrieval ────────────────────────────────────────────────
    seq1_group = parser.add_mutually_exclusive_group()
    seq1_group.add_argument(
        "--accession",
        type=str,
        help="NCBI accession number for sequence 1 (e.g. NM_001301717)",
    )
    seq1_group.add_argument(
        "--fasta",
        type=str,
        metavar="FILE",
        help=(
            "Path to a local FASTA file for sequence 1. "
            "The file must contain exactly one sequence. "
            "Cannot be used together with --accession."
        ),
    )

    seq2_group = parser.add_mutually_exclusive_group()
    seq2_group.add_argument(
        "--accession2",
        type=str,
        default=None,
        help=(
            "NCBI accession number for sequence 2. "
            "When provided, comparative analysis is performed between "
            "the two sequences."
        ),
    )
    seq2_group.add_argument(
        "--fasta2",
        type=str,
        default=None,
        metavar="FILE",
        help=(
            "Path to a local FASTA file for sequence 2. "
            "The file must contain exactly one sequence. "
            "When provided, comparative analysis is performed. "
            "Cannot be used together with --accession2."
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
        help="Start codons to search for. ATG is canonical; GTG and TTG are non-canonical alternatives. Default: ATG",
    )

    # ── Output options ────────────────────────────────────────────────────
    parser.add_argument(
        "--outdir",
        type=str,
        default="output",
        metavar="DIR",
        help="Directory for all output files (default: output/)",
    )

    args = parser.parse_args()

    # Create the output directory if it does not already exist
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # ── 1. Resolve sequence 1 source ──────────────────────────────────────
    # Determine whether the user is providing an accession or a local file.
    # If neither was passed on the command line, ask interactively.
    fasta_file  = args.fasta       # None if not supplied
    accession   = args.accession   # None if not supplied

    if fasta_file is None and accession is None:
        # Interactive fallback: ask which input method the user wants
        print("\nNo sequence 1 input provided. Choose an input method:")
        print("  1) NCBI accession number")
        print("  2) Local FASTA file path")
        choice = input("Enter 1 or 2: ").strip()
        if choice == "2":
            fasta_file = input("Enter path to FASTA file for sequence 1: ").strip()
            accession  = ""   # placeholder; overridden by the record ID in the file
        else:
            accession  = input("Enter NCBI accession number: ").strip()

    # ── 2. Resolve sequence 2 source ──────────────────────────────────────
    fasta_file2  = args.fasta2      # None if not supplied
    accession2   = args.accession2  # None if not supplied

    comparative  = (accession2 is not None) or (fasta_file2 is not None)

    # ── 3. Get email ──────────────────────────────────────────────────────
    email = args.email or input("Enter your email (required by NCBI): ").strip()

    # ── 4. Validate shared options ────────────────────────────────────────
    start_codons = validate_start_codons(args.start_codons)

    if args.min_length < 3:
        print("[ERROR] --min-length must be at least 3 (one codon).")
        sys.exit(1)

    # ── 5. Run pipeline for sequence 1 ───────────────────────────────────
    src1_label = f"(local file) {fasta_file}" if fasta_file else accession
    print(f"\n[ORCA] Processing sequence 1: {src1_label}")

    acc1, seq1, nested1, flat1 = _run_single_sequence(
        accession     = accession or "",
        email         = email,
        start_codons  = start_codons,
        min_length    = args.min_length,
        orf_csv       = os.path.join(outdir, "orfs.csv"),
        stats_csv     = os.path.join(outdir, "orf_stats.csv"),
        summary_txt   = os.path.join(outdir, "stats_summary.txt"),
        outdir        = outdir,
        label         = "Sequence 1",
        fasta_file    = fasta_file,
        comparative   = comparative,
        seq_num       = 1,
    )

    if acc1 is None:
        print("[ERROR] Pipeline failed: could not retrieve a valid sequence.")
        sys.exit(1)

    # ── 6. Run pipeline for sequence 2 (if requested) ────────────────────
    if comparative:
        src2_label = f"(local file) {fasta_file2}" if fasta_file2 else accession2
        print(f"\n[ORCA] Processing sequence 2: {src2_label}")

        acc2, seq2, nested2, flat2 = _run_single_sequence(
            accession     = accession2 or "",
            email         = email,
            start_codons  = start_codons,
            min_length    = args.min_length,
            orf_csv       = os.path.join(outdir, "orfs_seq2.csv"),
            stats_csv     = os.path.join(outdir, "orf_stats_seq2.csv"),
            summary_txt   = os.path.join(outdir, "stats_summary_seq2.txt"),
            outdir        = outdir,
            label         = "Sequence 2",
            fasta_file    = fasta_file2,
            comparative   = True,
            seq_num       = 2,
        )

        if acc2 is None:
            print("[ERROR] Pipeline failed for sequence 2.")
            sys.exit(1)

        write_combined_csv(
            acc1=acc1, flat1=flat1, seq1=seq1,
            output_path=os.path.join(outdir, "orfs.csv"),
            acc2=acc2, flat2=flat2, seq2=seq2,
        )

    else:
        write_combined_csv(
            acc1=acc1, flat1=flat1, seq1=seq1,
            output_path=os.path.join(outdir, "orfs.csv"),
        )

    # ── 7. ORF map ────────────────────────────────────────────────────────
    if comparative:
        plot_comparative_orf_map(
            flat1=flat1, seq_len1=len(seq1), acc1=acc1,
            flat2=flat2, seq_len2=len(seq2), acc2=acc2,
            output_path=os.path.join(outdir, "orf_map.png"),
        )
        plot_codon_usage_comparison(
            seq1=seq1, acc1=acc1,
            seq2=seq2, acc2=acc2,
            output_path=os.path.join(outdir, "codon_usage_comparison.png"),
        )
    else:
        plot_orf_map(
            flat_list=flat1, seq_len=len(seq1),
            accession=acc1, output_path=os.path.join(outdir, "orf_map.png"),
        )

    # ── 8. Enrich ORFs with sequence/GC/protein stats ────────────────────
    calculate_orf_stats(flat1, seq1)
    repeats1 = find_repeated_orfs(flat1)
    if repeats1:
        print(f"[INFO] Repeated ORF sequences in {acc1}: {len(repeats1)}")

    if comparative:
        calculate_orf_stats(flat2, seq2)
        repeats2 = find_repeated_orfs(flat2)
        if repeats2:
            print(f"[INFO] Repeated ORF sequences in {acc2}: {len(repeats2)}")
        write_orf_comparison_report(
            flat1=flat1, flat2=flat2,
            acc1=acc1,   acc2=acc2,
            filename=os.path.join(outdir, "orf_comparison_report.txt"),
            codon_csv_name="codon_comparison.csv",
        )
        write_comparative_csv(flat1, flat2, acc1=acc1, acc2=acc2,
                              filename=os.path.join(outdir, "codon_comparison.csv"))
    else:
        write_stats_to_file(flat1, filename=os.path.join(outdir, "orf_summary.txt"))

if __name__ == "__main__":
    main()
