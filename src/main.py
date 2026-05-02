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
    output including an RSCU codon-usage heatmap.

Input:
    DNA sequence(s) fetched from NCBI via accession number(s),
    OR loaded from local single-sequence FASTA files via --fasta / --fasta2.

Output (single-sequence mode):
    output/cleaned_sequence_1.fasta          : cleaned sequence
    output/<accession>.gff3                  : ORF annotations in GFF3 format
    output/orf_summary.txt                   : human-readable summary report
    output/orf_map.png                       : six-frame ORF map

Output (comparative mode, additional files):
    output/comp_cleaned_sequence_1.fasta     : cleaned sequence 1
    output/comp_cleaned_sequence_2.fasta     : cleaned sequence 2
    output/<accession1>.gff3                 : ORF annotations for sequence 1
    output/<accession2>.gff3                 : ORF annotations for sequence 2
    output/orf_comparison_report.txt         : side-by-side ORF report with codon usage
    output/codon_usage_comparison.png        : RSCU heatmap
"""

import argparse
import os
import sys

from src.input_validate import run as validate_run, validate_start_codons, validate_email
from src.graphics import (
    plot_orf_map,
    plot_comparative_orf_map,
    plot_codon_usage_comparison,
)
from src.orf_finder_lib.orf_finder import find_orfs
from src.analysis_lib.orf_analysis import calculate_orf_stats, find_repeated_orfs
from src.analysis_lib.statistics_summary import (
    write_stats_to_file,
    write_orf_comparison_report,
    write_gff3,
    print_summary,
)


class ORCAPipeline:
    """
    Orchestrates the full ORCA pipeline for one or two DNA sequences.

    Shared run parameters (email, output directory, start codons, minimum
    ORF length) are stored as instance attributes so they do not need to
    be threaded through every method call.  Per-sequence data is passed
    explicitly between methods.

    Typical usage
    -------------
    pipeline = ORCAPipeline(email=email, outdir=outdir,
                            start_codons=start_codons, min_length=min_length)
    pipeline.run(accession=accession, fasta_file=fasta_file,
                 accession2=accession2, fasta_file2=fasta_file2)
    """

    def __init__(
        self,
        email:        str,
        outdir:       str,
        start_codons: list[str],
        min_length:   int,
    ) -> None:
        """
        Initialise the pipeline and create the output directory.

        Parameters
        ----------
        email : str
            User email address required by NCBI Entrez.
        outdir : str
            Directory for all output files.  Created if it does not exist.
        start_codons : list[str]
            Validated start codons to search for, e.g. ``['ATG']``.
        min_length : int
            Minimum ORF length in nucleotides.
        """
        self.email        = email
        self.outdir       = outdir
        self.start_codons = start_codons
        self.min_length   = min_length
        os.makedirs(outdir, exist_ok=True)

    def process_sequence(
        self,
        accession:   str,
        fasta_file:  str | None,
        seq_num:     int,
        comparative: bool,
        label:       str,
    ) -> tuple[str, str, list, list] | tuple[None, None, None, None]:
        """
        Load, validate, and find ORFs for one sequence.

        Parameters
        ----------
        accession : str
            NCBI accession number.  Ignored when *fasta_file* is provided.
        fasta_file : str or None
            Path to a local FASTA file.  Overrides *accession* when given.
        seq_num : int
            Sequence number (1 or 2) used in the cleaned FASTA filename.
        comparative : bool
            When True the ``comp_`` filename prefix is applied.
        label : str
            Human-friendly label shown in terminal output.

        Returns
        -------
        tuple[str, str, list, list]
            ``(accession, clean_seq, nested_dict, flat_list)`` on success.
        tuple[None, None, None, None]
            On failure.
        """
        acc, clean_seq = validate_run(
            accession,
            self.email,
            outdir=self.outdir,
            fasta_file=fasta_file,
            comparative=comparative,
            seq_num=seq_num,
        )
        if clean_seq is None:
            src = f"file '{fasta_file}'" if fasta_file else f"accession '{accession}'"
            print(f"[ERROR] Pipeline failed for {src}.")
            return None, None, None, None

        nested, flat_list = find_orfs(
            clean_seq,
            start_codons=self.start_codons,
            min_length=self.min_length,
        )
        print_summary(nested, flat_list, label=label or acc)

        if not flat_list:
            src = f"file '{fasta_file}'" if fasta_file else f"accession '{accession}'"
            print(f"[WARNING] No ORFs found for {src}. No output files written.")

        return acc, clean_seq, nested, flat_list

    def plot(
        self,
        acc1:        str,
        seq1:        str,
        flat1:       list,
        comparative: bool,
        acc2:        str | None = None,
        seq2:        str | None = None,
        flat2:       list | None = None,
    ) -> None:
        """Generate the ORF map (single or comparative) and, in comparative
        mode, the RSCU codon-usage heatmap."""
        if comparative:
            plot_comparative_orf_map(
                flat1=flat1, seq_len1=len(seq1), acc1=acc1,
                flat2=flat2, seq_len2=len(seq2), acc2=acc2,
                output_path=os.path.join(self.outdir, "orf_map.png"),
            )
            plot_codon_usage_comparison(
                flat1=flat1, acc1=acc1, seq1=seq1,
                flat2=flat2, acc2=acc2, seq2=seq2,
                output_path=os.path.join(self.outdir, "codon_usage_comparison.png"),
            )
        else:
            plot_orf_map(
                flat_list=flat1, seq_len=len(seq1),
                accession=acc1,
                output_path=os.path.join(self.outdir, "orf_map.png"),
            )

    def write_reports(
        self,
        acc1:        str,
        seq1:        str,
        flat1:       list,
        comparative: bool,
        acc2:        str | None = None,
        seq2:        str | None = None,
        flat2:       list | None = None,
    ) -> None:
        """Write GFF3 annotation files and text reports."""
        calculate_orf_stats(flat1, seq1)
        repeats1 = find_repeated_orfs(flat1)
        if repeats1:
            print(f"[INFO] Repeated ORF sequences in {acc1}: {len(repeats1)}")

        write_gff3(flat1, acc1, len(seq1), outdir=self.outdir,
                   start_codons=self.start_codons, min_length=self.min_length)

        if comparative:
            calculate_orf_stats(flat2, seq2)
            repeats2 = find_repeated_orfs(flat2)
            if repeats2:
                print(f"[INFO] Repeated ORF sequences in {acc2}: {len(repeats2)}")

            write_gff3(flat2, acc2, len(seq2), outdir=self.outdir,
                       start_codons=self.start_codons, min_length=self.min_length)

            write_orf_comparison_report(
                flat1=flat1, flat2=flat2,
                acc1=acc1,   acc2=acc2,
                filename=os.path.join(self.outdir, "orf_comparison_report.txt"),
                start_codons=self.start_codons,
                min_length=self.min_length,
            )
        else:
            write_stats_to_file(
                flat1,
                filename=os.path.join(self.outdir, "orf_summary.txt"),
                accession=acc1,
                start_codons=self.start_codons,
                min_length=self.min_length,
            )

    # Public entry point
    def run(
        self,
        accession:   str,
        fasta_file:  str | None = None,
        accession2:  str | None = None,
        fasta_file2: str | None = None,
    ) -> None:
        """
        Run the full ORCA pipeline for one or two sequences.

        Parameters
        ----------
        accession : str
            NCBI accession number for sequence 1.  Ignored when
            *fasta_file* is provided.
        fasta_file : str or None, optional
            Local FASTA file for sequence 1.
        accession2 : str or None, optional
            NCBI accession number for sequence 2.  Triggers comparative
            mode when provided.
        fasta_file2 : str or None, optional
            Local FASTA file for sequence 2.  Triggers comparative mode
            when provided.
        """
        comparative = (accession2 is not None) or (fasta_file2 is not None)

        src1_label = f"(local file) {fasta_file}" if fasta_file else accession
        print(f"\n[ORCA] Processing sequence 1: {src1_label}")

        acc1, seq1, _, flat1 = self.process_sequence(
            accession=accession or "",
            fasta_file=fasta_file,
            seq_num=1,
            comparative=comparative,
            label="Sequence 1",
        )
        if acc1 is None:
            print("[ERROR] Pipeline failed: could not retrieve a valid sequence.")
            sys.exit(1)

        if not flat1:
            print("[ERROR] No ORFs found for sequence 1. Cannot produce output. "
                  "Try lowering --min-length or adding non-canonical start codons.")
            sys.exit(1)

        acc2 = seq2 = flat2 = None
        if comparative:
            src2_label = f"(local file) {fasta_file2}" if fasta_file2 else accession2
            print(f"\n[ORCA] Processing sequence 2: {src2_label}")

            acc2, seq2, _, flat2 = self.process_sequence(
                accession=accession2 or "",
                fasta_file=fasta_file2,
                seq_num=2,
                comparative=True,
                label="Sequence 2",
            )
            if acc2 is None:
                print("[ERROR] Pipeline failed for sequence 2.")
                sys.exit(1)

            if not flat2:
                print("[ERROR] No ORFs found for sequence 2. Cannot produce comparative output. "
                      "Try lowering --min-length or adding non-canonical start codons.")
                sys.exit(1)

        self.plot(acc1, seq1, flat1, comparative, acc2, seq2, flat2)
        self.write_reports(acc1, seq1, flat1, comparative, acc2, seq2, flat2)


# CLI entry point
def main() -> None:
    """Parse command-line arguments and run the ORCA pipeline."""
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

    seq1_group = parser.add_mutually_exclusive_group()
    seq1_group.add_argument(
        "--accession", type=str,
        help="NCBI accession number for sequence 1 (e.g. NM_001301717)",
    )
    seq1_group.add_argument(
        "--fasta", type=str, metavar="FILE",
        help=(
            "Path to a local FASTA file for sequence 1. "
            "The file must contain exactly one sequence. "
            "Cannot be used together with --accession."
        ),
    )

    seq2_group = parser.add_mutually_exclusive_group()
    seq2_group.add_argument(
        "--accession2", type=str, default=None,
        help=(
            "NCBI accession number for sequence 2. "
            "When provided, comparative analysis is performed between the two sequences."
        ),
    )
    seq2_group.add_argument(
        "--fasta2", type=str, default=None, metavar="FILE",
        help=(
            "Path to a local FASTA file for sequence 2. "
            "The file must contain exactly one sequence. "
            "When provided, comparative analysis is performed. "
            "Cannot be used together with --accession2."
        ),
    )

    parser.add_argument("--email", type=str, help="Email address required by NCBI Entrez")
    parser.add_argument("--min-length", type=int, default=30, metavar="NT",
                        help="Minimum ORF length in nucleotides")
    parser.add_argument("--start-codons", nargs="+", default=["ATG"], metavar="CODON",
                        help="Start codons to search for. Default: ATG")
    parser.add_argument("--outdir", type=str, default="output", metavar="DIR",
                        help="Directory for all output files (default: output/)")

    args = parser.parse_args()

    fasta_file = args.fasta
    accession  = args.accession

    if fasta_file is None and accession is None:
        print("\nNo sequence 1 input provided. Choose an input method:")
        print("  1) NCBI accession number")
        print("  2) Local FASTA file path")
        choice = input("Enter 1 or 2: ").strip()
        if choice == "2":
            fasta_file = input("Enter path to FASTA file for sequence 1: ").strip()
            accession  = ""
        else:
            accession = input("Enter NCBI accession number: ").strip()

    email = args.email or input("Enter your email (required by NCBI): ").strip()

    if not validate_email(email):
        sys.exit(1)

    try:
        start_codons = validate_start_codons(args.start_codons)
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    if args.min_length < 3:
        print("[ERROR] --min-length must be at least 3 (one codon).")
        sys.exit(1)

    pipeline = ORCAPipeline(
        email=email,
        outdir=args.outdir,
        start_codons=start_codons,
        min_length=args.min_length,
    )
    pipeline.run(
        accession=accession or "",
        fasta_file=fasta_file,
        accession2=args.accession2,
        fasta_file2=args.fasta2,
    )


if __name__ == "__main__":
    main()
