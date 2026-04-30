#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistics_summary.py

Amanda Yaworsky's part

Purpose:
    All file writing and reporting functions for the ORCA pipeline.
    Computation functions (gc_content, codon_usage, etc.) live in orf_analysis.py.

Exports:
    print_summary               -- Console summary of ORF counts.
    write_stats_to_file         -- Single-sequence orf_summary.txt.
    write_orf_comparison_report -- Combined report for comparative mode only.
    write_combined_csv          -- ORF table CSV (single or comparative mode).
"""

from __future__ import annotations

import csv
import os
from typing import Any, Dict, List, Optional, Tuple

from src.analysis_lib.orf_analysis import codon_usage, gc_content, protein_length, global_alignment_stats
from src.orf_finder_lib.frame_scanner import extract_orf_sequence
from src.orf_finder_lib.orf_finder import CSV_FIELDNAMES

OUTPUT_FIELDNAMES: List[str] = CSV_FIELDNAMES + ["sequence (5'->3')"]


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

def print_summary(nested: dict, flat_list: list, label: str = "") -> None:
    """Print a short summary of ORF counts to stdout."""
    canonical    = nested["canonical"]
    noncanonical = nested["noncanonical"]

    n_canonical    = len(canonical)
    n_noncanonical = sum(len(v) for v in noncanonical.values())
    total          = n_canonical + n_noncanonical
    plus_strand    = sum(1 for o in flat_list if o.get("strand") == "+")
    minus_strand   = sum(1 for o in flat_list if o.get("strand") == "-")

    header = f" ORF Summary{' — ' + label if label else ''} "
    print(f"\n{'-' * 10}{header}{'-' * 10}")
    print(f"  Total ORFs found            : {total}")
    print(f"  Forward strand (+)          : {plus_strand}")
    print(f"  Reverse strand (-)          : {minus_strand}")
    print(f"  Canonical   (ATG)           : {n_canonical}")
    if n_noncanonical > 0:
        print(f"  Non-canonical               : {n_noncanonical}")
        for sc in ("GTG", "TTG"):
            n = len(noncanonical.get(sc, {}))
            if n > 0:
                print(f"    {sc}                       : {n}")
    print("-" * (20 + len(header)))


# ---------------------------------------------------------------------------
# Single-sequence summary (unchanged)
# ---------------------------------------------------------------------------

def write_stats_to_file(
    flat_list: List[Dict[str, Any]],
    filename:  str = "output/orf_summary.txt",
) -> None:
    """
    Write a human-readable summary report for one sequence's ORF set.

    Sections
    --------
    1. Dataset-level stats (total ORFs, average GC)
    2. Longest ORF details
    3. Per-ORF table (index, length, GC%, protein length, strand, frame)
    4. Aggregate codon-usage table

    Parameters
    ----------
    flat_list : list of dict
        Flat ORF list enriched by calculate_orf_stats().
    filename : str
        Output file path.
    """
    total_orfs = len(flat_list)
    avg_gc     = (
        sum(o["gc_content"] for o in flat_list) / total_orfs
        if total_orfs else 0.0
    )
    longest: Optional[Dict[str, Any]] = (
        max(flat_list, key=lambda x: len(x.get("sequence", "")))
        if flat_list else None
    )

    with open(filename, "w") as fh:
        fh.write("=== ORF SUMMARY REPORT ===\n\n")
        fh.write(f"Total ORFs        : {total_orfs}\n")
        fh.write(f"Average GC Content: {avg_gc:.2f}%\n")

        if longest:
            fh.write("\nLongest ORF:\n")
            fh.write(f"  orf_id        : {longest.get('orf_id', 'N/A')}\n")
            fh.write(f"  Length (nt)   : {len(longest['sequence'])}\n")
            fh.write(f"  Strand        : {longest.get('strand', '?')}\n")
            fh.write(f"  Frame         : {longest.get('frame', '?')}\n")
            fh.write(f"  GC Content    : {longest['gc_content']:.2f}%\n")

        fh.write("\n--- Per-ORF Stats ---\n")
        header = (
            f"{'#':<5}{'orf_id':<12}{'Length':>8}"
            f"{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n"
        )
        fh.write(header)
        fh.write("-" * len(header.rstrip()) + "\n")

        for i, orf in enumerate(flat_list):
            fh.write(
                f"{i:<5}{orf.get('orf_id', ''):<12}"
                f"{len(orf.get('sequence', ''))!s:>8}"
                f"{orf['gc_content']:>8.2f}"
                f"{orf['protein_length']:>10}"
                f"{orf.get('strand', '?'):>8}"
                f"{orf.get('frame', '?')!s:>7}\n"
            )

        fh.write("\n--- Codon Usage (aggregate) ---\n")
        total_codons: Dict[str, int] = {}
        for orf in flat_list:
            for codon, count in codon_usage(orf.get("sequence", "")).items():
                total_codons[codon] = total_codons.get(codon, 0) + count
        for codon, count in sorted(total_codons.items()):
            fh.write(f"  {codon}: {count}\n")


# ---------------------------------------------------------------------------
# Private helpers for comparative report
# ---------------------------------------------------------------------------

def _write_sequence_section(
    fh,
    flat_list: List[Dict[str, Any]],
    label: str,
) -> None:
    """Write dataset summary and per-ORF table for one sequence."""
    total_orfs = len(flat_list)
    avg_gc     = (
        sum(o["gc_content"] for o in flat_list) / total_orfs
        if total_orfs else 0.0
    )
    longest: Optional[Dict[str, Any]] = (
        max(flat_list, key=lambda x: len(x.get("sequence", "")))
        if flat_list else None
    )

    fh.write(f"=== {label} ===\n\n")
    fh.write(f"Total ORFs        : {total_orfs}\n")
    fh.write(f"Average GC Content: {avg_gc:.2f}%\n")

    if longest:
        fh.write("\nLongest ORF:\n")
        fh.write(f"  orf_id        : {longest.get('orf_id', 'N/A')}\n")
        fh.write(f"  Length (nt)   : {len(longest['sequence'])}\n")
        fh.write(f"  Strand        : {longest.get('strand', '?')}\n")
        fh.write(f"  Frame         : {longest.get('frame', '?')}\n")
        fh.write(f"  GC Content    : {longest['gc_content']:.2f}%\n")

    fh.write("\n--- Per-ORF Stats ---\n")
    header = (
        f"{'#':<5}{'orf_id':<12}{'Length':>8}"
        f"{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n"
    )
    fh.write(header)
    fh.write("-" * len(header.rstrip()) + "\n")

    for i, orf in enumerate(flat_list):
        fh.write(
            f"{i:<5}{orf.get('orf_id', ''):<12}"
            f"{len(orf.get('sequence', ''))!s:>8}"
            f"{orf['gc_content']:>8.2f}"
            f"{orf['protein_length']:>10}"
            f"{orf.get('strand', '?'):>8}"
            f"{orf.get('frame', '?')!s:>7}\n"
        )


def _write_comparative_summary(
    fh,
    flat1: List[Dict[str, Any]],
    flat2: List[Dict[str, Any]],
    acc1:  str,
    acc2:  str,
) -> None:
    """Write the shared/unique ORF comparison block."""
    plus1  = sum(1 for o in flat1 if o.get("strand") == "+")
    minus1 = sum(1 for o in flat1 if o.get("strand") == "-")
    plus2  = sum(1 for o in flat2 if o.get("strand") == "+")
    minus2 = sum(1 for o in flat2 if o.get("strand") == "-")

    seqs1  = {o.get("sequence", "") for o in flat1}
    seqs2  = {o.get("sequence", "") for o in flat2}
    shared = seqs1 & seqs2

    fh.write("=== Comparative Summary ===\n\n")
    fh.write(f"{'Metric':<30} {acc1:>20} {acc2:>20}\n")
    fh.write("-" * 72 + "\n")
    fh.write(f"{'Total ORFs':<30} {len(flat1):>20} {len(flat2):>20}\n")
    fh.write(f"{'Forward strand (+)':<30} {plus1:>20} {plus2:>20}\n")
    fh.write(f"{'Reverse strand (-)':<30} {minus1:>20} {minus2:>20}\n")
    fh.write("\n")
    fh.write(f"Shared ORF sequences       : {len(shared)}\n")
    fh.write(f"Unique to {acc1:<22}: {len(seqs1 - seqs2)}\n")
    fh.write(f"Unique to {acc2:<22}: {len(seqs2 - seqs1)}\n")


def _write_alignment_section(
    fh,
    seq1: str,
    seq2: str,
    acc1: str,
    acc2: str,
) -> None:
    """Compute global alignment and write a summary block to an open file handle."""
    fh.write("=== Global Pairwise Alignment ===\n\n")
    fh.write(
        "  Algorithm : Needleman-Wunsch (global)\n"
        "  Scoring   : match +1 | mismatch 0 | gap open -2 | gap extend -0.5\n\n"
    )

    stats = global_alignment_stats(seq1, seq2)

    fh.write(f"  {acc1:<30}: {stats['seq1_len']:,} bp\n")
    fh.write(f"  {acc2:<30}: {stats['seq2_len']:,} bp\n")
    fh.write(f"  Alignment length            : {stats['alignment_length']:,} bp\n")
    fh.write(f"  Matches                     : {stats['matches']:,}\n")
    fh.write(f"  Mismatches                  : {stats['mismatches']:,}\n")
    fh.write(f"  Gaps                        : {stats['gaps']:,}\n")
    fh.write(f"  Sequence identity           : {stats['identity_pct']:.2f}%\n")
    fh.write(f"  Coverage (vs longer seq)    : {stats['coverage_pct']:.2f}%\n")
    fh.write(f"  Raw alignment score         : {stats['score']:.1f}\n")

    identity = stats["identity_pct"]
    if identity >= 95:
        interp = "highly conserved (≥95% identity)"
    elif identity >= 70:
        interp = "moderately conserved (70–94% identity)"
    elif identity >= 40:
        interp = "distantly related (40–69% identity)"
    else:
        interp = "low similarity (<40% identity) — sequences may be unrelated"

    fh.write(f"\n  Interpretation: {interp}\n")


# ---------------------------------------------------------------------------
# Comparative report (comparative mode only)
# ---------------------------------------------------------------------------

def write_orf_comparison_report(
    flat1:    List[Dict[str, Any]],
    flat2:    List[Dict[str, Any]],
    acc1:     str,
    acc2:     str,
    seq1:     str,
    seq2:     str,
    filename: str = "output/orf_comparison_report.txt",
) -> None:
    """
    Write a single combined report for a two-sequence comparative run.

    Sections
    --------
    1. Per-sequence ORF summaries (dataset stats, longest ORF, per-ORF table).
    2. Comparative ORF summary (shared / unique ORF sequences, strand counts).
    3. Global pairwise alignment (Needleman-Wunsch identity and coverage stats).

    Parameters
    ----------
    flat1 : list of dict
        Enriched ORF list for sequence 1 (requires calculate_orf_stats()).
    flat2 : list of dict
        Enriched ORF list for sequence 2.
    acc1 : str
        Accession or label for sequence 1.
    acc2 : str
        Accession or label for sequence 2.
    seq1 : str
        Cleaned DNA sequence for sequence 1; used for global alignment.
    seq2 : str
        Cleaned DNA sequence for sequence 2; used for global alignment.
    filename : str
        Output file path.
    """
    with open(filename, "w") as fh:
        fh.write("=== ORF COMPARISON REPORT ===\n\n")
        _write_sequence_section(fh, flat1, f"Sequence 1: {acc1}")
        fh.write("\n\n")
        _write_sequence_section(fh, flat2, f"Sequence 2: {acc2}")
        fh.write("\n\n")
        _write_comparative_summary(fh, flat1, flat2, acc1, acc2)
        fh.write("\n\n")
        _write_alignment_section(fh, seq1, seq2, acc1, acc2)


# ---------------------------------------------------------------------------
# Combined CSV output
# ---------------------------------------------------------------------------

def _write_sequence_block(
    fh,
    writer: csv.DictWriter,
    accession: str,
    flat_list: list,
    dna_sequence: str,
) -> None:
    """Write one accession header + ORF rows into an already-open CSV file."""
    fh.write(f"{accession}\n")
    writer.writeheader()
    for orf in flat_list:
        row = dict(orf)
        row["sequence (5'->3')"] = extract_orf_sequence(orf, dna_sequence)
        writer.writerow(row)


def write_combined_csv(
    acc1:        str,
    flat1:       list,
    seq1:        str,
    output_path: str,
    acc2:        Optional[str]  = None,
    flat2:       Optional[list] = None,
    seq2:        Optional[str]  = None,
) -> None:
    """Write one or two ORF tables into a single CSV file."""
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=OUTPUT_FIELDNAMES, extrasaction="ignore"
        )
        _write_sequence_block(fh, writer, acc1, flat1, seq1)

        if acc2 is not None and flat2 is not None and seq2 is not None:
            fh.write("\n\n")
            _write_sequence_block(fh, writer, acc2, flat2, seq2)
