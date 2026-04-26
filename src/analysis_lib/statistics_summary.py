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
    write_stats_to_file         -- Single-sequence orf_summary.txt (unchanged).
    write_orf_comparison_report -- Combined report for comparative mode only.
    write_comparative_csv       -- Codon usage CSV for two sequences.
    write_combined_csv          -- ORF table CSV (single or comparative mode).
"""

from __future__ import annotations

import csv
import os
from typing import Any, Dict, List, Optional, Tuple

from src.analysis_lib.orf_analysis import codon_usage
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
    codon_csv_name: str,
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
    fh.write(f"\nSee {codon_csv_name} for full codon usage breakdown.\n")


# ---------------------------------------------------------------------------
# Comparative report (comparative mode only)
# ---------------------------------------------------------------------------

def write_orf_comparison_report(
    flat1:          List[Dict[str, Any]],
    flat2:          List[Dict[str, Any]],
    acc1:           str,
    acc2:           str,
    filename:       str = "output/orf_comparison_report.txt",
    codon_csv_name: str = "codon_comparison.csv",
) -> None:
    """
    Write a single combined report for a two-sequence comparative run.

    Replaces the three separate files produced in comparative mode
    (orf_summary.txt, orf_summary_seq2.txt, comparison.txt) with one file.
    Single-sequence mode is unaffected — it still uses write_stats_to_file.

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
    filename : str
        Output file path.
    codon_csv_name : str
        Base filename of the codon comparison CSV, referenced at the end
        of the comparative summary section.
    """
    with open(filename, "w") as fh:
        fh.write("=== ORF COMPARISON REPORT ===\n\n")
        _write_sequence_section(fh, flat1, f"Sequence 1: {acc1}")
        fh.write("\n\n")
        _write_sequence_section(fh, flat2, f"Sequence 2: {acc2}")
        fh.write("\n\n")
        _write_comparative_summary(fh, flat1, flat2, acc1, acc2, codon_csv_name)


# ---------------------------------------------------------------------------
# Comparative codon usage CSV
# ---------------------------------------------------------------------------

def write_comparative_csv(
    flat1:    List[Dict[str, Any]],
    flat2:    List[Dict[str, Any]],
    acc1:     str = "Sequence 1",
    acc2:     str = "Sequence 2",
    filename: str = "output/codon_comparison.csv",
) -> None:
    """
    Write a CSV comparing codon-usage frequencies between two ORF sets.

    Columns: Codon, <acc1>_count, <acc2>_count, Delta (acc1 - acc2)

    Requires calculate_orf_stats() to have been called on both lists.
    """
    def total_codons(orfs: List[Dict[str, Any]]) -> Dict[str, int]:
        totals: Dict[str, int] = {}
        for orf in orfs:
            for codon, count in codon_usage(orf.get("sequence", "")).items():
                totals[codon] = totals.get(codon, 0) + count
        return totals

    codons1    = total_codons(flat1)
    codons2    = total_codons(flat2)
    all_codons = sorted(set(codons1) | set(codons2))

    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Codon", f"{acc1}_count", f"{acc2}_count", "Delta"])
        for codon in all_codons:
            c1 = codons1.get(codon, 0)
            c2 = codons2.get(codon, 0)
            writer.writerow([codon, c1, c2, c1 - c2])


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
