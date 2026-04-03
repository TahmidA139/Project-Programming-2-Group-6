#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistics_summary.py

Amanda Yaworsky's part 

Purpose:
    File writing and reporting functions for the ORCA pipeline.
    Computation functions (gc_content, codon_usage, etc.) live in orf_analysis.py.
"""

from __future__ import annotations

import csv
from typing import Any, Dict, List, Optional, Tuple
from src.analysis_lib.orf_analysis import codon_usage


# ---------------------------------------------------------------------------
# Single sequence summary
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

        # Per-ORF table
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

        # Aggregate codon usage
        fh.write("\n--- Codon Usage (aggregate) ---\n")
        total_codons: Dict[str, int] = {}
        for orf in flat_list:
            for codon, count in codon_usage(orf.get("sequence", "")).items():
                total_codons[codon] = total_codons.get(codon, 0) + count
        for codon, count in sorted(total_codons.items()):
            fh.write(f"  {codon}: {count}\n")


# ---------------------------------------------------------------------------
# Comparative report
# ---------------------------------------------------------------------------

def avg_gc(orfs: List[Dict[str, Any]]) -> float:
    """Return average GC content across a list of ORFs."""
    return sum(o["gc_content"] for o in orfs) / len(orfs) if orfs else 0.0


def strand_counts(orfs: List[Dict[str, Any]]) -> Tuple[int, int]:
    """Return (plus_count, minus_count) for a list of ORFs."""
    plus  = sum(1 for o in orfs if o.get("strand") == "+")
    minus = sum(1 for o in orfs if o.get("strand") == "-")
    return plus, minus


def write_comparative_report(
    flat1:    List[Dict[str, Any]],
    flat2:    List[Dict[str, Any]],
    acc1:     str = "Sequence 1",
    acc2:     str = "Sequence 2",
    filename: str = "output/comparison.txt",
) -> None:
    """
    Write a side-by-side text comparison of two ORF sets.
    
    Parameters
    ----------
    flat1 : list of dict
        Enriched ORF list for sequence 1.
    flat2 : list of dict
        Enriched ORF list for sequence 2.
    acc1 : str
        Label for sequence 1.
    acc2 : str
        Label for sequence 2.
    filename : str
        Output file path.
    """
    with open(filename, "w") as fh:
        fh.write("=== COMPARATIVE ORF REPORT ===\n\n")

        for label, orfs in ((acc1, flat1), (acc2, flat2)):
            plus, minus = strand_counts(orfs)
            fh.write(f"[{label}]\n")
            fh.write(f"  Total ORFs      : {len(orfs)}\n")
            fh.write(f"  Average GC      : {avg_gc(orfs):.2f}%\n")
            fh.write(f"  Forward strand  : {plus}\n")
            fh.write(f"  Reverse strand  : {minus}\n\n")

        # Shared / unique sequences
        seqs1  = {o.get("sequence", "") for o in flat1}
        seqs2  = {o.get("sequence", "") for o in flat2}
        shared = seqs1 & seqs2
        fh.write(f"Shared ORF sequences  : {len(shared)}\n")
        fh.write(f"Unique to {acc1:<12}: {len(seqs1 - seqs2)}\n")
        fh.write(f"Unique to {acc2:<12}: {len(seqs2 - seqs1)}\n")

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

    Parameters
    ----------
    flat1 : list of dict
        Enriched ORF list for sequence 1.
    flat2 : list of dict
        Enriched ORF list for sequence 2.
    acc1 : str
        Label for sequence 1 (used as column header prefix).
    acc2 : str
        Label for sequence 2 (used as column header prefix).
    filename : str
        Output file path.
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

    
