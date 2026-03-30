#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistics_summary.py

Purpose:
    Write ORF analysis results to files. No calculations are performed here —
    all computation is handled by orf_analysis.py, which this module imports.

Output (single-sequence mode):
    output/stats_summary.txt  : human-readable summary report
    output/orf_stats.csv      : per-ORF statistics table

Output (comparative mode, additional files):
    output/comparative_report.txt  : side-by-side human-readable comparison
    output/comparative_stats.csv   : codon-usage delta table
"""

from __future__ import annotations

import csv
import os
from typing import Any, Dict, List

from src.analysis_lib.orf_analysis import (
    calculate_orf_stats,
    _gc_content,
    _protein_length,
    ORF_STATS_FIELDNAMES,
    COMPARATIVE_CSV_FIELDNAMES,
)


# ---------------------------------------------------------------------------
# Single-sequence: human-readable summary report
# ---------------------------------------------------------------------------

def write_stats_to_file(
    flat_list:    List[Dict[str, Any]],
    dna_sequence: str,
    accession:    str,
    outfile:      str = "output/stats_summary.txt",
) -> None:
    """
    Write a human-readable summary report to a text file.

    Includes dataset-level counts, genomic GC content, longest ORF details,
    and a per-ORF table with GC content, protein length, and codon usage.

    Parameters
    ----------
    flat_list : list of dict
        ORF records from ORF_finder.find_orfs().
    dna_sequence : str
        Original forward-strand DNA sequence.
    accession : str
        NCBI accession number (for report header).
    outfile : str
        Path to write the summary report.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True) \
        if os.path.dirname(outfile) else None

    per_orf_stats = calculate_orf_stats(flat_list, dna_sequence)

    total      = len(flat_list)
    complete   = sum(1 for o in flat_list if o["status"].split("|")[0] == "complete")
    incomplete = total - complete
    nested     = sum(1 for o in flat_list if "nested" in o["status"])
    plus       = sum(1 for o in flat_list if o["strand"] == "+")
    minus      = sum(1 for o in flat_list if o["strand"] == "-")
    longest    = max(flat_list, key=lambda o: o["length_nt"]) if flat_list else None
    genomic_gc = _gc_content(dna_sequence)

    with open(outfile, "w") as fh:
        fh.write("=" * 60 + "\n")
        fh.write(f"  ORF Statistics Summary\n")
        fh.write(f"  Accession        : {accession}\n")
        fh.write(f"  Sequence length  : {len(dna_sequence):,} nt\n")
        fh.write(f"  Genomic GC       : {genomic_gc}%\n")
        fh.write("=" * 60 + "\n\n")

        fh.write("--- Dataset Overview ---\n")
        fh.write(f"  Total ORFs found    : {total}\n")
        fh.write(f"  Complete ORFs       : {complete}\n")
        fh.write(f"  Incomplete ORFs     : {incomplete}\n")
        fh.write(f"  Nested ORFs         : {nested}\n")
        fh.write(f"  Forward strand (+)  : {plus}\n")
        fh.write(f"  Reverse strand (-)  : {minus}\n\n")

        if longest:
            prot_len = _protein_length(longest["length_nt"], longest["status"])
            fh.write("--- Longest ORF ---\n")
            fh.write(f"  ID             : {longest['orf_id']}\n")
            fh.write(f"  Strand         : {longest['strand']}\n")
            fh.write(f"  Frame          : {longest['frame']}\n")
            fh.write(f"  Location       : {longest['start']}..{longest['end']}\n")
            fh.write(f"  Span           : {longest['length_nt']} nt\n")
            fh.write(f"  Protein length : {prot_len} aa\n")
            fh.write(f"  Start codon    : {longest['start_codon']}\n")
            fh.write(f"  Status         : {longest['status']}\n\n")

        fh.write("--- Per-ORF Statistics ---\n")
        fh.write(
            f"{'ORF ID':<20} {'Strand':>6} {'Frame':>5} "
            f"{'Start':>7} {'End':>7} {'Length(nt)':>10} "
            f"{'Protein(aa)':>11} {'GC%':>6}\n"
        )
        fh.write("-" * 80 + "\n")

        for rec in per_orf_stats:
            end_str = str(rec["end"]) if rec["end"] is not None else "N/A"
            fh.write(
                f"{rec['orf_id']:<20} {rec['strand']:>6} {rec['frame']:>5} "
                f"{rec['start']:>7} {end_str:>7} {rec['length_nt']:>10} "
                f"{rec['protein_length_aa']:>11} {rec['gc_content_pct']:>6}\n"
            )

        fh.write("\n--- Codon Usage Per ORF ---\n")
        for rec in per_orf_stats:
            fh.write(
                f"\n  {rec['orf_id']} "
                f"({rec['strand']} strand, frame {rec['frame']}):\n"
            )
            items = list(rec["codon_usage"].items())
            for idx in range(0, len(items), 8):
                row = items[idx: idx + 8]
                fh.write("    " + "  ".join(f"{c}:{n}" for c, n in row) + "\n")

    print(f"[INFO] Statistics report written to: {outfile}")


# ---------------------------------------------------------------------------
# Comparative: human-readable side-by-side report
# ---------------------------------------------------------------------------

def write_comparative_report(
    comparison: Dict[str, Any],
    outfile:    str = "output/comparative_report.txt",
) -> None:
    """
    Write a human-readable side-by-side comparative report for two sequences.

    Parameters
    ----------
    comparison : dict
        The dict returned by orf_analysis.compare_orf_sets().
    outfile : str
        Destination path for the text report.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True) \
        if os.path.dirname(outfile) else None

    acc1, acc2 = comparison["accessions"]
    L1,   L2   = comparison["seq_lengths"]
    gc1,  gc2  = comparison["genomic_gc"]
    col        = 16

    def _side(val1: Any, val2: Any) -> str:
        return f"{str(val1):>{col}}  {str(val2):>{col}}"

    with open(outfile, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("  ORCA — Comparative ORF Analysis Report\n")
        fh.write(f"  Sequence 1 : {acc1}\n")
        fh.write(f"  Sequence 2 : {acc2}\n")
        fh.write("=" * 70 + "\n\n")

        # Sequence overview
        fh.write(f"{'Metric':<30} {acc1:>{col}}  {acc2:>{col}}\n")
        fh.write("-" * 66 + "\n")
        fh.write(f"{'Sequence length (nt)':<30} {_side(f'{L1:,}', f'{L2:,}')}\n")
        fh.write(f"{'Genomic GC content (%)':<30} {_side(gc1, gc2)}\n\n")

        # ORF counts
        fh.write("--- ORF Counts ---\n")
        for label, key in [
            ("Total ORFs",         "total_orfs"),
            ("Complete ORFs",      "complete_orfs"),
            ("Incomplete ORFs",    "incomplete_orfs"),
            ("Nested ORFs",        "nested_orfs"),
            ("Forward strand (+)", "plus_strand_orfs"),
            ("Reverse strand (-)", "minus_strand_orfs"),
        ]:
            v1, v2 = comparison[key]
            fh.write(f"  {label:<28} {_side(v1, v2)}\n")
        fh.write("\n")

        # ORF length statistics
        fh.write("--- ORF Length Statistics (nt) ---\n")
        ls1, ls2 = comparison["length_stats"]
        for stat in ("min", "max", "mean"):
            fh.write(f"  {stat.capitalize():<28} {_side(ls1[stat], ls2[stat])}\n")
        fh.write("\n")

        # Mean GC per ORF
        fh.write("--- Mean GC Content per ORF (%) ---\n")
        mg1, mg2 = comparison["mean_gc_per_orf"]
        fh.write(f"  {'Mean ORF GC%':<28} {_side(mg1, mg2)}\n\n")

        # Frame distribution
        fh.write("--- ORFs per Reading Frame ---\n")
        fd1, fd2 = comparison["frame_distribution"]
        for frame in sorted(set(fd1) | set(fd2)):
            fh.write(
                f"  {frame:<28} "
                f"{_side(fd1.get(frame, 0), fd2.get(frame, 0))}\n"
            )
        fh.write("\n")

        # Top-10 codons
        fh.write("--- Top 10 Codons (aggregate across all ORFs) ---\n")
        tc1, tc2 = comparison["top_codons"]
        fh.write(f"\n  {acc1}:\n")
        for codon, count in tc1:
            fh.write(f"    {codon}  {count}\n")
        fh.write(f"\n  {acc2}:\n")
        for codon, count in tc2:
            fh.write(f"    {codon}  {count}\n")
        fh.write("\n")

        # Codon usage delta
        fh.write("--- Codon Usage Differences (seq1 count − seq2 count) ---\n")
        fh.write(f"  {'Codon':<6} {'Seq1':>8} {'Seq2':>8} {'Delta':>8}\n")
        fh.write("  " + "-" * 32 + "\n")
        for codon, (c1, c2, delta) in sorted(
            comparison["codon_usage_delta"].items(),
            key=lambda x: abs(x[1][2]),
            reverse=True,
        ):
            sign = "+" if delta > 0 else ""
            fh.write(f"  {codon:<6} {c1:>8} {c2:>8} {sign}{delta:>7}\n")
        fh.write("\n")

        # Shared / unique start sites
        shared  = comparison["shared_start_sites"]
        unique1 = comparison["unique_to_seq1"]
        unique2 = comparison["unique_to_seq2"]
        fh.write("--- ORF Start-Site Conservation ---\n")
        fh.write(f"  Shared start positions   : {len(shared)}\n")
        fh.write(f"  Unique to {acc1:<12}: {len(unique1)}\n")
        fh.write(f"  Unique to {acc2:<12}: {len(unique2)}\n")
        if shared:
            fh.write(f"\n  Shared positions (first 20): {shared[:20]}\n")

    print(f"[INFO] Comparative report written to: {outfile}")


# ---------------------------------------------------------------------------
# Comparative: codon-usage delta CSV
# ---------------------------------------------------------------------------

def write_comparative_csv(
    comparison: Dict[str, Any],
    outfile:    str = "output/comparative_stats.csv",
) -> None:
    """
    Write a codon-usage delta table to a CSV file.

    Columns: codon, count_<acc1>, count_<acc2>, delta

    Parameters
    ----------
    comparison : dict
        The dict returned by orf_analysis.compare_orf_sets().
    outfile : str
        Destination CSV path.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True) \
        if os.path.dirname(outfile) else None

    acc1, acc2 = comparison["accessions"]
    fieldnames = ["codon", f"count_{acc1}", f"count_{acc2}", "delta"]

    rows = [
        {
            "codon":           codon,
            f"count_{acc1}":   c1,
            f"count_{acc2}":   c2,
            "delta":           delta,
        }
        for codon, (c1, c2, delta) in sorted(
            comparison["codon_usage_delta"].items(),
            key=lambda x: abs(x[1][2]),
            reverse=True,
        )
    ]

    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INFO] Comparative codon-usage CSV written to: {outfile}")
