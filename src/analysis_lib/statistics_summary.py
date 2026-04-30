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
from datetime import datetime
from typing import Any, Dict, List, Optional

from src.analysis_lib.orf_analysis import (
    codon_usage, gc_content, protein_length,
    global_alignment_stats, local_alignment_stats,
)
from src.orf_finder_lib.frame_scanner import extract_orf_sequence
from src.orf_finder_lib.orf_finder import CSV_FIELDNAMES

OUTPUT_FIELDNAMES: List[str] = CSV_FIELDNAMES + ["sequence (5'->3')"]

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

_W = 72   # report width

def _rule(char: str = "─") -> str:
    return char * _W

def _header(title: str, char: str = "═") -> str:
    pad = (_W - len(title) - 2) // 2
    return f"{char * pad} {title} {char * (_W - pad - len(title) - 2)}"

def _subheader(title: str) -> str:
    return f"  {title}\n  {'─' * (len(title))}"


# ---------------------------------------------------------------------------
# Run metadata helper
# ---------------------------------------------------------------------------

def _write_run_metadata(
    fh,
    accessions:   List[str],
    start_codons: List[str],
    min_length:   int,
) -> None:
    """Write a run parameters block at the top of a report file.

    Parameters
    ----------
    fh :
        Open file handle to write into.
    accessions : List[str]
        One or two accession numbers / labels for the run.
    start_codons : List[str]
        Start codons used in this run, e.g. ['ATG', 'GTG'].
    min_length : int
        Minimum ORF length in nucleotides used in this run.
    """
    fh.write(_subheader("Run Parameters") + "\n")
    fh.write(f"  Generated    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    for i, acc in enumerate(accessions, 1):
        fh.write(f"  Sequence {i}   : {acc}\n")
    fh.write(f"  Start codons : {', '.join(start_codons)}\n")
    fh.write(f"  Min length   : {min_length} nt\n\n")


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

    title = f"ORF Summary — {label}" if label else "ORF Summary"
    print(f"\n{_header(title)}")
    print(f"  Total ORFs found  : {total}")
    print(f"  Forward strand (+): {plus_strand}")
    print(f"  Reverse strand (-): {minus_strand}")
    print(f"  Canonical (ATG)   : {n_canonical}")
    if n_noncanonical > 0:
        print(f"  Non-canonical     : {n_noncanonical}")
        for sc in ("GTG", "TTG"):
            n = len(noncanonical.get(sc, {}))
            if n > 0:
                print(f"    {sc}               : {n}")
    print(_rule())


# ---------------------------------------------------------------------------
# Single-sequence summary
# ---------------------------------------------------------------------------

def write_stats_to_file(
    flat_list:    List[Dict[str, Any]],
    filename:     str       = "output/orf_summary.txt",
    accession:    str       = "N/A",
    start_codons: List[str] = None,
    min_length:   int       = 30,
) -> None:
    """Write a human-readable summary report for one sequence's ORF set.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Enriched ORF list (requires calculate_orf_stats() to have been called).
    filename : str
        Output file path.
    accession : str
        Accession number or label for the sequence.
    start_codons : List[str]
        Start codons used in this run.  Defaults to ['ATG'].
    min_length : int
        Minimum ORF length in nucleotides used in this run.
    """
    if start_codons is None:
        start_codons = ["ATG"]

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
        fh.write(_header("ORF SUMMARY REPORT") + "\n")
        _write_run_metadata(fh, [accession], start_codons, min_length)

        fh.write(_subheader("Dataset Statistics") + "\n")
        fh.write(f"  Total ORFs         : {total_orfs}\n")
        fh.write(f"  Average GC Content : {avg_gc:.2f}%\n\n")

        if longest:
            fh.write(_subheader("Longest ORF") + "\n")
            fh.write(f"  ID         : {longest.get('orf_id', 'N/A')}\n")
            fh.write(f"  Length     : {len(longest['sequence'])} nt\n")
            fh.write(f"  Strand     : {longest.get('strand', '?')}\n")
            fh.write(f"  Frame      : {longest.get('frame', '?')}\n")
            fh.write(f"  GC Content : {longest['gc_content']:.2f}%\n\n")

        fh.write(_subheader("Per-ORF Statistics") + "\n")
        fh.write(f"  {'#':<5}{'ID':<14}{'Length':>8}{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n")
        fh.write(f"  {_rule('─')}\n")
        for i, orf in enumerate(flat_list):
            fh.write(
                f"  {i:<5}{orf.get('orf_id', ''):<14}"
                f"{len(orf.get('sequence', ''))!s:>8}"
                f"{orf['gc_content']:>8.2f}"
                f"{orf['protein_length']:>10}"
                f"{orf.get('strand', '?'):>8}"
                f"{orf.get('frame', '?')!s:>7}\n"
            )

        fh.write(f"\n{_subheader('Aggregate Codon Usage')}\n")
        total_codons: Dict[str, int] = {}
        for orf in flat_list:
            for codon, count in codon_usage(orf.get("sequence", "")).items():
                total_codons[codon] = total_codons.get(codon, 0) + count
        for codon, count in sorted(total_codons.items()):
            fh.write(f"  {codon} : {count}\n")


# ---------------------------------------------------------------------------
# Private helpers for comparative report
# ---------------------------------------------------------------------------

def _write_sequence_section(
    fh,
    flat_list: List[Dict[str, Any]],
    label:     str,
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

    fh.write(_header(label) + "\n\n")

    fh.write(_subheader("Dataset Statistics") + "\n")
    fh.write(f"  Total ORFs         : {total_orfs}\n")
    fh.write(f"  Average GC Content : {avg_gc:.2f}%\n\n")

    if longest:
        fh.write(_subheader("Longest ORF") + "\n")
        fh.write(f"  ID         : {longest.get('orf_id', 'N/A')}\n")
        fh.write(f"  Length     : {len(longest['sequence'])} nt\n")
        fh.write(f"  Strand     : {longest.get('strand', '?')}\n")
        fh.write(f"  Frame      : {longest.get('frame', '?')}\n")
        fh.write(f"  GC Content : {longest['gc_content']:.2f}%\n\n")

    fh.write(_subheader("Per-ORF Statistics") + "\n")
    fh.write(f"  {'#':<5}{'ID':<14}{'Length':>8}{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n")
    fh.write(f"  {_rule('─')}\n")
    for i, orf in enumerate(flat_list):
        fh.write(
            f"  {i:<5}{orf.get('orf_id', ''):<14}"
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

    fh.write(_header("Comparative Summary") + "\n\n")
    col = 24
    fh.write(f"  {'Metric':<28} {acc1:>{col}} {acc2:>{col}}\n")
    fh.write(f"  {_rule('─')}\n")
    fh.write(f"  {'Total ORFs':<28} {len(flat1):>{col}} {len(flat2):>{col}}\n")
    fh.write(f"  {'Forward strand (+)':<28} {plus1:>{col}} {plus2:>{col}}\n")
    fh.write(f"  {'Reverse strand (-)':<28} {minus1:>{col}} {minus2:>{col}}\n")
    fh.write(f"\n  {'Shared ORF sequences':<28} {len(shared):>{col}}\n")
    fh.write(f"  {f'Unique to {acc1}':<28} {len(seqs1 - seqs2):>{col}}\n")
    fh.write(f"  {f'Unique to {acc2}':<28} {len(seqs2 - seqs1):>{col}}\n")


def _interpret_identity(pct: float) -> str:
    if pct >= 95:  return "highly conserved (≥95%)"
    if pct >= 70:  return "moderately conserved (70–94%)"
    if pct >= 40:  return "distantly related (40–69%)"
    return "low similarity (<40%) — may be unrelated"


def _write_alignment_section(
    fh,
    seq1: str,
    seq2: str,
    acc1: str,
    acc2: str,
) -> None:
    """Compute global and local alignments and write a dual summary block."""
    g = global_alignment_stats(seq1, seq2)
    l = local_alignment_stats(seq1, seq2)

    fh.write(_header("Sequence Alignment") + "\n\n")

    # Global
    fh.write(_subheader("Global Alignment  (Needleman-Wunsch)") + "\n")
    fh.write("  Scoring: match +1 | mismatch 0 | gap open −2 | gap extend −0.5\n")
    fh.write(f"  Note: global alignment spans the full length of both sequences.\n")
    fh.write(f"  High gap counts are expected when sequence lengths differ greatly.\n\n")
    fh.write(f"  {'Sequence 1 length':<28}: {g['seq1_len']:,} bp  ({acc1})\n")
    fh.write(f"  {'Sequence 2 length':<28}: {g['seq2_len']:,} bp  ({acc2})\n")
    fh.write(f"  {'Alignment length':<28}: {g['alignment_length']:,} bp\n")
    fh.write(f"  {'Matches':<28}: {g['matches']:,}\n")
    fh.write(f"  {'Mismatches':<28}: {g['mismatches']:,}\n")
    fh.write(f"  {'Gaps':<28}: {g['gaps']:,}\n")
    fh.write(f"  {'Identity':<28}: {g['identity_pct']:.2f}%\n")
    fh.write(f"  {'Coverage (vs longer seq)':<28}: {g['coverage_pct']:.2f}%\n")
    fh.write(f"  {'Raw score':<28}: {g['score']:.1f}\n")
    fh.write(f"  {'Interpretation':<28}: {_interpret_identity(g['identity_pct'])}\n\n")

    # Local
    fh.write(_subheader("Local Alignment  (Smith-Waterman)") + "\n")
    fh.write("  Scoring: match +2 | mismatch −1 | gap open −2 | gap extend −0.5\n")
    fh.write(f"  Note: reports only the best-matching subsequence region.\n")
    fh.write(f"  A high identity here with low global identity indicates a conserved domain.\n\n")
    fh.write(f"  {'Best region length':<28}: {l['alignment_length']:,} bp\n")
    fh.write(f"  {'Matches':<28}: {l['matches']:,}\n")
    fh.write(f"  {'Mismatches':<28}: {l['mismatches']:,}\n")
    fh.write(f"  {'Gaps':<28}: {l['gaps']:,}\n")
    fh.write(f"  {'Identity (local region)':<28}: {l['identity_pct']:.2f}%\n")
    fh.write(f"  {'Raw score':<28}: {l['score']:.1f}\n")
    fh.write(f"  {'Interpretation':<28}: {_interpret_identity(l['identity_pct'])}\n")


# ---------------------------------------------------------------------------
# Comparative report (comparative mode only)
# ---------------------------------------------------------------------------

def write_orf_comparison_report(
    flat1:        List[Dict[str, Any]],
    flat2:        List[Dict[str, Any]],
    acc1:         str,
    acc2:         str,
    seq1:         str,
    seq2:         str,
    filename:     str       = "output/orf_comparison_report.txt",
    start_codons: List[str] = None,
    min_length:   int       = 30,
) -> None:
    """
    Write a single combined report for a two-sequence comparative run.

    Sections
    --------
    1. Run parameters (accessions, start codons, min length, timestamp).
    2. Per-sequence ORF summaries (dataset stats, longest ORF, per-ORF table).
    3. Comparative ORF summary (shared / unique ORF sequences, strand counts).
    4. Pairwise alignment (global Needleman-Wunsch + local Smith-Waterman).

    Parameters
    ----------
    flat1 : List[Dict[str, Any]]
        Enriched ORF list for sequence 1 (requires calculate_orf_stats()).
    flat2 : List[Dict[str, Any]]
        Enriched ORF list for sequence 2.
    acc1 : str
        Accession or label for sequence 1.
    acc2 : str
        Accession or label for sequence 2.
    seq1 : str
        Cleaned DNA sequence for sequence 1; used for alignment.
    seq2 : str
        Cleaned DNA sequence for sequence 2; used for alignment.
    filename : str
        Output file path.
    start_codons : List[str]
        Start codons used in this run.  Defaults to ['ATG'].
    min_length : int
        Minimum ORF length in nucleotides used in this run.
    """
    if start_codons is None:
        start_codons = ["ATG"]

    with open(filename, "w") as fh:
        fh.write(_header("ORF COMPARISON REPORT", "═") + "\n")
        _write_run_metadata(fh, [acc1, acc2], start_codons, min_length)
        fh.write(_rule() + "\n\n")
        _write_sequence_section(fh, flat1, f"Sequence 1 — {acc1}")
        fh.write("\n")
        _write_sequence_section(fh, flat2, f"Sequence 2 — {acc2}")
        fh.write("\n")
        _write_comparative_summary(fh, flat1, flat2, acc1, acc2)
        fh.write("\n")
        _write_alignment_section(fh, seq1, seq2, acc1, acc2)


# ---------------------------------------------------------------------------
# Combined CSV output
# ---------------------------------------------------------------------------

def _write_sequence_block(
    fh,
    writer:       csv.DictWriter,
    accession:    str,
    flat_list:    list,
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
