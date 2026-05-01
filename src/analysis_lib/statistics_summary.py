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
from typing import Any, Dict, List, Optional, Tuple

from src.analysis_lib.orf_analysis import codon_usage
from src.orf_finder_lib.frame_scanner import extract_orf_sequence
from src.orf_finder_lib.orf_finder import CSV_FIELDNAMES

OUTPUT_FIELDNAMES: List[str] = CSV_FIELDNAMES + ["sequence (5'->3')"]

_W = 72   # report width


def rule(char: str = "─") -> str:
    return char * _W


def header(title: str, char: str = "═") -> str:
    pad = (_W - len(title) - 2) // 2
    return f"{char * pad} {title} {char * (_W - pad - len(title) - 2)}"


def subheader(title: str) -> str:
    return f"  {title}\n  {'─' * len(title)}"


def write_run_metadata(
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
    fh.write(subheader("Run Parameters") + "\n")
    fh.write(f"  Generated    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    for i, acc in enumerate(accessions, 1):
        fh.write(f"  Sequence {i}   : {acc}\n")
    fh.write(f"  Start codons : {', '.join(start_codons)}\n")
    fh.write(f"  Min length   : {min_length} nt\n\n")



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
    print(f"\n{header(title)}")
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
    print(rule())


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
        fh.write(header("ORF SUMMARY REPORT") + "\n")
        write_run_metadata(fh, [accession], start_codons, min_length)

        fh.write(subheader("Dataset Statistics") + "\n")
        fh.write(f"  Total ORFs         : {total_orfs}\n")
        fh.write(f"  Average GC Content : {avg_gc:.2f}%\n\n")

        if longest:
            fh.write(subheader("Longest ORF") + "\n")
            fh.write(f"  ID         : {longest.get('orf_id', 'N/A')}\n")
            fh.write(f"  Length     : {len(longest['sequence'])} nt\n")
            fh.write(f"  Strand     : {longest.get('strand', '?')}\n")
            fh.write(f"  Frame      : {longest.get('frame', '?')}\n")
            fh.write(f"  GC Content : {longest['gc_content']:.2f}%\n\n")

        fh.write(subheader("Per-ORF Statistics") + "\n")
        fh.write(f"  {'#':<5}{'ID':<14}{'Length':>8}{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n")
        fh.write(f"  {rule('─')}\n")
        for i, orf in enumerate(flat_list):
            fh.write(
                f"  {i:<5}{orf.get('orf_id', ''):<14}"
                f"{len(orf.get('sequence', ''))!s:>8}"
                f"{orf['gc_content']:>8.2f}"
                f"{orf['protein_length']:>10}"
                f"{orf.get('strand', '?'):>8}"
                f"{orf.get('frame', '?')!s:>7}\n"
            )

        fh.write(f"\n{subheader('Aggregate Codon Usage')}\n")
        total_codons: Dict[str, int] = {}
        for orf in flat_list:
            for codon, count in codon_usage(orf.get("sequence", "")).items():
                total_codons[codon] = total_codons.get(codon, 0) + count
        for codon, count in sorted(total_codons.items()):
            fh.write(f"  {codon} : {count}\n")



def avg_gc(orfs: List[Dict[str, Any]]) -> float:
    """Return average GC content across a list of ORFs."""
    return sum(o["gc_content"] for o in orfs) / len(orfs) if orfs else 0.0


def strand_counts(orfs: List[Dict[str, Any]]) -> Tuple[int, int]:
    """Return (plus_count, minus_count) for a list of ORFs."""
    plus  = sum(1 for o in orfs if o.get("strand") == "+")
    minus = sum(1 for o in orfs if o.get("strand") == "-")
    return plus, minus


def write_sequence_section(
    fh,
    flat_list: List[Dict[str, Any]],
    label:     str,
) -> None:
    """Write dataset summary and per-ORF table for one sequence."""
    total_orfs = len(flat_list)
    avg_gc_val = (
        sum(o["gc_content"] for o in flat_list) / total_orfs
        if total_orfs else 0.0
    )
    longest: Optional[Dict[str, Any]] = (
        max(flat_list, key=lambda x: len(x.get("sequence", "")))
        if flat_list else None
    )

    fh.write(header(label) + "\n\n")

    fh.write(subheader("Dataset Statistics") + "\n")
    fh.write(f"  Total ORFs         : {total_orfs}\n")
    fh.write(f"  Average GC Content : {avg_gc_val:.2f}%\n\n")

    if longest:
        fh.write(subheader("Longest ORF") + "\n")
        fh.write(f"  ID         : {longest.get('orf_id', 'N/A')}\n")
        fh.write(f"  Length     : {len(longest['sequence'])} nt\n")
        fh.write(f"  Strand     : {longest.get('strand', '?')}\n")
        fh.write(f"  Frame      : {longest.get('frame', '?')}\n")
        fh.write(f"  GC Content : {longest['gc_content']:.2f}%\n\n")

    fh.write(subheader("Per-ORF Statistics") + "\n")
    fh.write(f"  {'#':<5}{'ID':<14}{'Length':>8}{'GC%':>8}{'Prot_len':>10}{'Strand':>8}{'Frame':>7}\n")
    fh.write(f"  {rule('─')}\n")
    for i, orf in enumerate(flat_list):
        fh.write(
            f"  {i:<5}{orf.get('orf_id', ''):<14}"
            f"{len(orf.get('sequence', ''))!s:>8}"
            f"{orf['gc_content']:>8.2f}"
            f"{orf['protein_length']:>10}"
            f"{orf.get('strand', '?'):>8}"
            f"{orf.get('frame', '?')!s:>7}\n"
        )


def write_comparative_summary(
    fh,
    flat1: List[Dict[str, Any]],
    flat2: List[Dict[str, Any]],
    acc1:  str,
    acc2:  str,
) -> None:
    """Write the shared/unique ORF comparison block."""
    plus1,  minus1 = strand_counts(flat1)
    plus2,  minus2 = strand_counts(flat2)

    seqs1  = {o.get("sequence", "") for o in flat1}
    seqs2  = {o.get("sequence", "") for o in flat2}
    shared = seqs1 & seqs2

    fh.write(header("Comparative Summary") + "\n\n")
    col = 24
    fh.write(f"  {'Metric':<28} {acc1:>{col}} {acc2:>{col}}\n")
    fh.write(f"  {rule('─')}\n")
    fh.write(f"  {'Total ORFs':<28} {len(flat1):>{col}} {len(flat2):>{col}}\n")
    fh.write(f"  {'Forward strand (+)':<28} {plus1:>{col}} {plus2:>{col}}\n")
    fh.write(f"  {'Reverse strand (-)':<28} {minus1:>{col}} {minus2:>{col}}\n")
    fh.write(f"\n  {'Shared ORF sequences':<28} {len(shared):>{col}}\n")
    fh.write(f"  {f'Unique to {acc1}':<28} {len(seqs1 - seqs2):>{col}}\n")
    fh.write(f"  {f'Unique to {acc2}':<28} {len(seqs2 - seqs1):>{col}}\n")


def collect_codons(orfs: List[Dict[str, Any]]) -> Dict[str, int]:
    """Return aggregate codon counts across all ORF sequences in *orfs*.

    Requires ``calculate_orf_stats()`` to have been called first so that
    each ORF dict contains a ``"sequence"`` key.  Only counts codons from
    within detected ORFs, not the full background sequence.
    """
    totals: Dict[str, int] = {}
    for orf in orfs:
        for codon, count in codon_usage(orf.get("sequence", "")).items():
            totals[codon] = totals.get(codon, 0) + count
    return totals


def write_codon_section(
    fh,
    flat1: List[Dict[str, Any]],
    flat2: List[Dict[str, Any]],
    acc1:  str,
    acc2:  str,
) -> None:
    """Write a side-by-side relative-frequency codon table for two ORF sets.

    Frequencies are expressed as a percentage of all valid (ACGT-only) codons
    across each sequence's ORFs.  The Delta column is acc1 minus acc2 in
    percentage points, making it easy to spot codons that are enriched or
    depleted in one sequence relative to the other.
    """
    codons1    = collect_codons(flat1)
    codons2    = collect_codons(flat2)
    total1     = sum(codons1.values()) or 1   # guard against empty ORF sets
    total2     = sum(codons2.values()) or 1
    all_codons = sorted(set(codons1) | set(codons2))

    fh.write(header("Aggregate Codon Usage (% of ORF codons)") + "\n\n")
    col = 12
    fh.write(f"  {'Codon':<8} {acc1:>{col}} {acc2:>{col}} {'Delta (pp)':>{col}}\n")
    fh.write(f"  {rule('─')}\n")
    for codon in all_codons:
        f1    = codons1.get(codon, 0) / total1 * 100
        f2    = codons2.get(codon, 0) / total2 * 100
        delta = f1 - f2
        fh.write(
            f"  {codon:<8} {f1:>{col}.2f}% {f2:>{col}.2f}% {delta:>+{col}.2f}pp\n"
        )


def write_orf_comparison_report(
    flat1:        List[Dict[str, Any]],
    flat2:        List[Dict[str, Any]],
    acc1:         str,
    acc2:         str,
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
    4. Aggregate codon usage (side-by-side counts and delta for every codon).

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
        fh.write(header("ORF COMPARISON REPORT", "═") + "\n")
        write_run_metadata(fh, [acc1, acc2], start_codons, min_length)
        fh.write(rule() + "\n\n")
        write_sequence_section(fh, flat1, f"Sequence 1 — {acc1}")
        fh.write("\n")
        write_sequence_section(fh, flat2, f"Sequence 2 — {acc2}")
        fh.write("\n")
        write_comparative_summary(fh, flat1, flat2, acc1, acc2)
        fh.write("\n")
        write_codon_section(fh, flat1, flat2, acc1, acc2)


def write_sequence_block(
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
        write_sequence_block(fh, writer, acc1, flat1, seq1)

        if acc2 is not None and flat2 is not None and seq2 is not None:
            fh.write("\n\n")
            write_sequence_block(fh, writer, acc2, flat2, seq2)
