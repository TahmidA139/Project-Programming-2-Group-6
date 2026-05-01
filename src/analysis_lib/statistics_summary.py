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
    write_gff3                  -- GFF3 annotation file for one sequence.
"""

from __future__ import annotations

import os
import re
from datetime import datetime
from typing import Any, Dict, IO, List, Optional, Tuple

from src.analysis_lib.orf_analysis import codon_usage

_W = 72   # report width


def rule(char: str = "─") -> str:
    """Return a horizontal rule of width ``_W`` built from *char*."""
    return char * _W


def header(title: str, char: str = "═") -> str:
    """Return a centred section header of width ``_W`` using *char* as the border."""
    pad = (_W - len(title) - 2) // 2
    return f"{char * pad} {title} {char * (_W - pad - len(title) - 2)}"


def subheader(title: str) -> str:
    """Return a two-line subheading with a plain underline matching *title*'s length."""
    return f"  {title}\n  {'─' * len(title)}"


def write_run_metadata(
    fh:           IO[str],
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

    with open(filename, "w") as fh:
        fh.write(header("ORF SUMMARY REPORT") + "\n")
        write_run_metadata(fh, [accession], start_codons, min_length)
        write_sequence_section(fh, flat_list, accession)

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
    fh:        IO[str],
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
    fh:    IO[str],
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


# ---------------------------------------------------------------------------
# Comparative report
# ---------------------------------------------------------------------------

def collect_codons(orfs: List[Dict[str, Any]]) -> Dict[str, int]:
    """Return aggregate codon counts across all ORFs in *orfs*."""
    totals: Dict[str, int] = {}
    for orf in orfs:
        for codon, count in codon_usage(orf.get("sequence", "")).items():
            totals[codon] = totals.get(codon, 0) + count
    return totals


def write_codon_section(
    fh:    IO[str],
    flat1: List[Dict[str, Any]],
    flat2: List[Dict[str, Any]],
    acc1:  str,
    acc2:  str,
) -> None:
    """Write a side-by-side aggregate codon-usage table for two ORF sets."""
    codons1    = collect_codons(flat1)
    codons2    = collect_codons(flat2)
    all_codons = sorted(set(codons1) | set(codons2))

    fh.write(header("Aggregate Codon Usage") + "\n\n")
    col = 12
    fh.write(f"  {'Codon':<8} {acc1:>{col}} {acc2:>{col}} {'Delta':>{col}}\n")
    fh.write(f"  {rule('─')}\n")
    for codon in all_codons:
        c1    = codons1.get(codon, 0)
        c2    = codons2.get(codon, 0)
        delta = c1 - c2
        fh.write(f"  {codon:<8} {c1:>{col}} {c2:>{col}} {delta:>{col}}\n")


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


_GFF3_SOURCE = "ORCA"


def _safe_filename(accession: str) -> str:
    """Replace characters that are unsafe in filenames with underscores.

    NCBI accession numbers (e.g. NM_001838.4) are preserved as-is.
    Unusual record IDs from local FASTA headers (e.g. ``lcl|chromosome 1``)
    have spaces, pipes, and slashes replaced so the resulting filename is
    always valid on Linux, macOS, and Windows.

    Parameters
    ----------
    accession : str
        Raw accession string or FASTA record ID.

    Returns
    -------
    str
        Filename-safe version of *accession*.
    """
    return re.sub(r"[^\w.\-]", "_", accession)


def write_gff3(
    flat_list:    List[Dict[str, Any]],
    accession:    str,
    seq_len:      int,
    outdir:       str       = "output",
    start_codons: List[str] = None,
    min_length:   int       = 30,
) -> str:
    """
    Write a GFF3 annotation file for one sequence's ORF set.

    The output filename is derived from *accession* (e.g. ``NM_001838.4.gff3``).
    Coordinates are converted from the pipeline's internal 0-based half-open
    representation to GFF3's 1-based fully-closed convention:
    ``gff_start = orf["start"] + 1``, ``gff_end = orf["end"]`` (unchanged).

    Requires ``calculate_orf_stats()`` to have been called on *flat_list*
    beforehand so that ``gc_content`` and ``protein_length`` are present
    in each ORF record.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Enriched ORF list for one sequence.
    accession : str
        Accession number or record ID; used as both the ``seqid`` column
        value and the output filename stem.
    seq_len : int
        Total sequence length in bp; written into the ``##sequence-region``
        pragma.
    outdir : str, optional
        Output directory.  Created if it does not exist.  Defaults to
        ``'output'``.
    start_codons : List[str], optional
        Start codons used in this run; written as a ``##parameter`` pragma
        for provenance.  Defaults to ``['ATG']``.
    min_length : int, optional
        Minimum ORF length used in this run; written as a ``##parameter``
        pragma for provenance.  Defaults to ``30``.

    Returns
    -------
    str
        Absolute or relative path of the file written.
    """
    if start_codons is None:
        start_codons = ["ATG"]

    os.makedirs(outdir, exist_ok=True)
    filepath = os.path.join(outdir, f"{_safe_filename(accession)}.gff3")

    with open(filepath, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {accession} 1 {seq_len}\n")
        fh.write(f"##source-version ORCA\n")
        fh.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        fh.write(f"##parameter start_codons={','.join(start_codons)}\n")
        fh.write(f"##parameter min_length={min_length}\n")
        for orf in flat_list:
            gff_start = orf["start"] + 1   # 0-based inclusive -> 1-based inclusive
            gff_end   = orf["end"]          # 0-based exclusive == 1-based inclusive
            attrs = (
                f"ID={orf.get('orf_id', '.')};"
                f"start_codon={orf.get('start_codon', '.')};"
                f"length_nt={orf.get('length_nt', '.')};"
                f"gc_content={orf.get('gc_content', 0.0):.2f};"
                f"protein_length={orf.get('protein_length', '.')}"
            )
            fh.write(
                f"{accession}\t{_GFF3_SOURCE}\tORF\t"
                f"{gff_start}\t{gff_end}\t.\t"
                f"{orf.get('strand', '.')}\t0\t{attrs}\n"
            )

    return filepath
