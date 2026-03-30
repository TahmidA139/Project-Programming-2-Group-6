#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistics_summary.py

Purpose:
    Generate summary statistics from ORF datasets for reporting.
    Includes GC content per ORF, codon usage, protein length,
    and overall dataset statistics.

    Also provides comparative reporting between two sequences when
    --accession2 is supplied via main.py.

Input:
    flat_list   : list of ORF dicts from ORF_finder.find_orfs()
    dna_sequence: original forward-strand DNA sequence string

Output (single-sequence mode):
    output/stats_summary.txt  : human-readable summary report
    output/orf_stats.csv      : per-ORF statistics table

Output (comparative mode, additional files):
    output/comparative_report.txt  : side-by-side human-readable comparison
    output/comparative_stats.csv   : codon-usage delta table
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
import numpy as np

ORF_STATS_FIELDNAMES: List[str] = [
    "orf_id",
    "strand",
    "frame",
    "start",
    "end",
    "length_nt",
    "protein_length_aa",
    "gc_content_pct",
    "status",
    "start_codon",
]

COMPARATIVE_CSV_FIELDNAMES: List[str] = [
    "codon",
    f"count_seq1",
    f"count_seq2",
    "delta",          # seq1 - seq2
]


# ---------------------------------------------------------------------------
# Internal helpers (unchanged from original)
# ---------------------------------------------------------------------------

def _extract_sequence(orf: Dict[str, Any], dna_sequence: str) -> str:
    """
    Extract the nucleotide sequence of an ORF from the full DNA sequence.

    For '+' strand ORFs the slice is straightforward.
    For '-' strand ORFs the stored coordinates are in forward-strand space
    so we slice and reverse complement.
    """
    start = orf["start"]
    end   = orf["end"] if orf["end"] is not None else len(dna_sequence)

    if orf["strand"] == "+":
        return dna_sequence[start:end]
    else:
        fwd_slice  = dna_sequence[start:end]
        char_arr   = np.array(list(fwd_slice), dtype="<U1")
        _comp      = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        complement = np.vectorize(_comp.get)(char_arr, char_arr)
        return "".join(complement[::-1])


def _gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence using NumPy boolean masking.

    Returns GC content as a percentage (0.0 – 100.0).
    """
    if not sequence:
        return 0.0
    arr = np.array(list(sequence), dtype="<U1")
    gc  = np.isin(arr, ["G", "C"]).sum()
    return round(float(gc / len(sequence) * 100), 2)


def _codon_usage(sequence: str) -> Dict[str, int]:
    """
    Count the occurrence of every codon in a nucleotide sequence using NumPy.

    Codons are extracted in-frame from position 0.
    The stop codon at the end of complete ORFs is included in the count.
    """
    n_codons = len(sequence) // 3
    if n_codons == 0:
        return {}

    char_arr    = np.array(list(sequence[: n_codons * 3]), dtype="<U1")
    codon_chars = char_arr.reshape(n_codons, 3)
    codons      = np.char.add(
        np.char.add(codon_chars[:, 0], codon_chars[:, 1]),
        codon_chars[:, 2],
    )

    unique, counts = np.unique(codons, return_counts=True)
    return dict(sorted(zip(unique.tolist(), counts.tolist())))


def _protein_length(length_nt: int, status: str) -> int:
    """
    Calculate protein length in amino acids from ORF nucleotide length.

    For complete ORFs the stop codon is subtracted (not translated).
    For incomplete ORFs the full length is used as an estimate.
    """
    base_status = status.split("|")[0]
    if base_status == "complete":
        return max(0, (length_nt - 3) // 3)
    else:
        return length_nt // 3


# ---------------------------------------------------------------------------
# Public: per-ORF statistics (unchanged)
# ---------------------------------------------------------------------------

def calculate_orf_stats(
    flat_list:    List[Dict[str, Any]],
    dna_sequence: str,
) -> List[Dict[str, Any]]:
    """
    Compute per-ORF statistics including GC content, codon usage,
    and protein length.
    """
    results: List[Dict[str, Any]] = []

    for orf in flat_list:
        seq      = _extract_sequence(orf, dna_sequence)
        gc       = _gc_content(seq)
        codons   = _codon_usage(seq)
        prot_len = _protein_length(orf["length_nt"], orf["status"])

        record = dict(orf)
        record["gc_content_pct"]    = gc
        record["protein_length_aa"] = prot_len
        record["codon_usage"]       = codons
        results.append(record)

    return results


# ---------------------------------------------------------------------------
# Public: single-sequence summary report (unchanged)
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
    """
    import os
    os.makedirs(os.path.dirname(outfile), exist_ok=True) if os.path.dirname(outfile) else None

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
        fh.write(f"  Accession : {accession}\n")
        fh.write(f"  Sequence length : {len(dna_sequence):,} nt\n")
        fh.write(f"  Genomic GC content : {genomic_gc}%\n")
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
# Public: comparative report  (NEW)
# ---------------------------------------------------------------------------

def write_comparative_report(
    comparison:  Dict[str, Any],
    outfile:     str = "output/comparative_report.txt",
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
    import os
    os.makedirs(os.path.dirname(outfile), exist_ok=True) if os.path.dirname(outfile) else None

    acc1, acc2 = comparison["accessions"]
    L1,   L2   = comparison["seq_lengths"]
    gc1,  gc2  = comparison["genomic_gc"]

    def _side(val1: Any, val2: Any, width: int = 12) -> str:
        """Format two values side by side."""
        return f"{str(val1):>{width}}  {str(val2):>{width}}"

    with open(outfile, "w") as fh:
        # ── Header ──────────────────────────────────────────────────────────
        fh.write("=" * 70 + "\n")
        fh.write("  ORCA — Comparative ORF Analysis Report\n")
        fh.write(f"  Sequence 1 : {acc1}\n")
        fh.write(f"  Sequence 2 : {acc2}\n")
        fh.write("=" * 70 + "\n\n")

        # ── Sequence overview ────────────────────────────────────────────────
        col = 16
        fh.write(f"{'Metric':<30} {acc1:>{col}}  {acc2:>{col}}\n")
        fh.write("-" * 66 + "\n")
        fh.write(f"{'Sequence length (nt)':<30} {_side(f'{L1:,}', f'{L2:,}', col)}\n")
        fh.write(f"{'Genomic GC content (%)':<30} {_side(gc1, gc2, col)}\n")
        fh.write("\n")

        # ── ORF counts ───────────────────────────────────────────────────────
        fh.write("--- ORF Counts ---\n")
        rows = [
            ("Total ORFs",         "total_orfs"),
            ("Complete ORFs",      "complete_orfs"),
            ("Incomplete ORFs",    "incomplete_orfs"),
            ("Nested ORFs",        "nested_orfs"),
            ("Forward strand (+)", "plus_strand_orfs"),
            ("Reverse strand (-)", "minus_strand_orfs"),
        ]
        for label, key in rows:
            v1, v2 = comparison[key]
            fh.write(f"  {label:<28} {_side(v1, v2, col)}\n")
        fh.write("\n")

        # ── ORF length statistics ────────────────────────────────────────────
        fh.write("--- ORF Length Statistics (nt) ---\n")
        ls1, ls2 = comparison["length_stats"]
        for stat in ("min", "max", "mean"):
            fh.write(
                f"  {stat.capitalize():<28} "
                f"{_side(ls1[stat], ls2[stat], col)}\n"
            )
        fh.write("\n")

        # ── GC content per ORF ──────────────────────────────────────────────
        fh.write("--- Mean GC Content per ORF (%) ---\n")
        mg1, mg2 = comparison["mean_gc_per_orf"]
        fh.write(f"  {'Mean ORF GC%':<28} {_side(mg1, mg2, col)}\n\n")

        # ── Reading-frame distribution ───────────────────────────────────────
        fh.write("--- ORFs per Reading Frame ---\n")
        fd1, fd2 = comparison["frame_distribution"]
        all_frames = sorted(set(fd1) | set(fd2))
        for frame in all_frames:
            fh.write(
                f"  {frame:<28} "
                f"{_side(fd1.get(frame, 0), fd2.get(frame, 0), col)}\n"
            )
        fh.write("\n")

        # ── Top-10 codons ────────────────────────────────────────────────────
        fh.write("--- Top 10 Codons (aggregate across all ORFs) ---\n")
        tc1, tc2 = comparison["top_codons"]

        fh.write(f"\n  {acc1}:\n")
        for codon, count in tc1:
            fh.write(f"    {codon}  {count}\n")

        fh.write(f"\n  {acc2}:\n")
        for codon, count in tc2:
            fh.write(f"    {codon}  {count}\n")
        fh.write("\n")

        # ── Codon usage differences ──────────────────────────────────────────
        fh.write("--- Codon Usage Differences (seq1 count − seq2 count) ---\n")
        fh.write(
            f"  {'Codon':<6} {'Seq1':>8} {'Seq2':>8} {'Delta':>8}\n"
        )
        fh.write("  " + "-" * 32 + "\n")

        # Sort by absolute delta (largest difference first)
        delta_rows = sorted(
            comparison["codon_usage_delta"].items(),
            key=lambda x: abs(x[1][2]),
            reverse=True,
        )
        for codon, (c1, c2, delta) in delta_rows:
            sign = "+" if delta > 0 else ""
            fh.write(
                f"  {codon:<6} {c1:>8} {c2:>8} {sign}{delta:>7}\n"
            )
        fh.write("\n")

        # ── Shared / unique start sites ─────────────────────────────────────
        shared   = comparison["shared_start_sites"]
        unique1  = comparison["unique_to_seq1"]
        unique2  = comparison["unique_to_seq2"]

        fh.write("--- ORF Start-Site Conservation ---\n")
        fh.write(f"  Shared start positions   : {len(shared)}\n")
        fh.write(f"  Unique to {acc1:<12}: {len(unique1)}\n")
        fh.write(f"  Unique to {acc2:<12}: {len(unique2)}\n")

        if shared:
            fh.write(
                f"\n  Shared positions (first 20 shown): "
                f"{shared[:20]}\n"
            )

    print(f"[INFO] Comparative report written to: {outfile}")


def write_comparative_csv(
    comparison: Dict[str, Any],
    outfile:    str = "output/comparative_stats.csv",
) -> None:
    """
    Write a codon-usage delta table to a CSV file.

    Columns: codon, count_seq1, count_seq2, delta

    Parameters
    ----------
    comparison : dict
        The dict returned by orf_analysis.compare_orf_sets().
    outfile : str
        Destination CSV path.
    """
    import csv
    import os
    os.makedirs(os.path.dirname(outfile), exist_ok=True) if os.path.dirname(outfile) else None

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
