#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
output_writer.py

Purpose:
    Output helpers for the ORCA pipeline.
    Handles terminal summary printing and writing ORF tables (with extracted
    nucleotide sequences) to a single CSV file, optionally combining results
    from two sequences in comparative mode.

Location:
    src/orf_finder_lib/output_writer.py

Public API
----------
    print_summary(nested, flat_list, nested_count, label)
    write_combined_csv(acc1, flat1, seq1, output_path, acc2, flat2, seq2)
"""

from __future__ import annotations

import csv
import os
from typing import List, Optional

from src.orf_finder_lib.frame_scanner import extract_orf_sequence
from src.orf_finder_lib.orf_finder import CSV_FIELDNAMES

# Fieldnames with the sequence column appended
OUTPUT_FIELDNAMES: List[str] = CSV_FIELDNAMES + ["sequence (5'->3')"]


# ---------------------------------------------------------------------------
# Summary printing
# ---------------------------------------------------------------------------

def print_summary(
    nested:       dict,
    flat_list:    list,
    nested_count: int  = 0,
    label:        str  = "",
) -> None:
    """
    Print a short summary of ORF counts to stdout.

    Parameters
    ----------
    nested : dict
        Nested dict returned by find_orfs().
    flat_list : list
        Flat ORF list returned by find_orfs() (may already be filtered).
    nested_count : int
        Number of nested ORFs detected before any ignore-nested filtering.
        Pass the third return value of find_orfs() here.
    label : str
        Optional label shown in the header (e.g. accession number).
    """
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

    print(f"  Nested ORFs detected        : {nested_count}")
    print("-" * (20 + len(header)))


# ---------------------------------------------------------------------------
# CSV writing
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
    """
    Write one or two ORF tables into a single CSV file.

    Each sequence block begins with a row containing just the accession number,
    followed by the column header row and then one row per ORF. In comparative
    mode the two blocks are separated by two blank rows.

    Parameters
    ----------
    acc1 : str
        Accession number for the first sequence (used as a section header).
    flat1 : list of dict
        Flat ORF list returned by find_orfs() for the first sequence.
    seq1 : str
        Forward-strand DNA sequence for sequence 1 (used to extract ORF seqs).
    output_path : str
        Destination file path for the CSV.
    acc2 : str, optional
        Accession number for the second sequence (comparative mode only).
    flat2 : list of dict, optional
        Flat ORF list for the second sequence (comparative mode only).
    seq2 : str, optional
        Forward-strand DNA sequence for sequence 2 (comparative mode only).
    """
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=OUTPUT_FIELDNAMES, extrasaction="ignore"
        )

        # Sequence 1 block
        _write_sequence_block(fh, writer, acc1, flat1, seq1)

        # Sequence 2 block (comparative mode only)
        if acc2 is not None and flat2 is not None and seq2 is not None:
            fh.write("\n\n")
            _write_sequence_block(fh, writer, acc2, flat2, seq2)
