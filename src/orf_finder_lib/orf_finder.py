#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_finder.py

Purpose:
    Public API for ORF detection across all six reading frames.
    Low-level scanning helpers live in _frame_scanner.py.

Public API
----------
    find_orfs(dna_sequence, start_codons, min_length, ignore_nested)
        → (nested_dict, flat_list)

    CSV_FIELDNAMES : list of str
        Column order for CSV export (used by main.py).

Constants
---------
    ALL_START_CODONS, STOP_CODONS, CANONICAL_START, NONCANONICAL_STARTS
    DEFAULT_START_CODONS, DEFAULT_MIN_LENGTH, DEFAULT_IGNORE_NESTED
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from src.orf_finder_lib.frame_scanner import (
    _reverse_complement,
    _mark_nested,
    scan_frame,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALL_START_CODONS:    List[str] = ["ATG", "GTG", "TTG"]
STOP_CODONS:         List[str] = ["TAA", "TAG", "TGA"]
CANONICAL_START:     str       = "ATG"
NONCANONICAL_STARTS: List[str] = ["GTG", "TTG"]

DEFAULT_START_CODONS:  List[str] = ["ATG"]
DEFAULT_MIN_LENGTH:    int       = 30
DEFAULT_IGNORE_NESTED: bool      = False

CSV_FIELDNAMES: List[str] = [
    "orf_id", "status", "strand", "start_codon",
    "frame", "start", "end", "length_nt",
]


# ---------------------------------------------------------------------------
# Step 1 — scan all six frames
# ---------------------------------------------------------------------------

def _scan_all_frames(
    dna_sequence: str,
    start_codons: List[str],
    min_length:   int,
) -> List[Dict[str, Any]]:
    """Scan all six reading frames and return a flat list of raw ORF records."""
    seq_len  = len(dna_sequence)
    rev_comp = _reverse_complement(dna_sequence)
    orfs: List[Dict[str, Any]] = []
    for frame in range(3):
        orfs.extend(scan_frame(dna_sequence, frame, start_codons, min_length, "+", seq_len))
        orfs.extend(scan_frame(rev_comp,     frame, start_codons, min_length, "-", seq_len))
    return orfs


# ---------------------------------------------------------------------------
# Step 2 — nesting filter and status encoding
# ---------------------------------------------------------------------------

def _apply_nesting(
    all_orfs:      List[Dict[str, Any]],
    ignore_nested: bool,
) -> List[Dict[str, Any]]:
    """Annotate is_nested, optionally remove nested ORFs, encode into status."""
    all_orfs = _mark_nested(all_orfs)
    if ignore_nested:
        all_orfs = [o for o in all_orfs if not o["is_nested"]]
    for orf in all_orfs:
        if orf["is_nested"]:
            orf["status"] = f"{orf['status']}|nested"
    return all_orfs


# ---------------------------------------------------------------------------
# Step 3 — build the nested output dictionary
# ---------------------------------------------------------------------------

def _make_nested_dict(
    active_noncanonical: List[str],
) -> Dict[str, Any]:
    """Return an empty nested output dictionary with the correct structure."""
    return {
        "complete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in active_noncanonical},
        },
        "incomplete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in active_noncanonical},
        },
    }


def _label_and_insert(
    orf:          Dict[str, Any],
    nested_dict:  Dict[str, Any],
    counts:       Dict[str, int],
    active_nc:    List[str],
) -> str:
    """
    Assign a label to one ORF, insert it into nested_dict, and return the label.

    counts keys: 'canonical_complete', 'canonical_incomplete',
                 '<CODON>_complete', '<CODON>_incomplete'
    """
    sc          = orf["start_codon"]
    base_status = orf["status"].split("|")[0]
    complete    = base_status == "complete"

    if sc == CANONICAL_START:
        key   = f"canonical_{base_status}"
        counts[key] = counts.get(key, 0) + 1
        label = f"ORF{counts[key]}" if complete else f"Incomplete_ORF{counts[key]}"
        nested_dict["complete" if complete else "incomplete"]["canonical"][label] = orf

    elif sc in active_nc:
        key   = f"{sc}_{base_status}"
        counts[key] = counts.get(key, 0) + 1
        n     = counts[key]
        label = f"{sc}_ORF{n}" if complete else f"{sc}_Incomplete_ORF{n}"
        nested_dict["complete" if complete else "incomplete"]["noncanonical"][sc][label] = orf

    else:
        label = "unknown"

    return label


def _build_outputs(
    all_orfs:     List[Dict[str, Any]],
    start_codons: List[str],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """Build the nested dict and flat list from the processed ORF records."""
    active_nc   = [sc for sc in NONCANONICAL_STARTS if sc in start_codons]
    nested_dict = _make_nested_dict(active_nc)
    flat_list:  List[Dict[str, Any]] = []
    counts:     Dict[str, int]       = {}

    for orf in all_orfs:
        label               = _label_and_insert(orf, nested_dict, counts, active_nc)
        flat_record         = dict(orf)
        flat_record["orf_id"] = label
        flat_list.append(flat_record)

    return nested_dict, flat_list


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def find_orfs(
    dna_sequence:  str,
    start_codons:  List[str] = DEFAULT_START_CODONS,
    min_length:    int       = DEFAULT_MIN_LENGTH,
    ignore_nested: bool      = DEFAULT_IGNORE_NESTED,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Find all ORFs in all six reading frames of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string (A, T, G, C).
    start_codons : list of str, optional
        Start codons to search for. Default: ['ATG'].
    min_length : int, optional
        Minimum ORF length in nucleotides. Default: 30.
    ignore_nested : bool, optional
        If True, nested ORFs are removed from the output. Default: False.

    Returns
    -------
    nested_dict : dict
        Hierarchical dict organised by complete/incomplete → canonical/noncanonical.
    flat_list : list of dict
        Flat list of all ORF records with orf_id, strand, status, is_nested fields.
    """
    dna_sequence = dna_sequence.upper().strip()

    all_orfs = _scan_all_frames(dna_sequence, start_codons, min_length)
    all_orfs = _apply_nesting(all_orfs, ignore_nested)

    return _build_outputs(all_orfs, start_codons)
