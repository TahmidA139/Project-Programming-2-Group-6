#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_finder.py

Purpose:
    Public API for ORF detection across all six reading frames.
    Low-level scanning helpers live in frame_scanner.py.

Public API
----------
    find_orfs(dna_sequence, start_codons, min_length, ignore_nested)
        → (nested_dict, flat_list)
    find_nested(flat_list)
        → list of nested ORF records

    CSV_FIELDNAMES : list of str
        Column order for CSV export (used by output_writer.py).
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
    "orf_id", "strand", "start_codon",
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
# Step 2 — nesting filter
# ---------------------------------------------------------------------------

def _apply_nesting(
    all_orfs:      List[Dict[str, Any]],
    ignore_nested: bool,
) -> List[Dict[str, Any]]:
    """Annotate is_nested and optionally remove nested ORFs."""
    all_orfs = _mark_nested(all_orfs)
    if ignore_nested:
        all_orfs = [o for o in all_orfs if not o["is_nested"]]
    return all_orfs


# ---------------------------------------------------------------------------
# Step 3 — build the nested output dictionary
# ---------------------------------------------------------------------------

def _make_nested_dict(
    active_noncanonical: List[str],
) -> Dict[str, Any]:
    """Return an empty nested output dictionary with the correct structure."""
    return {
        "canonical":    {},
        "noncanonical": {sc: {} for sc in active_noncanonical},
    }


def _label_and_insert(
    orf:         Dict[str, Any],
    nested_dict: Dict[str, Any],
    counts:      Dict[str, int],
    active_nc:   List[str],
) -> str:
    """
    Assign a label to one ORF, insert it into nested_dict, and return the label.

    counts keys: 'canonical', '<CODON>'
    """
    sc = orf["start_codon"]

    if sc == CANONICAL_START:
        counts["canonical"] = counts.get("canonical", 0) + 1
        label = f"ORF{counts['canonical']}"
        nested_dict["canonical"][label] = orf

    elif sc in active_nc:
        counts[sc] = counts.get(sc, 0) + 1
        label = f"{sc}_ORF{counts[sc]}"
        nested_dict["noncanonical"][sc][label] = orf

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
        label             = _label_and_insert(orf, nested_dict, counts, active_nc)
        flat_record       = dict(orf)
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
    Find all complete ORFs in all six reading frames of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string (A, T, G, C).
    start_codons : list of str, optional
        Start codons to search for. Default: ['ATG'].
    min_length : int, optional
        Minimum ORF length in nucleotides. Default: 30.
    ignore_nested : bool, optional
        If True, overlapping shorter ORFs are removed, keeping only the
        longest non-overlapping ORF from each overlapping group per strand.
        Matches NCBI ORF Finder's 'Nested ORFs removed' behaviour.
        Default: False.

    Returns
    -------
    nested_dict : dict
        Hierarchical dict organised by canonical/noncanonical.
    flat_list : list of dict
        Flat list of all ORF records with orf_id, strand, is_nested fields.
    """
    dna_sequence = dna_sequence.upper().strip()

    all_orfs = _scan_all_frames(dna_sequence, start_codons, min_length)
    all_orfs = _apply_nesting(all_orfs, ignore_nested)

    return _build_outputs(all_orfs, start_codons)


def find_nested(flat_list: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Return ORFs that overlap with a longer ORF on the same strand.

    Matches NCBI ORF Finder's definition: an ORF is considered nested if
    any longer ORF on the same strand overlaps it (even partially).
    Frame is intentionally NOT considered — this is purely coordinate-based.
    """
    nested_orfs = []
    for i, orf in enumerate(flat_list):
        orf_s = min(orf["start"], orf["end"])
        orf_e = max(orf["start"], orf["end"])
        for j, other in enumerate(flat_list):
            if i == j:
                continue
            if orf["strand"] != other["strand"]:
                continue
            other_s = min(other["start"], other["end"])
            other_e = max(other["start"], other["end"])
            # overlaps AND the other ORF is strictly longer
            if (other_s < orf_e
                    and other_e > orf_s
                    and other["length_nt"] > orf["length_nt"]):
                nested_orfs.append(orf)
                break
    return nested_orfs
