#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ORF_finder.py

Purpose:
    Detect open reading frames (ORFs) in a DNA sequence across all three
    forward reading frames using NumPy vectorization.

Features:
    - NumPy-accelerated codon extraction via stride tricks
    - Detects nested ORFs (every start codon is evaluated independently)
    - Separates ORFs by completeness (complete / incomplete)
    - Separates ORFs by start codon type (ATG canonical; GTG, TTG non-canonical)
    - Returns a structured nested dictionary and a flat list of dicts for CSV export

Output dictionary schema
------------------------
{
    "complete": {
        "canonical": {          # ATG starts
            "ORF1": {frame, start, end, length, start_codon, status},
            ...
        },
        "noncanonical": {
            "GTG": {"GTG_ORF1": {...}, ...},
            "TTG": {"TTG_ORF1": {...}, ...},
        },
    },
    "incomplete": {
        "canonical":    {...},
        "noncanonical": {"GTG": {...}, "TTG": {...}},
    },
}
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

START_CODONS: List[str] = ["ATG", "GTG", "TTG"]
STOP_CODONS:  List[str] = ["TAA", "TAG", "TGA"]
CANONICAL_START: str    = "ATG"
NONCANONICAL_STARTS: List[str] = ["GTG", "TTG"]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _sequence_to_codon_array(dna_sequence: str, frame: int) -> np.ndarray:
    """
    Convert a DNA string into a 1-D NumPy array of 3-character codon strings
    for the given reading frame.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string (A, T, G, C only).
    frame : int
        Reading frame offset (0, 1, or 2).

    Returns
    -------
    np.ndarray of dtype '<U3'
        Array of codons in reading order; shape (n_codons,).
    """
    # Slice the sequence to start at the correct frame, then trim to a
    # multiple of 3 so we can reshape cleanly.
    trimmed = dna_sequence[frame:]
    n_codons = len(trimmed) // 3
    if n_codons == 0:
        return np.array([], dtype="<U3")

    # Build a character array and reshape into (n_codons, 3), then join each
    # row into a 3-character string using np.char operations.
    char_arr = np.array(list(trimmed[: n_codons * 3]), dtype="<U1")
    codon_chars = char_arr.reshape(n_codons, 3)

    # Concatenate the three columns: col0 + col1 + col2
    codons = (
        np.char.add(
            np.char.add(codon_chars[:, 0], codon_chars[:, 1]),
            codon_chars[:, 2],
        )
    )
    return codons  # shape (n_codons,), dtype '<U3'


def _codon_indices_in_sequence(frame: int, codon_index: int) -> int:
    """
    Convert a codon index (within a frame-sliced array) back to the
    nucleotide index in the *original* full sequence.

    Parameters
    ----------
    frame : int
    codon_index : int
        0-based index into the codon array for this frame.

    Returns
    -------
    int
        0-based nucleotide position in the original sequence.
    """
    return frame + codon_index * 3


def _find_stop_codon_index(
    codons: np.ndarray, start_codon_idx: int
) -> int | None:
    """
    Find the index of the first stop codon *after* `start_codon_idx` within
    the codon array `codons`.

    Uses NumPy boolean masking: build a boolean array for each stop codon,
    OR them together, zero out positions at or before the start, then find
    the first True.

    Parameters
    ----------
    codons : np.ndarray
        Full codon array for one reading frame.
    start_codon_idx : int
        Index of the start codon in `codons`.

    Returns
    -------
    int or None
        Index of the first downstream stop codon, or None if none found.
    """
    stop_mask = np.zeros(len(codons), dtype=bool)
    for sc in STOP_CODONS:
        stop_mask |= codons == sc

    # Only look at positions strictly after the start codon
    stop_mask[: start_codon_idx + 1] = False

    candidates = np.nonzero(stop_mask)[0]
    return int(candidates[0]) if candidates.size > 0 else None


# ---------------------------------------------------------------------------
# Core scanning function
# ---------------------------------------------------------------------------

def scan_frame(dna_sequence: str, frame: int) -> List[Dict[str, Any]]:
    """
    Scan a DNA sequence in a single reading frame and return all ORFs,
    including nested ones.

    Every start codon encountered is treated as an independent ORF candidate.
    This means that if a shorter ORF begins inside a longer one (nested ORF),
    both are reported.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string.
    frame : int
        Reading frame offset (0, 1, or 2).

    Returns
    -------
    list of dict
        Each dict has keys:
            frame       : int
            start       : int  (0-based nt index, inclusive)
            end         : int or None  (0-based nt index, exclusive)
            length_nt   : int  (nucleotide length; None → length to sequence end)
            start_codon : str
            status      : 'complete' or 'incomplete'
    """
    codons = _sequence_to_codon_array(dna_sequence, frame)
    if codons.size == 0:
        return []

    # Boolean mask: True where codon is *any* start codon
    start_mask = np.zeros(len(codons), dtype=bool)
    for sc in START_CODONS:
        start_mask |= codons == sc

    start_indices: np.ndarray = np.nonzero(start_mask)[0]  # codon-level indices

    results: List[Dict[str, Any]] = []

    for ci in start_indices:
        start_nt = _codon_indices_in_sequence(frame, int(ci))
        start_codon = str(codons[ci])

        stop_ci = _find_stop_codon_index(codons, int(ci))

        if stop_ci is not None:
            end_nt = _codon_indices_in_sequence(frame, stop_ci) + 3  # exclusive
            length_nt = end_nt - start_nt
            status = "complete"
        else:
            end_nt = None
            length_nt = len(dna_sequence) - start_nt
            status = "incomplete"

        results.append(
            {
                "frame":        frame,
                "start":        start_nt,
                "end":          end_nt,
                "length_nt":    length_nt,
                "start_codon":  start_codon,
                "status":       status,
            }
        )

    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def find_orfs(dna_sequence: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Find all ORFs in all three forward reading frames of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string (A, T, G, C). Non-standard characters are
        tolerated but will not match any start/stop codon.

    Returns
    -------
    nested_dict : dict
        Hierarchical dictionary organised as shown in the module docstring.
    flat_list : list of dict
        Flat list of all ORF records, each augmented with an 'orf_id' field
        suitable for CSV export via csv.DictWriter or pandas.

    Notes
    -----
    Nested ORFs (start codons that fall inside an outer ORF) are fully
    reported — they appear as separate entries in both outputs.
    """
    dna_sequence = dna_sequence.upper().strip()

    all_orfs: List[Dict[str, Any]] = []
    for frame in range(3):
        all_orfs.extend(scan_frame(dna_sequence, frame))

    # ------------------------------------------------------------------
    # Build the nested output dictionary
    # ------------------------------------------------------------------
    # Counters for labelling
    canonical_complete_count:      int = 0
    canonical_incomplete_count:    int = 0
    noncanonical_complete_counts:  Dict[str, int] = {sc: 0 for sc in NONCANONICAL_STARTS}
    noncanonical_incomplete_counts: Dict[str, int] = {sc: 0 for sc in NONCANONICAL_STARTS}

    nested: Dict[str, Any] = {
        "complete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in NONCANONICAL_STARTS},
        },
        "incomplete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in NONCANONICAL_STARTS},
        },
    }

    flat_list: List[Dict[str, Any]] = []

    for orf in all_orfs:
        sc     = orf["start_codon"]
        status = orf["status"]

        if sc == CANONICAL_START:
            if status == "complete":
                canonical_complete_count += 1
                label = f"ORF{canonical_complete_count}"
                nested["complete"]["canonical"][label] = orf
            else:
                canonical_incomplete_count += 1
                label = f"Incomplete_ORF{canonical_incomplete_count}"
                nested["incomplete"]["canonical"][label] = orf
        elif sc in NONCANONICAL_STARTS:
            if status == "complete":
                noncanonical_complete_counts[sc] += 1
                n = noncanonical_complete_counts[sc]
                label = f"{sc}_ORF{n}"
                nested["complete"]["noncanonical"][sc][label] = orf
            else:
                noncanonical_incomplete_counts[sc] += 1
                n = noncanonical_incomplete_counts[sc]
                label = f"{sc}_Incomplete_ORF{n}"
                nested["incomplete"]["noncanonical"][sc][label] = orf

        flat_record = dict(orf)
        flat_record["orf_id"] = label
        flat_list.append(flat_record)

    return nested, flat_list


# ---------------------------------------------------------------------------
# Convenience: CSV field order for export (used by main.py)
# ---------------------------------------------------------------------------

CSV_FIELDNAMES: List[str] = [
    "orf_id",
    "status",
    "start_codon",
    "frame",
    "start",
    "end",
    "length_nt",
]
