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
    - Supports filtering by start codon type, minimum length, and nested status
    - Returns a structured nested dictionary and a flat list of dicts for CSV export

Output dictionary schema
------------------------
{
    "complete": {
        "canonical": {          # ATG starts
            "ORF1": {frame, start, end, length_nt, start_codon, status, is_nested},
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

CSV status values
-----------------
    "complete"             - has start and stop codon, not nested
    "incomplete"           - no stop codon found, not nested
    "complete|nested"      - complete ORF whose start falls inside another ORF
    "incomplete|nested"    - incomplete ORF whose start falls inside another ORF
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALL_START_CODONS:     List[str] = ["ATG", "GTG", "TTG"]
STOP_CODONS:          List[str] = ["TAA", "TAG", "TGA"]
CANONICAL_START:      str       = "ATG"
NONCANONICAL_STARTS:  List[str] = ["GTG", "TTG"]

# Default parameter values (mirrors argparse defaults in main.py)
DEFAULT_START_CODONS:  List[str] = ["ATG"]
DEFAULT_MIN_LENGTH:    int       = 30
DEFAULT_IGNORE_NESTED: bool      = False


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
    trimmed  = dna_sequence[frame:]
    n_codons = len(trimmed) // 3
    if n_codons == 0:
        return np.array([], dtype="<U3")

    char_arr    = np.array(list(trimmed[: n_codons * 3]), dtype="<U1")
    codon_chars = char_arr.reshape(n_codons, 3)

    codons = np.char.add(
        np.char.add(codon_chars[:, 0], codon_chars[:, 1]),
        codon_chars[:, 2],
    )
    return codons  # shape (n_codons,), dtype '<U3'


def _codon_index_to_nt(frame: int, codon_index: int) -> int:
    """
    Convert a codon index (within a frame-sliced array) back to the
    nucleotide index in the original full sequence.
    """
    return frame + codon_index * 3


def _find_stop_codon_index(
    codons: np.ndarray, start_codon_idx: int
) -> Optional[int]:
    """
    Find the index of the first stop codon strictly after `start_codon_idx`
    within `codons` using NumPy boolean masking.

    Returns None if no stop codon is found downstream.
    """
    stop_mask = np.zeros(len(codons), dtype=bool)
    for sc in STOP_CODONS:
        stop_mask |= codons == sc

    stop_mask[: start_codon_idx + 1] = False

    candidates = np.nonzero(stop_mask)[0]
    return int(candidates[0]) if candidates.size > 0 else None


def _mark_nested(all_orfs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Add an 'is_nested' boolean field to every ORF record in-place.

    An ORF is nested if its start position falls strictly inside another
    ORF's start-end range within the same reading frame.

    Parameters
    ----------
    all_orfs : list of dict
        Raw ORF records (must have 'frame', 'start', 'end' keys).

    Returns
    -------
    The same list with 'is_nested' set on every record.
    """
    for i, orf in enumerate(all_orfs):
        is_nested = False
        for j, other in enumerate(all_orfs):
            if i == j:
                continue
            if orf["frame"] != other["frame"]:
                continue
            if other["end"] is None:
                continue
            if other["start"] < orf["start"] < other["end"]:
                is_nested = True
                break
        orf["is_nested"] = is_nested
    return all_orfs


# ---------------------------------------------------------------------------
# Core scanning function
# ---------------------------------------------------------------------------

def scan_frame(
    dna_sequence: str,
    frame:        int,
    start_codons: List[str],
    min_length:   int,
) -> List[Dict[str, Any]]:
    """
    Scan a DNA sequence in a single reading frame and return all ORFs,
    including nested ones, that pass the start-codon and length filters.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string.
    frame : int
        Reading frame offset (0, 1, or 2).
    start_codons : list of str
        Only ORFs beginning with one of these codons are reported.
    min_length : int
        Minimum ORF length in nucleotides. ORFs shorter than this are dropped.

    Returns
    -------
    list of dict
        Keys: frame, start, end, length_nt, start_codon, status
        ('is_nested' is added later by _mark_nested)
    """
    codons = _sequence_to_codon_array(dna_sequence, frame)
    if codons.size == 0:
        return []

    # Build start mask using only the caller-requested start codons
    start_mask = np.zeros(len(codons), dtype=bool)
    for sc in start_codons:
        start_mask |= codons == sc

    start_indices: np.ndarray = np.nonzero(start_mask)[0]

    results: List[Dict[str, Any]] = []

    for ci in start_indices:
        start_nt    = _codon_index_to_nt(frame, int(ci))
        start_codon = str(codons[ci])
        stop_ci     = _find_stop_codon_index(codons, int(ci))

        if stop_ci is not None:
            end_nt    = _codon_index_to_nt(frame, stop_ci) + 3
            length_nt = end_nt - start_nt
            status    = "complete"
        else:
            end_nt    = None
            length_nt = len(dna_sequence) - start_nt
            status    = "incomplete"

        # Apply minimum length filter
        if length_nt < min_length:
            continue

        results.append(
            {
                "frame":       frame,
                "start":       start_nt,
                "end":         end_nt,
                "length_nt":   length_nt,
                "start_codon": start_codon,
                "status":      status,
            }
        )

    return results


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
    Find all ORFs in all three forward reading frames of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        Upper-case DNA string (A, T, G, C). Non-standard characters are
        tolerated but will not match any start/stop codon.
    start_codons : list of str, optional
        Start codons to search for. Default: ["ATG"].
    min_length : int, optional
        Minimum ORF length in nucleotides. Default: 30.
    ignore_nested : bool, optional
        If True, ORFs whose start falls inside another ORF in the same
        frame are removed from the output. Default: False.

    Returns
    -------
    nested_dict : dict
        Hierarchical dictionary organised as shown in the module docstring.
    flat_list : list of dict
        Flat list of all ORF records. The 'status' field encodes both
        completeness and nesting:
            "complete", "incomplete",
            "complete|nested", "incomplete|nested"
    """
    dna_sequence = dna_sequence.upper().strip()

    # ------------------------------------------------------------------
    # 1. Scan all three frames with the user-supplied filters
    # ------------------------------------------------------------------
    all_orfs: List[Dict[str, Any]] = []
    for frame in range(3):
        all_orfs.extend(
            scan_frame(dna_sequence, frame, start_codons, min_length)
        )

    # ------------------------------------------------------------------
    # 2. Annotate every ORF with is_nested
    # ------------------------------------------------------------------
    all_orfs = _mark_nested(all_orfs)

    # ------------------------------------------------------------------
    # 3. Optionally remove nested ORFs
    # ------------------------------------------------------------------
    if ignore_nested:
        all_orfs = [o for o in all_orfs if not o["is_nested"]]

    # ------------------------------------------------------------------
    # 4. Encode nesting into the status field for the CSV
    #    "complete"          -> "complete"
    #    "complete" + nested -> "complete|nested"
    # ------------------------------------------------------------------
    for orf in all_orfs:
        if orf["is_nested"]:
            orf["status"] = f"{orf['status']}|nested"

    # ------------------------------------------------------------------
    # 5. Build the nested output dictionary
    # ------------------------------------------------------------------
    active_noncanonical = [sc for sc in NONCANONICAL_STARTS if sc in start_codons]

    canonical_complete_count:       int            = 0
    canonical_incomplete_count:     int            = 0
    noncanonical_complete_counts:   Dict[str, int] = {sc: 0 for sc in active_noncanonical}
    noncanonical_incomplete_counts: Dict[str, int] = {sc: 0 for sc in active_noncanonical}

    nested_dict: Dict[str, Any] = {
        "complete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in active_noncanonical},
        },
        "incomplete": {
            "canonical":    {},
            "noncanonical": {sc: {} for sc in active_noncanonical},
        },
    }

    flat_list: List[Dict[str, Any]] = []

    for orf in all_orfs:
        sc          = orf["start_codon"]
        # Use only the base status (strip |nested) for routing into the dict
        base_status = orf["status"].split("|")[0]

        if sc == CANONICAL_START:
            if base_status == "complete":
                canonical_complete_count += 1
                label = f"ORF{canonical_complete_count}"
                nested_dict["complete"]["canonical"][label] = orf
            else:
                canonical_incomplete_count += 1
                label = f"Incomplete_ORF{canonical_incomplete_count}"
                nested_dict["incomplete"]["canonical"][label] = orf

        elif sc in active_noncanonical:
            if base_status == "complete":
                noncanonical_complete_counts[sc] += 1
                n     = noncanonical_complete_counts[sc]
                label = f"{sc}_ORF{n}"
                nested_dict["complete"]["noncanonical"][sc][label] = orf
            else:
                noncanonical_incomplete_counts[sc] += 1
                n     = noncanonical_incomplete_counts[sc]
                label = f"{sc}_Incomplete_ORF{n}"
                nested_dict["incomplete"]["noncanonical"][sc][label] = orf

        flat_record           = dict(orf)
        flat_record["orf_id"] = label
        flat_list.append(flat_record)

    return nested_dict, flat_list


# ---------------------------------------------------------------------------
# Convenience: CSV field order for export (used by main.py)
# ---------------------------------------------------------------------------

CSV_FIELDNAMES: List[str] = [
    "orf_id",
    "status",
    "is_nested",
    "start_codon",
    "frame",
    "start",
    "end",
    "length_nt",
]
