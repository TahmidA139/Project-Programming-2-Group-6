#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
frame_scanner.py

Nicole Decocker's part

Purpose:
    Low-level reading-frame scanner for the ORCA pipeline.
    Converts a DNA sequence into codon arrays (via NumPy) and locates all
    ORFs within a single reading frame that pass the minimum-length filter.

    High-level orchestration (scanning all six frames, labelling ORFs,
    building output structures) lives in orf_finder.py.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

STOP_CODONS: List[str] = ["TAA", "TAG", "TGA"]
_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def _reverse_complement(dna_sequence: str) -> str:
    """Return the reverse complement of a DNA sequence using NumPy."""
    char_arr   = np.array(list(dna_sequence), dtype="<U1")
    complement = np.vectorize(_COMPLEMENT.get)(char_arr, char_arr)
    return "".join(complement[::-1])


def _sequence_to_codon_array(dna_sequence: str, frame: int) -> np.ndarray:
    """Convert a DNA string into a 1-D array of 3-character codon strings."""
    trimmed  = dna_sequence[frame:]
    n_codons = len(trimmed) // 3
    if n_codons == 0:
        return np.array([], dtype="<U3")
    char_arr    = np.array(list(trimmed[: n_codons * 3]), dtype="<U1")
    codon_chars = char_arr.reshape(n_codons, 3)
    return np.char.add(
        np.char.add(codon_chars[:, 0], codon_chars[:, 1]),
        codon_chars[:, 2],
    )


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def _codon_index_to_nt(frame: int, codon_index: int) -> int:
    """Convert a codon index within a frame-sliced array to a nucleotide index."""
    return frame + codon_index * 3


def _rc_coords_to_forward(rc_start: int, rc_end: int, seq_len: int) -> Tuple[int, int]:
    """Convert reverse-complement start/end coordinates to forward-strand positions."""
    return seq_len - rc_end, seq_len - rc_start


def _resolve_coords(
    strand: str, rc_start: int, rc_end: int, seq_len: int
) -> Tuple[int, int]:
    """Return forward-strand (start, end) from raw scan coordinates."""
    if strand == "-":
        return _rc_coords_to_forward(rc_start, rc_end, seq_len)
    return rc_start, rc_end


def _find_stop_codon_index(
    codons: np.ndarray, start_codon_idx: int
) -> Optional[int]:
    """Return the index of the first stop codon strictly after start_codon_idx."""
    stop_mask = np.zeros(len(codons), dtype=bool)
    for sc in STOP_CODONS:
        stop_mask |= codons == sc
    stop_mask[: start_codon_idx + 1] = False
    candidates = np.nonzero(stop_mask)[0]
    return int(candidates[0]) if candidates.size > 0 else None

def _process_start_codon(
    ci:         int,
    codons:     np.ndarray,
    frame:      int,
    strand:     str,
    seq_len:    int,
    min_length: int,
) -> Tuple[Optional[Dict[str, Any]], Optional[int]]:
    """
    Build one ORF record for a single start codon index.

    Returns ``(None, None)`` if no downstream stop codon exists or the
    resulting ORF is shorter than *min_length*.  Otherwise returns both
    the ORF record dict and the stop codon index so the caller can use
    the stop index as a deduplication key without recomputing it.

    Parameters
    ----------
    ci : int
        Codon index of the start codon within *codons*.
    codons : np.ndarray
        Full codon array for this frame.
    frame : int
        Reading frame offset (0, 1, or 2).
    strand : str
        ``"+"`` or ``"-"``.
    seq_len : int
        Total length of the forward-strand sequence.
    min_length : int
        Minimum ORF length in nucleotides.

    Returns
    -------
    Tuple[Dict[str, Any], int]
        ``(orf_record, stop_ci)`` on success.
    Tuple[None, None]
        When no valid ORF can be built from this start codon.
    """
    rc_start = _codon_index_to_nt(frame, ci)
    stop_ci  = _find_stop_codon_index(codons, ci)

    if stop_ci is None:
        return None, None

    rc_end    = _codon_index_to_nt(frame, stop_ci) + 3
    length_nt = rc_end - rc_start

    if length_nt < min_length:
        return None, None

    start, end = _resolve_coords(strand, rc_start, rc_end, seq_len)
    record = {
        "strand":      strand,
        "frame":       frame,
        "start":       start,
        "end":         end,
        "length_nt":   length_nt,
        "start_codon": str(codons[ci]),
        "status":      "complete",
    }
    return record, stop_ci


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def scan_frame(
    dna_sequence: str,
    frame:        int,
    start_codons: List[str],
    min_length:   int,
    strand:       str,
    seq_len:      int,
) -> List[Dict[str, Any]]:
    """
    Scan one reading frame and return all complete ORFs passing filters.

    For each stop codon position, only the longest ORF (i.e. the one with
    the earliest start codon) is kept, matching standard longest-ORF
    convention.  The stop codon index returned by ``_process_start_codon``
    is reused directly as the deduplication key — ``_find_stop_codon_index``
    is therefore called exactly once per start codon candidate.

    Parameters
    ----------
    dna_sequence : str
        The sequence to scan (forward or reverse-complement strand).
    frame : int
        Reading-frame offset: 0, 1, or 2.
    start_codons : List[str]
        Start codons to search for.  Already validated by ``find_orfs``.
    min_length : int
        Minimum ORF length in nucleotides.
    strand : str
        ``"+"`` for the forward strand, ``"-"`` for the reverse complement.
    seq_len : int
        Length of the original forward-strand sequence; used for coordinate
        conversion on the ``"-"`` strand.

    Returns
    -------
    List[Dict[str, Any]]
        One record per unique stop-codon position (longest ORF wins).
    """
    codons = _sequence_to_codon_array(dna_sequence, frame)
    if codons.size == 0:
        return []

    start_mask = np.zeros(len(codons), dtype=bool)
    for sc in start_codons:
        start_mask |= codons == sc

    # Iterating start codons earliest-first means the first record stored for
    # each stop position is always the longest — later (shorter) ORFs for the
    # same stop codon are simply skipped.
    best_per_stop: Dict[int, Dict[str, Any]] = {}
    for ci in np.nonzero(start_mask)[0]:
        record, stop_ci = _process_start_codon(
            int(ci), codons, frame, strand, seq_len, min_length
        )
        if record is None:
            continue
        stop_key = _codon_index_to_nt(frame, stop_ci)
        if stop_key not in best_per_stop:
            best_per_stop[stop_key] = record   # first (longest) wins

    return list(best_per_stop.values())


def extract_orf_sequence(orf: dict, forward_seq: str) -> str:
    """
    Extract the ORF nucleotide sequence in 5'→3' direction.

    For ``'+'`` strand ORFs, slices directly from *forward_seq*.
    For ``'-'`` strand ORFs, returns the reverse complement of the relevant
    slice so the result always reads 5'→3'.

    Parameters
    ----------
    orf : dict
        ORF record containing ``"start"``, ``"end"``, and ``"strand"`` keys.
    forward_seq : str
        The original forward-strand DNA sequence.

    Returns
    -------
    str
        Nucleotide sequence of the ORF, 5'→3'.
    """
    start  = orf["start"]
    end    = orf["end"]
    strand = orf["strand"]

    if strand == "+":
        return forward_seq[start:end]
    return _reverse_complement(forward_seq[start:end])
