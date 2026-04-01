#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
frame_scanner.py

Purpose:
     
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

STOP_CODONS: List[str] = ["TAA", "TAG", "TGA"]
_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",}


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


def _codon_index_to_nt(frame: int, codon_index: int) -> int:
    """Convert a codon index within a frame-sliced array to a nucleotide index."""
    return frame + codon_index * 3


def _rc_coords_to_forward(
    rc_start: int, rc_end: int, seq_len: int
) -> Tuple[int, int]:
    """Convert reverse complement start/end coordinates to forward-strand positions."""
    fwd_end   = seq_len - rc_start
    fwd_start = seq_len - rc_end
    return fwd_start, fwd_end


# ---------------------------------------------------------------------------
# Stop codon search
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Nesting annotation
# ---------------------------------------------------------------------------

def _mark_nested(all_orfs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Mark ORFs that overlap with a longer ORF on the same strand.

    An ORF is considered nested/overlapping (and marked is_nested=True) if
    there exists any other ORF on the same strand that:
      1. Overlaps it (even partially), AND
      2. Is strictly longer.

    This matches NCBI ORF Finder's 'Nested ORFs removed' behaviour, which
    removes shorter ORFs that overlap with any longer ORF on the same strand —
    not just ORFs that are fully contained.
    """
    for i, orf in enumerate(all_orfs):
        orf_s = min(orf["start"], orf["end"])
        orf_e = max(orf["start"], orf["end"])
        orf["is_nested"] = any(
            orf["strand"] == other["strand"]
            and other["length_nt"] > orf["length_nt"]
            and min(other["start"], other["end"]) < orf_e   # intervals overlap
            and max(other["start"], other["end"]) > orf_s
            for j, other in enumerate(all_orfs)
            if i != j and other.get("end") is not None and orf.get("end") is not None
        )
    return all_orfs


# ---------------------------------------------------------------------------
# Core frame scanner
# ---------------------------------------------------------------------------

def _resolve_coords(
    strand: str, rc_start: int, rc_end: int, seq_len: int
) -> Tuple[int, int]:
    """Return forward-strand (start, end) from raw rc coordinates."""
    if strand == "-":
        return _rc_coords_to_forward(rc_start, rc_end, seq_len)
    return rc_start, rc_end


def _process_start_codon(
    ci: int, codons: np.ndarray, frame: int,
    strand: str, seq_len: int, min_length: int,
) -> Optional[Dict[str, Any]]:
    """
    Build one ORF record for a single start codon index, or None if no stop
    codon is found or the ORF is too short.
    """
    rc_start = _codon_index_to_nt(frame, ci)
    stop_ci  = _find_stop_codon_index(codons, ci)

    if stop_ci is None:
        return None

    rc_end    = _codon_index_to_nt(frame, stop_ci) + 3
    length_nt = rc_end - rc_start

    if length_nt < min_length:
        return None

    start, end = _resolve_coords(strand, rc_start, rc_end, seq_len)
    return {
        "strand": strand, "frame": frame,
        "start": start, "end": end,
        "length_nt": length_nt, "start_codon": str(codons[ci]),
        "status": "complete",
    }


def scan_frame(
    dna_sequence: str, frame: int, start_codons: List[str],
    min_length: int, strand: str, seq_len: int,
) -> List[Dict[str, Any]]:
    """Scan one reading frame and return all complete ORFs passing filters."""
    codons = _sequence_to_codon_array(dna_sequence, frame)
    if codons.size == 0:
        return []

    start_mask = np.zeros(len(codons), dtype=bool)
    for sc in start_codons:
        start_mask |= codons == sc

    results = []
    for ci in np.nonzero(start_mask)[0]:
        record = _process_start_codon(
            int(ci), codons, frame, strand, seq_len, min_length
        )
        if record is not None:
            results.append(record)
    return results


# ---------------------------------------------------------------------------
# ORF sequence extraction
# ---------------------------------------------------------------------------

def extract_orf_sequence(orf: dict, forward_seq: str) -> str:
    """
    Extract the ORF nucleotide sequence in 5'->3' direction.

    For '+' strand ORFs, slices directly from forward_seq.
    For '-' strand ORFs, takes the reverse complement of the relevant slice
    so the result always reads 5'->3'.
    """
    start  = orf["start"]
    end    = orf["end"]
    strand = orf["strand"]

    if strand == "+":
        return forward_seq[start:end]
    else:
        return _reverse_complement(forward_seq[start:end])
