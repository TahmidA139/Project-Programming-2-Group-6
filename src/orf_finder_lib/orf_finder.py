#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_finder.py

Nicole Decocker's Part

Purpose:
    Im loosing my mind trying to get the code to work correclty so when everything is figured out i will add docstrings. 
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from src.orf_finder_lib.frame_scanner import (
    _reverse_complement,
    _mark_nested,
    scan_frame,
)

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

def _apply_nesting(
    all_orfs:      List[Dict[str, Any]],
    ignore_nested: bool,
) -> List[Dict[str, Any]]:
    """Annotate is_nested and optionally remove nested ORFs."""
    all_orfs = _mark_nested(all_orfs)
    if ignore_nested:
        all_orfs = [o for o in all_orfs if not o["is_nested"]]
    return all_orfs

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

def find_orfs(
    dna_sequence:  str,
    start_codons:  List[str] = DEFAULT_START_CODONS,
    min_length:    int       = DEFAULT_MIN_LENGTH,
    ignore_nested: bool      = DEFAULT_IGNORE_NESTED,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Find all complete ORFs in all six reading frames of a DNA sequence.
    """
    dna_sequence = dna_sequence.upper().strip()

    all_orfs = _scan_all_frames(dna_sequence, start_codons, min_length)
    all_orfs = _apply_nesting(all_orfs, ignore_nested)

    return _build_outputs(all_orfs, start_codons)


def find_nested(flat_list: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Return the subset of ORFs that are nested inside another ORF."""
    nested_orfs = []
    for i, orf in enumerate(flat_list):
        for j, other in enumerate(flat_list):
            if i == j:
                continue
            if orf["strand"] != other["strand"]:
                continue
            if orf["frame"] != other["frame"]:
                continue
            if other["start"] < orf["start"] < other["end"]:
                nested_orfs.append(orf)
                break
    return nested_orfs
