#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_finder.py

Nicole Decocker's Part

Purpose:
    High-level ORF detection orchestrator for the ORCA pipeline.
    Coordinates scanning all six reading frames (three forward, three reverse
    complement) of a DNA sequence and assembles the results into a nested
    dictionary (organised by start-codon type) and a flat list suitable for
    CSV output.

    Low-level frame scanning is delegated to frame_scanner.py.

Exports:
    find_orfs        -- Primary public interface; returns (nested_dict, flat_list).

Constants:
    ALL_START_CODONS     -- All recognised start codons: ATG, GTG, TTG.
    STOP_CODONS          -- Universal stop codons: TAA, TAG, TGA.
                           Imported from frame_scanner.py (single source of truth).
    CANONICAL_START      -- The canonical start codon (ATG).
    NONCANONICAL_STARTS  -- Non-canonical start codons (GTG, TTG).
    DEFAULT_START_CODONS -- Default start-codon list used by find_orfs().
    DEFAULT_MIN_LENGTH   -- Default minimum ORF length in nucleotides (30).
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from src.orf_finder_lib.frame_scanner import (
    _reverse_complement,
    scan_frame,
    STOP_CODONS,
)

ALL_START_CODONS:    List[str] = ["ATG", "GTG", "TTG"]
# STOP_CODONS is the single source of truth in frame_scanner.py; imported above.
CANONICAL_START:     str       = "ATG"
NONCANONICAL_STARTS: List[str] = ["GTG", "TTG"]

DEFAULT_START_CODONS: List[str] = ["ATG"]
DEFAULT_MIN_LENGTH:   int       = 30

def _codon_category(codon: str) -> str:
    """
    Return the classification key for a start codon.

    Returns ``"canonical"`` for ATG, or the codon string itself (e.g.
    ``"GTG"``) for a recognised non-canonical codon.  Raises ``ValueError``
    for anything not in ``ALL_START_CODONS`` — this path is unreachable in
    normal use because ``find_orfs`` validates all codons at the entry point.

    Parameters
    ----------
    codon : str
        A three-letter start codon string (uppercase).

    Returns
    -------
    str
        ``"canonical"`` or the non-canonical codon string.

    Raises
    ------
    ValueError
        If *codon* is not in ``ALL_START_CODONS``.
    """
    if codon == CANONICAL_START:
        return "canonical"
    if codon in NONCANONICAL_STARTS:
        return codon
    raise ValueError(
        f"Unrecognised start codon: {codon!r}. Must be one of {ALL_START_CODONS}"
    )


def _scan_all_frames(
    dna_sequence: str,
    start_codons: List[str],
    min_length:   int,
) -> List[Dict[str, Any]]:
    """
    Scan all six reading frames and return a flat list of raw ORF records.

    Iterates over frames 0, 1, and 2 for both the forward strand and the
    reverse complement, collecting every ORF that passes the minimum-length
    filter.  The reverse complement is computed once and reused for all three
    reverse frames.

    Parameters
    ----------
    dna_sequence : str
        Forward-strand DNA sequence (uppercase, ACGT only).
    start_codons : List[str]
        One or more start codons to search for (e.g. ["ATG", "GTG"]).
    min_length : int
        Minimum ORF length in nucleotides; shorter ORFs are discarded.

    Returns
    -------
    List[Dict[str, Any]]
        Flat list of raw ORF record dicts.  Each dict contains at minimum:
        ``strand``, ``frame``, ``start``, ``end``, ``length_nt``,
        ``start_codon``, ``status``.  Coordinates are in the forward-strand
        (0-based) reference frame for both strands.
    """
    seq_len  = len(dna_sequence)
    rev_comp = _reverse_complement(dna_sequence)
    orfs: List[Dict[str, Any]] = []
    for frame in range(3):
        orfs.extend(scan_frame(dna_sequence, frame, start_codons, min_length, "+", seq_len))
        orfs.extend(scan_frame(rev_comp,     frame, start_codons, min_length, "-", seq_len))
    return orfs


def _make_nested_dict(start_codons: List[str]) -> Dict[str, Any]:
    """
    Return an empty nested output dictionary with the correct structure.

    Non-canonical keys are derived directly from *start_codons*, so only
    codons that were actually requested appear in the output — no empty
    buckets for unused codons.

    Parameters
    ----------
    start_codons : List[str]
        The start codons for the current run.

    Returns
    -------
    Dict[str, Any]
        A dict with the shape::

            {
                "canonical":    {},
                "noncanonical": {
                    "GTG": {},   # only present if GTG was requested
                    "TTG": {},   # only present if TTG was requested
                }
            }
    """
    active_nc = [sc for sc in NONCANONICAL_STARTS if sc in start_codons]
    return {
        "canonical":    {},
        "noncanonical": {sc: {} for sc in active_nc},
    }


def _label_and_insert(
    orf:         Dict[str, Any],
    nested_dict: Dict[str, Any],
    counts:      Dict[str, int],
) -> str:
    """
    Assign a human-readable label to one ORF, insert it into nested_dict,
    and return the label.

    Canonical ORFs are numbered sequentially as ``ORF1``, ``ORF2``, …
    Non-canonical ORFs are labelled ``<CODON>_ORF1``, ``<CODON>_ORF2``, …
    (e.g. ``GTG_ORF1``).

    The ``nested_dict["noncanonical"]`` keys are the single source of truth
    for which non-canonical codons are active — no separate ``active_nc``
    list is needed.

    Parameters
    ----------
    orf : Dict[str, Any]
        A single raw ORF record as returned by ``_scan_all_frames``.
    nested_dict : Dict[str, Any]
        The nested output dict that will be mutated in-place.
    counts : Dict[str, int]
        Running counter dict keyed on ``"canonical"`` or the codon string
        (e.g. ``"GTG"``).  Mutated in-place; must be shared across all calls
        within the same ``find_orfs`` invocation.

    Returns
    -------
    str
        The label assigned to this ORF (e.g. ``"ORF3"`` or ``"GTG_ORF1"``).
    """
    sc       = orf["start_codon"]
    category = _codon_category(sc)

    if category == "canonical":
        counts["canonical"] = counts.get("canonical", 0) + 1
        label = f"ORF{counts['canonical']}"
        nested_dict["canonical"][label] = orf
    else:
        counts[sc] = counts.get(sc, 0) + 1
        label = f"{sc}_ORF{counts[sc]}"
        nested_dict["noncanonical"][sc][label] = orf

    return label


def _build_outputs(
    all_orfs:     List[Dict[str, Any]],
    start_codons: List[str],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Build the nested dict and flat list from the processed ORF records.

    Parameters
    ----------
    all_orfs : List[Dict[str, Any]]
        Flat list of raw ORF records from ``_scan_all_frames``.
    start_codons : List[str]
        The start codons used in the current run.

    Returns
    -------
    Tuple[Dict[str, Any], List[Dict[str, Any]]]
        A two-element tuple:

        - **nested_dict** (``Dict[str, Any]``): ORFs grouped by start-codon
          type.  Shape::

              {
                  "canonical":    {"ORF1": {...}, "ORF2": {...}, ...},
                  "noncanonical": {"GTG": {...}, "TTG": {...}},
              }

        - **flat_list** (``List[Dict[str, Any]]``): One dict per ORF, each
          containing all fields from the raw record plus ``"orf_id"``.
    """
    nested_dict = _make_nested_dict(start_codons)
    flat_list:  List[Dict[str, Any]] = []
    counts:     Dict[str, int]       = {}

    for orf in all_orfs:
        label             = _label_and_insert(orf, nested_dict, counts)
        flat_record       = dict(orf)
        flat_record["orf_id"] = label
        flat_list.append(flat_record)

    return nested_dict, flat_list


def find_orfs(
    dna_sequence: str,
    start_codons: List[str] = None,
    min_length:   int       = DEFAULT_MIN_LENGTH,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Find all complete ORFs in all six reading frames of a DNA sequence.

    This is the primary public interface for ORF detection.  It validates
    the requested start codons, normalises the input sequence, delegates
    frame scanning to ``_scan_all_frames``, and assembles the results into
    two complementary data structures.

    Parameters
    ----------
    dna_sequence : str
        Raw DNA sequence string.  May be mixed case and contain leading/
        trailing whitespace — both are normalised internally.  Only ACGT
        characters should appear after cleaning (handled upstream by
        input_validate.py).
    start_codons : List[str], optional
        One or more start codons to search for.  Defaults to ``["ATG"]``.
        Supported values: ``"ATG"``, ``"GTG"``, ``"TTG"``.
    min_length : int, optional
        Minimum ORF length in nucleotides (inclusive).  Defaults to ``30``.
        Must be ≥ 3 (one complete codon).

    Returns
    -------
    Tuple[Dict[str, Any], List[Dict[str, Any]]]
        A two-element tuple:

        - **nested_dict** (``Dict[str, Any]``): ORFs grouped under
          ``"canonical"`` (ATG) and ``"noncanonical"`` (GTG/TTG) keys.
        - **flat_list** (``List[Dict[str, Any]]``): One dict per ORF with
          fields: ``orf_id``, ``strand``, ``start_codon``, ``frame``,
          ``start``, ``end``, ``length_nt``, ``status``.

    Raises
    ------
    ValueError
        If any element of *start_codons* is not in ``ALL_START_CODONS``.

    Examples
    --------
    >>> grouped, flat = find_orfs("ATGAAATAA", min_length=3)
    >>> flat[0]["start_codon"]
    'ATG'
    """
    if start_codons is None:
        start_codons = DEFAULT_START_CODONS

    invalid = [sc for sc in start_codons if sc not in ALL_START_CODONS]
    if invalid:
        raise ValueError(
            f"Unrecognised start codon(s): {invalid}. "
            f"Must be a subset of {ALL_START_CODONS}."
        )

    dna_sequence = dna_sequence.upper().strip()
    all_orfs     = _scan_all_frames(dna_sequence, start_codons, min_length)
    return _build_outputs(all_orfs, start_codons)
