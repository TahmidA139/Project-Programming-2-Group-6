#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_analysis.py

Amanda Yaworsky's part

Purpose:
    Per-ORF statistics and repeated-ORF detection for the ORCA pipeline.
    All file writing and console reporting functions live in statistics_summary.py.

Exports:
    gc_content          -- Return GC content of a sequence as a percentage.
    protein_length      -- Return the number of translated codons in a sequence.
    codon_usage         -- Return a codon-frequency dict for a sequence.
    calculate_orf_stats -- Enrich each ORF dict with sequence-derived statistics.
    find_repeated_orfs  -- Return ORF sequences that appear more than once.
"""

from __future__ import annotations
from typing import Any, Dict, List
from src.orf_finder_lib.frame_scanner import extract_orf_sequence


def gc_content(sequence: str) -> float:
    """
    Return the GC content of *sequence* as a percentage. Counts uppercase G and C only;
    call ``sequence.upper()`` beforehand if the input may contain lowercase bases.

    When called from ``calculate_orf_stats``, *sequence* is the coding region
    with the stop codon already stripped, so the percentage reflects only the
    translated portion of the ORF.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence.  May be empty.

    Returns
    -------
    float
        GC percentage in ``[0.0, 100.0]``, or ``0.0`` for an empty sequence.
    """
    if not sequence:
        return 0.0
    gc: int = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100


def protein_length(sequence: str) -> int:
    """
    Return the number of translated codons (amino acids) in *sequence*.

    Expects the *coding* region only — i.e. the stop codon must be stripped
    by the caller before passing the sequence here (see ``calculate_orf_stats``).
    Incomplete trailing nucleotides (``len(sequence) % 3 != 0``) are silently
    discarded.

    Parameters
    ----------
    sequence : str
        Coding nucleotide sequence with the stop codon already removed.
        May be empty.

    Returns
    -------
    int
        Number of translated codons, or ``0`` for sequences shorter than 3 nt.
    """
    return len(sequence) // 3


def codon_usage(sequence: str) -> Dict[str, int]:
    """
    Return a codon-frequency dict for *sequence*.
    Iterates non-overlapping triplets from position 0.
    Incomplete trailing codons are ignored.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence.  May be empty.

    Returns
    -------
    Dict[str, int]
        Mapping of three-character codon to its occurrence count.
        Returns an empty dict for sequences shorter than 3 nt.
    """
    counts: Dict[str, int] = {}
    for i in range(0, len(sequence) - 2, 3):
        codon: str = sequence[i : i + 3]
        if all(base in "ACGT" for base in codon):   # skip codons with N
            counts[codon] = counts.get(codon, 0) + 1
    return counts


def calculate_orf_stats(
    flat_list: List[Dict[str, Any]],
    dna_sequence: str,
) -> List[Dict[str, Any]]:
    """
    Enrich each ORF dict in *flat_list* with sequence-derived statistics.

    Adds three keys to every record in-place:
        ``sequence``       -- full nucleotide sequence, 5'->3', strand-corrected,
                             including the stop codon (start codon to stop codon
                             inclusive).  Retained in full so coordinate checks
                             and repeated-ORF comparisons stay consistent.
        ``gc_content``     -- GC percentage of the *coding* region only
                             (stop codon excluded).
        ``protein_length`` -- number of translated codons in the *coding* region
                             (stop codon excluded).

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Flat ORF list returned by ``find_orfs()``.
    dna_sequence : str
        Forward-strand DNA used to find the ORFs.

    Returns
    -------
    List[Dict[str, Any]]
        The same list, mutated in-place and returned.
    """
    for orf in flat_list:
        seq: str = extract_orf_sequence(orf, dna_sequence)
        orf["sequence"] = seq
        coding_seq: str = seq[:-3]          # strip stop codon before computing stats
        orf["gc_content"] = gc_content(coding_seq)
        orf["protein_length"] = protein_length(coding_seq)
    return flat_list


def find_repeated_orfs(flat_list: List[Dict[str, Any]]) -> Dict[str, int]:
    """
    Return ORF nucleotide sequences that appear more than once.

    Requires ``calculate_orf_stats()`` to have been called first so that
    each ORF dict contains a ``"sequence"`` key.  Dicts missing the key
    or with an empty sequence are silently skipped.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Flat ORF list enriched by ``calculate_orf_stats()``.

    Returns
    -------
    Dict[str, int]
        Mapping of repeated sequence to its occurrence count.
    """
    counts: Dict[str, int] = {}
    for orf in flat_list:
        seq: str = orf.get("sequence", "")
        if seq:
            counts[seq] = counts.get(seq, 0) + 1
    return {seq: n for seq, n in counts.items() if n > 1}
