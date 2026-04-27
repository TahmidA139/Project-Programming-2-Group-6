#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_analysis.py

Amanda Yaworsky's part

Purpose:
    Per-ORF statistics and repeated-ORF detection for the ORCA pipeline.
    All file writing and console reporting functions live in statistics_summary.py.
"""

from __future__ import annotations
from typing import Any, Dict, List
from src.orf_finder_lib.frame_scanner import extract_orf_sequence

def gc_content(sequence: str) -> float:
    """
    Return the GC content of *sequence* as a percentage. Counts uppercase G and C only;
    call ``sequence.upper()`` beforehand if the input may contain lowercase bases.

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
    Return the number of complete codons encoded by *sequence*.
    Incomplete trailing nucleotides (``len(sequence) % 3 != 0``) are
    silently discarded.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence.  May be empty.

    Returns
    -------
    int
        Number of complete codons, or ``0`` for sequences shorter than 3 nt.
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
        counts[codon] = counts.get(codon, 0) + 1
    return counts

def calculate_orf_stats(
    flat_list: List[Dict[str, Any]],
    dna_sequence: str,
) -> List[Dict[str, Any]]:
    """
    Enrich each ORF dict in *flat_list* with sequence-derived statistics.

    Adds three keys to every record in-place:
        ``sequence``       -- nucleotide sequence, 5'->3', strand-corrected.
        ``gc_content``     -- GC percentage of that sequence.
        ``protein_length`` -- number of complete codons.

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
        orf["gc_content"] = gc_content(seq)
        orf["protein_length"] = protein_length(seq)
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
