#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orf_analysis.py

Amanda Yaworsky's part 

Purpose:
    Per-ORF statistics and repeated ORF detection for the ORCA pipeline.
    Writing/reporting functions live in statistics_summary.py.

"""

from __future__ import annotations

from typing import Any, Dict, List

from src.orf_finder_lib.frame_scanner import extract_orf_sequence


# ---------------------------------------------------------------------------
# Sequence-level statistics
# ---------------------------------------------------------------------------

def gc_content(sequence: str) -> float:
    """Return the GC content of a nucleotide sequence as a percentage."""
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100


def protein_length(sequence: str) -> int:
    """Return the number of complete codons (amino acids) in a nucleotide sequence."""
    return len(sequence) // 3


def codon_usage(sequence: str) -> Dict[str, int]:
    """
    Return a codon-frequency dict for a nucleotide sequence.
    Incomplete trailing codons are ignored.
    """
    counts: Dict[str, int] = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        counts[codon] = counts.get(codon, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# ORF-level stats
# ---------------------------------------------------------------------------

def calculate_orf_stats(
    flat_list:    List[Dict[str, Any]],
    dna_sequence: str,
) -> List[Dict[str, Any]]:
    """
    Enrich each ORF dict with sequence-derived statistics.

    Adds the following keys to each ORF record in-place:
        - sequence       : nucleotide sequence, 5'->3', correct for strand
        - gc_content     : GC% of that sequence
        - protein_length : number of complete codons

    Parameters
    ----------
    flat_list : list of dict
        Flat ORF list returned by find_orfs().
    dna_sequence : str
        Forward-strand DNA used to find the ORFs.

    Returns
    -------
    The same list, mutated in-place and returned.
    """
    for orf in flat_list:
        seq                  = extract_orf_sequence(orf, dna_sequence)
        orf["sequence"]      = seq
        orf["gc_content"]    = gc_content(seq)
        orf["protein_length"] = protein_length(seq)
    return flat_list


# ---------------------------------------------------------------------------
# Repeated-ORF detection
# ---------------------------------------------------------------------------

def find_repeated_orfs(flat_list: List[Dict[str, Any]]) -> Dict[str, int]:
    """
    Identify ORF nucleotide sequences that appear more than once.

    Requires calculate_orf_stats() to have been called first so that
    each ORF dict contains a 'sequence' key.

    Returns
    -------
    dict mapping repeated sequence -> occurrence count
    """
    counts: Dict[str, int] = {}
    for orf in flat_list:
        seq = orf.get("sequence", "")
        if seq:
            counts[seq] = counts.get(seq, 0) + 1
    return {seq: n for seq, n in counts.items() if n > 1}
