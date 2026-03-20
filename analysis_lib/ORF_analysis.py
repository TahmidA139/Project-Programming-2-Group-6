#!/usr/bin/env python3
"""
(Amanda portion)
ORF_analysis.py
Purpose:
    This module analyzes detected ORFs to identify repeated sequences and similarity scores.
Role in Project:
    It performs downstream ORF comparison after ORFs are identified.
Input:
    A list of ORFs produced by the ORF detection module.
Output:
    Grouped repeated ORFs and similarity scores.
"""

from collections import defaultdict
from itertools import combinations


def _extract_sequence(orf):
    """
    Helper function to get the sequence from an ORF.

    Supports:
      - raw sequence strings
      - dicts with a 'sequence' key

    Raises:
        ValueError: if the ORF format is unsupported
    """
    if isinstance(orf, str):
        return orf
    if isinstance(orf, dict) and "sequence" in orf:
        return orf["sequence"]

    raise ValueError(
        "Each ORF must be either a sequence string or a dict containing a 'sequence' key."
    )


def find_repeated_orfs(orfs):
    """
    Objective:
        Identify ORFs with identical sequences.
    Input:
        orfs (list): Detected ORFs.
    Output:
        repeated_orfs (dict): Groups of matching ORFs.
    High-Level Steps:
        - Compare ORF sequences and group matches
    """
    grouped_orfs = defaultdict(list)

    for index, orf in enumerate(orfs):
        sequence = _extract_sequence(orf)
        grouped_orfs[sequence].append({
            "index": index,
            "orf": orf
        })

    # Keep only sequences that appear more than once
    repeated_orfs = {
        sequence: entries
        for sequence, entries in grouped_orfs.items()
        if len(entries) > 1
    }

    return repeated_orfs


def _pairwise_similarity(seq1, seq2):
    """
    Compute a simple similarity score between two sequences.

    Similarity is defined as:
        matches at aligned positions / length of longer sequence

    Returns:
        float: similarity score between 0.0 and 1.0
    """
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 1.0

    matches = sum(base1 == base2 for base1, base2 in zip(seq1, seq2))
    return matches / max_len


def calculate_similarity_scores(orfs):
    """
    Objective:
        Compute similarity scores between ORF sequences.
    Input:
        orfs (list): Detected ORFs.
    Output:
        similarity_scores (dict): Pairwise similarity values.
    High-Level Steps:
        - Compare ORF sequences to compute similarity
    """
    similarity_scores = {}

    for i, j in combinations(range(len(orfs)), 2):
        seq1 = _extract_sequence(orfs[i])
        seq2 = _extract_sequence(orfs[j])

        similarity_scores[(i, j)] = _pairwise_similarity(seq1, seq2)

    return similarity_scores
orfs = [
    "ATGAAATAG",
    "ATGCCCTAG",
    "ATGAAATAG",
    "ATGAAAtaa".upper()
]

print(find_repeated_orfs(orfs))
print(calculate_similarity_scores(orfs))


