#!/usr/bin/env python3
"""
ORF_analysis.py
Amandas portion 

This module finds repeated ORFs and calculates simple similarity scores.
"""

def find_repeated_orfs(orfs):
    """
    Find ORFs that appear more than once.

    Input:
        orfs (list): List of ORF sequences
    Output:
        dict: Repeated ORFs and how many times they appear
    """
    counts = {}

    for orf in orfs:
        if orf in counts:
            counts[orf] += 1
        else:
            counts[orf] = 1

    repeated_orfs = {}

    for orf in counts:
        if counts[orf] > 1:
            repeated_orfs[orf] = counts[orf]

    return repeated_orfs


def calculate_similarity_scores(orfs):
    """
    Compare each ORF to every other ORF.

    Input:
        orfs (list): List of ORF sequences
    Output:
        dict: Similarity scores for each pair
    """
    similarity_scores = {}

    for i in range(len(orfs)):
        for j in range(i + 1, len(orfs)):
            seq1 = orfs[i]
            seq2 = orfs[j]

            matches = 0
            length = min(len(seq1), len(seq2))

            for k in range(length):
                if seq1[k] == seq2[k]:
                    matches += 1

            score = matches / length
            similarity_scores[(i, j)] = score

    return similarity_scores
orfs = ["ATGAAA", "ATGAAA", "ATGCCA", "ATGAAA"]

print(find_repeated_orfs(orfs))
print(calculate_similarity_scores(orfs))


