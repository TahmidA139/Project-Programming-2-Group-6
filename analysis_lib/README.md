ORF Analysis Module

This Python module provides basic tools for analyzing Open Reading Frames (ORFs). It focuses on identifying repeated ORFs and calculating simple similarity scores between ORF sequences.

Features
Find Repeated ORFs
Identifies ORF sequences that appear more than once in a list.
Returns a dictionary with ORFs and their counts.
Calculate Similarity Scores
Compares each ORF to every other ORF.
Computes a similarity score based on matching characters at the same positions.
Uses the length of the shorter sequence for comparison.
Functions
find_repeated_orfs(orfs)

Input:

orfs (list): A list of ORF sequences (strings)

Output:

dict: ORFs that appear more than once, with their counts
calculate_similarity_scores(orfs)

Input:

orfs (list): A list of ORF sequences (strings)

Output:

dict: Pairwise similarity scores using index pairs as keys

How it works:

Compares each pair of sequences
Counts matching positions
Divides by the length of the shorter sequence

Notes
Similarity is position-based and does not perform sequence alignment.
Scores range from 0.0 (no matches) to 1.0 (identical sequences).
The function assumes sequences are non-empty; empty sequences may cause errors.
Requirements
Python 3.x
Usage

Run the script directly:

python ORF_analysis.py
