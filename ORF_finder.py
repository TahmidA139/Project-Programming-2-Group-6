#!/usr/bin/env python3
"""
ORF_finder.py
Purpose:
    This module analyzes detected ORFs to identify repeated sequences and similarity scores.
Role in Project:
    It performs downstream ORF comparison after ORFs are identified.
Input:
    DNA sequence strings.
Output:
    Collection of ORFs and their metadata.
"""

def find_orfs(dna sequence):
    """
    Objective:
        Identify ORFs with identical sequences.
    Input:
        Cleaned DNA sequence from input_validate.py
    Output:
        A list of dictionaries where each dictionary contains:
            - frame: reading frame numer (0,1,2)
            - start: starting index of the ORF
            - end: ending index of the ORF
            
    High-Level Steps:
        - Define start and stop codons
        - Iterate through all three reading frames
        - Scan the sequence codon by codon
        - Detect start codon
        - Continue iterating for the whole sequence and find ORFs
        - Return all detected ORFs list
        
    """
    pass

def orfs_metadata(dna_sequence):
    """
    Objective:
         Identify the genomic positions and associated metadata for each detected ORF in the DNA sequence.
    Input:
        orfs (list): Detected ORFs.
    Output:
        similarity_scores (dict): Pairwise similarity values.
    High-Level Steps:
        - Compare ORF sequences to compute similarity
    """
    pass
