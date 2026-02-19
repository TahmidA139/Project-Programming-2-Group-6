#!/usr/bin/env python3
"""
ORF_finder.py
Purpose:
    Detect all open reading frames (ORFs) within a DNA sequences across multiple reading frames and their corresponding metadata.
Role in Project:
    Converts raw DNA sequences into candidate protein-coding regions.
Input:
    DNA sequence strings.
Output:
    Collection of ORFs with positions and sequences
"""

def find_orfs(dna sequence):
    """
    Objective:
        Identify ORFs with identical sequences. ##(TO identify the start and stop condons of a sequnce stop codons -UAA, UAG, and UGA) (Start Codon AUG)
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
       Cleaned DNA sequence from input_validate.py.
    Output:
        ORF list with their corresponding metadata
    High-Level Steps:
        - For each ORF in the ORF list, find ORF sequence and location
        - Output ORF sequence and location for each corresponding ORF.
    """
    pass
