#!/usr/bin/env python3

"""
Input_Validation.py
Purpose:
    This module is responsible for obtaining a DNA sequence from NCBI and preparing it for downstream ORF detection.
Role in Project:
    It serves as the data acquisition and preprocessing step in the ORF analysis pipeline.
Input:
    An NCBI nucleotide accession number provided by the user.
Output:
    A validated and cleaned DNA sequence written to a FASTA file.
"""

def fetch_fasta_from_ncbi(accession):
    """
    Objective:
        Download a DNA sequence in FASTA format from NCBI.
    Input:
        accession (str): NCBI nucleotide accession number.
    Output:
        sequence (str): Raw DNA sequence retrieved from NCBI.
    High-Level Steps:
        - Query NCBI using the accession number
        - Retrieve FASTA record
        - Extract DNA sequence
    """
    pass


def validate_dna_sequence(sequence):
    """
    Objective:
        Verify that the DNA sequence contains only valid nucleotides.
    Input:
        sequence (str): DNA sequence string.
    Output:
        is_valid (bool): Uppercase DNA strings with invalid characters removed.
    High-Level Steps:
        - Check characters against valid bases
        - Flag and remove invalid characters if present
    """
    pass
    

