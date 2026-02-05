#!/usr/bin/env python3

"""
Purpose:
    Downloads a FASTA DNA sequence from NCBI and prepares it for ORF detection.

Input:
    NCBI nucleotide accession number (provided by the user)

Output:
  A cleaned DNA sequence file in FASTA format

High-level steps:
    1. Prompt the user for an NCBI accession number.
    2. Download the FASTA sequence from NCBI.
    4. Validate that the sequence contains only valid DNA bases.
    5. Write the cleaned sequence to an output file.
"""
