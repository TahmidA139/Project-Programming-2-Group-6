#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Nicole's portion 
# this was made before Tahmid worked on the input file (so there could be naming inconsistencies)
# also thinking about adding a funciton/script that would support looking at reverse complement ORFs
# right now I have ORF1, ORF2, ORF3 ect. and Incomplete_ORF1, Incomplete_ORF2, ect. 
# Stored in dictionaries separately.
# but I think I want to have it like this. I would love to know your idea on if this is a reasonable set up?
# {
#   "complete": {
#       "canonical": {  # ATG
#           "ORF1": {...},
#           "ORF2": {...}
#       },
#       "noncanonical": {
#           "GTG": {
#               "GTG_ORF1": {...}
#           },
#           "TTG": {
#               "TTG_ORF1": {...}
#           }
#       }
#   },
#   "incomplete": {
#       "canonical": {...},
#       "noncanonical": {
#           "GTG": {...},
#           "TTG": {...}
#       }
#   }
# }
"""
ORF_finder.py

Purpose:
    Detect open reading frames (ORFs) in a DNA sequence
    across all three forward reading frames.

Features:
    - Separates ORFs by start codon type (ATG, GTG, TTG)
    - Separates complete and incomplete ORFs
    - Labels ORFs numerically
"""
from typing import Dict

START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

def _scan_frame(dna_sequence: str, frame: int):
    """
    Scans a DNA sequence in a single reading frame for ORFs.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence to scan (assumes only A, T, G, C characters).
    frame : int
        Reading frame (0, 1, or 2) to scan.

    Returns
    -------
    generator of dict
        Each dictionary contains:
        - frame : int
        - start : int (start index of ORF)
        - end : int or None (end index of ORF if stop codon found)
        - start_type : str (start codon used)
        - status : str ('complete' or 'incomplete')
    """
    seq_len = len(dna_sequence)
    i = frame

    while i <= seq_len - 3:
        codon = dna_sequence[i:i + 3]
        
        if codon in START_CODONS:
            start_index = i
            start_type = codon
            j = i + 3
            found_stop = False

            while j <= seq_len - 3:
                stop_codon = dna_sequence[j:j + 3]

                if stop_codon in STOP_CODONS:
                    yield {
                        "frame": frame,
                        "start": start_index,
                        "end": j + 3,
                        "start_type": start_type,
                        "status": "complete"}
                    found_stop = True
                    break

                j += 3

            if not found_stop:
                yield {
                    "frame": frame,
                    "start": start_index,
                    "end": None,
                    "start_type": start_type,
                    "status": "incomplete"}

        i += 3

