#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Nicole's portion 
# this was made before Tahmid worked on the input file (so there could be naming inconsistencies)
# also thinking about adding a funciton that would support looking at reverse complement ORFs
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
    Generator that yields raw ORF info including start codon type.
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


def find_orfs(dna_sequence: str) -> Dict:
    """
    Detect ORFs and organize them by start codon type
    and completion status.
    """
    dna_sequence = dna_sequence.upper()

    # Complete ORFs
    atg_orfs = {}
    gtg_orfs = {}
    ttg_orfs = {}

    # Incomplete ORFs
    incomplete_atg = {}
    incomplete_gtg = {}
    incomplete_ttg = {}

    counters = {
        "ATG": 1,
        "GTG": 1,
        "TTG": 1,
        "I_ATG": 1,
        "I_GTG": 1,
        "I_TTG": 1,}

    for frame in range(3):
        for orf in _scan_frame(dna_sequence, frame):

            start = orf["start"]
            end = orf["end"]
            start_type = orf["start_type"]

            if end is not None:
                sequence = dna_sequence[start:end]
                length = len(sequence)

                label = f"{start_type}_ORF{counters[start_type]}"
                counters[start_type] += 1

                data = {
                    "frame": frame,
                    "start": start,
                    "end": end,
                    "length": length,
                    "sequence": sequence
                }

                if start_type == "ATG":
                    atg_orfs[label] = data
                elif start_type == "GTG":
                    gtg_orfs[label] = data
                elif start_type == "TTG":
                    ttg_orfs[label] = data

            else:
                label = f"Incomplete_{start_type}_ORF{counters['I_' + start_type]}"
                counters["I_" + start_type] += 1

                data = {
                    "frame": frame,
                    "start": start,
                    "end": None,
                    "length": None,
                    "sequence": None,
                    "note": "No in-frame stop codon detected"
                }

                if start_type == "ATG":
                    incomplete_atg[label] = data
                elif start_type == "GTG":
                    incomplete_gtg[label] = data
                elif start_type == "TTG":
                    incomplete_ttg[label] = data

    return {
        "ATG_complete": atg_orfs,
        "GTG_complete": gtg_orfs,
        "TTG_complete": ttg_orfs,
        "ATG_incomplete": incomplete_atg,
        "GTG_incomplete": incomplete_gtg,
        "TTG_incomplete": incomplete_ttg,
    }
