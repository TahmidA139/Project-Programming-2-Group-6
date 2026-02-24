#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ORF_finder.py

Purpose:
    Detect open reading frames (ORFs) in a DNA sequence
    across all three forward reading frames.

Features:
    - Supports bacterial alternative start codons (ATG, GTG, TTG)
    - Detects valid in-frame stop codons (TAA, TAG, TGA)
    - Separates complete and incomplete ORFs
    - Labels ORFs as ORF1, ORF2, etc.
"""

from typing import Dict


START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def _scan_frame(dna_sequence: str, frame: int):
    """
    Scan a single reading frame and return raw ORF data.
    """
    seq_len = len(dna_sequence)
    i = frame

    while i <= seq_len - 3:
        codon = dna_sequence[i:i + 3]

        if codon in START_CODONS:
            start_index = i
            j = i + 3
            found_stop = False

            while j <= seq_len - 3:
                stop_codon = dna_sequence[j:j + 3]

                if stop_codon in STOP_CODONS:
                    yield {
                        "frame": frame,
                        "start": start_index,
                        "end": j + 3,
                        "status": "complete"
                    }
                    found_stop = True
                    break

                j += 3

            if not found_stop:
                yield {
                    "frame": frame,
                    "start": start_index,
                    "end": None,
                    "status": "incomplete"
                }

        i += 3


def find_orfs(dna_sequence: str) -> Dict:
    """
    Detect ORFs and organize them into labeled dictionaries.
    """
    dna_sequence = dna_sequence.upper()

    complete_orfs = {}
    incomplete_orfs = {}

    complete_count = 1
    incomplete_count = 1

    for frame in range(3):
        for orf in _scan_frame(dna_sequence, frame):

            start = orf["start"]
            end = orf["end"]

            if end is not None:
                sequence = dna_sequence[start:end]
                length = len(sequence)

                complete_orfs[f"ORF{complete_count}"] = {
                    "frame": frame,
                    "start": start,
                    "end": end,
                    "length": length,
                    "sequence": sequence
                }

                complete_count += 1

            else:
                incomplete_orfs[f"Incomplete_ORF{incomplete_count}"] = {
                    "frame": frame,
                    "start": start,
                    "end": None,
                    "length": None,
                    "sequence": None,
                    "note": "No in-frame stop codon detected"
                }

                incomplete_count += 1

    return {"complete_orfs": complete_orfs, "incomplete_orfs": incomplete_orfs}
