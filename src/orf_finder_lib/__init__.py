"""
orf_finder_lib — ORF detection package.

Nicole Decocker's part 

What this folder contains:
    frame_scanner.py — low-level DNA scanning using NumPy vectorization.
                       Handles reverse complement, codon arrays, and
                       finding stop codons. Called by orf_finder.py only.
    orf_finder.py    — high-level ORF detection orchestrator. Scans all
                       six reading frames and assembles results into a
                       nested dict and flat list for CSV output.

Tells Python that orf_finder_lib/ is a package so other modules can import from it:
    from src.orf_finder_lib.orf_finder import find_orfs
    from src.orf_finder_lib.frame_scanner import scan_frame
"""

# ORF finder
# High-level ORF detection from orf_finder.py
# find_orfs() — main entry point; scans all six reading frames and returns
#               a nested dict (by start-codon type) and a flat list of ORFs
from src.orf_finder_lib.orf_finder import find_orfs

# Frame scanner
# Low-level scanning utilities from frame_scanner.py
# scan_frame()         — scans a single reading frame for ORFs above the length threshold
# reverse_complement() — returns the reverse complement of a DNA sequence
from src.orf_finder_lib.frame_scanner import (
    scan_frame,
    reverse_complement,
)

# __all__ controls what gets exported when someone does: from src.orf_finder_lib import *
__all__ = [
    "find_orfs",
    "scan_frame",
    "reverse_complement",
]
