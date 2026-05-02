"""
Nicole Decocker part 
orf_finder_lib — ORF detection package.

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
