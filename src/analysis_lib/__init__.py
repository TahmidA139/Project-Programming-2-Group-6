"""
analysis_lib — ORF analysis and statistics package.


Amanda Yaworsky part  


What this folder contains:
    orf_analysis.py      — identifies repeated ORFs and calculates pairwise
                           similarity scores, GC content, and codon usage
                           for each detected ORF.
    statistics_summary.py— generates a comprehensive summary report of all
                           ORFs found, including total counts, longest ORF,
                           strand distribution, and writes orf_summary.txt.

Tells Python that analysis_lib/ is a package so other modules can import from it:
    from src.analysis_lib.orf_analysis import analyse_orfs
    from src.analysis_lib.statistics_summary import summarise
"""

# ORF analysis
# Per-ORF statistics and repeated-ORF detection from orf_analysis.py
# gc_content()          — GC content of a sequence as a percentage
# protein_length()      — number of translated codons in a sequence
# codon_usage()         — codon-frequency dict for a sequence
# calculate_orf_stats() — bundles gc_content, protein_length, and codon_usage per ORF
# find_repeated_orfs()  — detects ORF sequences that appear more than once
from src.analysis_lib.orf_analysis import (
    gc_content,
    protein_length,
    codon_usage,
    calculate_orf_stats,
    find_repeated_orfs,
)

# Statistics and reporting
# File writing and console reporting functions from statistics_summary.py
# print_summary()              — prints ORF counts to the console
# write_stats_to_file()        — writes orf_summary.txt for a single sequence
# write_orf_comparison_report()— writes the combined report for comparative mode
# write_gff3()                 — writes a GFF3 annotation file for a sequence
from src.analysis_lib.statistics_summary import (
    print_summary,
    write_stats_to_file,
    write_orf_comparison_report,
    write_gff3,
)

# __all__ controls what gets exported when someone does: from src.analysis_lib import *
__all__ = [
    "gc_content",
    "protein_length",
    "codon_usage",
    "calculate_orf_stats",
    "find_repeated_orfs",
    "print_summary",
    "write_stats_to_file",
    "write_orf_comparison_report",
    "write_gff3",
]
