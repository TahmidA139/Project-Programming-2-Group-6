"""
analysis_lib — ORF analysis and statistics package.

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
