#!/usr/bin/env python3
"""
(Whoever gets done the soonest and all contributes to this file)
statistics_summary.py


Purpose:
    Generate summary statistics from ORF datasets for reporting.
Input:
    Lists of ORF records and flagged ORFs.
Output:
    stats_summary.out containing counts and the longest ORF.
"""

def calculate_orf_stats(orf_list, flagged_orfs):
    """
    Purpose:
        Compute total ORFs, repeated ORFs, and the longest ORF.
    Input:
        orf_list (list), flagged_orfs (list)
    Output:
        Dictionary with 'total', 'repeated', 'longest' keys.
    High-level steps:
        1. Count all ORFs in orf_list.
        2. Count unique repeated ORFs in flagged_orfs.
        3. Identify ORF with maximum length in orf_list.
        4. Return all three stats as a dictionary.
    """
    pass

def write_stats_to_file(stats, outfile="stats_summary.out"):
    """
    Purpose:
        Save ORF statistics to a file.
    Input:
        stats (dict), outfile (str)
    Output:
        File containing total ORFs, repeated ORFs, and longest ORF.
    High-level steps:
        1. Open outfile for writing.
        2. Write each stat in a readable format.
        3. Close the file.
    """
    pass

