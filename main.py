#!/usr/bin/env python3
"""
(Everyone in the group will participate in this portion)
main.py
Purpose:
    This script serves as the main driver for the ORF analysis pipeline.
Role in Project:
    This file functions by calling other modules to handle sequence loading, ORF detection, analysis, and reporting.
Input:
    DNA sequence file in FASTA format from NCBI
Output:
    Summary of ORFs found and their statistics.
"""

from pprint import pprint
from 20260226_input_validate import run as validate_run
from ORF_finder import find_orfs


def main():
    accession = "NM_001301717"

    acc, clean_seq = validate_run(accession)

    if clean_seq is None:
        print("Pipeline failed.")
        return

    print("\n[INFO] Running ORF finder...\n")
    orfs = find_orfs(clean_seq, include_reverse=True)

    pprint(orfs)


if __name__ == "__main__":
    main()


