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
import argparse

from input_lib.input_validate import run as validate_run
from orf_finder_lib.ORF_finder import _scan_frame


def main():
    parser = argparse.ArgumentParser(description="Main ORF analysis pipeline")
    parser.add_argument("--accession", type=str, help="NCBI accession number")
    parser.add_argument("--email", type=str, help="Email for NCBI Entrez")
    args = parser.parse_args()
    
    # Prompt if not provided
    accession = args.accession or input("Enter NCBI accession number: ")
    email = args.email or input("Enter your email (required by NCBI): ")
    acc, clean_seq = validate_run(accession, email)

    if clean_seq is None:
        print("Pipeline failed.")
        return

    print("\n[INFO] Running ORF finder...\n")
    orfs = find_orfs(clean_seq, include_reverse=True)
    pprint(orfs)


if __name__ == "__main__":
    main()


