# ORF Analysis Pipeline
A Python program that detects Open Reading Frames (ORFs) in DNA sequences, analyzes repeated or similar ORFs, and produces summary statistics.

## Objective 
This project exists to automate ORF detection and analysis in DNA sequences. It helps researchers identify potential protein-coding regions, repeated genes, and conserved sequences efficiently, while generating other ORF statistics for further study.

## Project-Programming-2-Group-6
Github Repo Link: https://github.com/TahmidA139/Project-Programming-2-Group-6.git

# Project Overview
This project analyzes DNA sequences to detect Open Reading Frames (ORFs) and generate insights about their properties. 
The workflow:

* Fetch DNA sequences from NCBI in FASTA format.
* Scan sequences to detect all ORFs across reading frames (start codon → stop codon).
* Analyze ORFs to identify repeated or similar sequences, revealing potential gene duplications or conserved regions.
* Generate summary statistics, including total ORFs, repeated ORFs, and the longest ORF.

## Files and Descriptions
* main.py – Main driver for the ORF analysis pipeline. Holds main function which calls other modules and displays results.
* input_validate.py – Fetches DNA sequences from NCBI, validates, and writes cleaned FASTA files.
* orf_finder.py – Detects ORFs across reading frames and parses sequences.
* orf_analysis.py – Finds repeated ORFs and calculates pairwise similarity scores.
* Statistics_summary.py – Generates summary statistics (total ORFs, repeated, longest) and writes to a file.

# Installation Instructions
1. Clone the repository:
```bash
git clone https://github.com/TahmidA139/Project-Programming-2-Group-6.git
```
2. Ensure Python 3 is installed (≥3 recommended).

# License
This project is licensed under the GNU GPL v2.1. Chosen for open collaboration, ease of edits, and public use.

# Authors:
Tahmid Anwar
- Email: tanwar@charlotte.edu
- Student ID: 801501080

Amanda Yaworsky
- Student ID:801489950 
- Email: ayaworsk@charlotte.edu

Erin Nicole Decocker
- Student ID:801442694
- Email: edecocke@charlotte.edu
