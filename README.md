# ORCA (ORF Recognition and Classification Analyzer)
**GitHub:** https://github.com/TahmidA139/Project-Programming-2-Group-6

## Objective 
The project will develop a program that automate ORF detection and analysis in DNA sequences. It helps researchers identify potential protein-coding regions, repeated genes, and conserved sequences efficiently, while generating other ORF statistics for further study.

## Features
- Downloads DNA sequences in FASTA format from NCBI
- Scans sequences to identify start and stop codons across all reading frames
- Detects repeated or similar ORF sequences
- Generates summary statistics: total ORFs, repeated ORFs, and longest ORF

### Files and Descriptions
* main.py
 – This holds main function which calls other modules and displays results.
* input_validate.py
– This fetches DNA sequences from NCBI, validates, and writes cleaned FASTA files.
* ORF_finder.py
– Detects ORFs across reading frames and parses sequences.
* ORF_analysis.py
– Finds repeated ORFs and calculates pairwise similarity scores.
* statistics_summary.py
– Generates summary statistics (total ORFs, repeated, longest) and writes to a output file.

# Task distribution among members:
  - Libraries: 
      * Input_validate_lib/ (Tahmid Anwar)  — Input validation and FASTA file loading
           - init.py 
           - fetch_fasta_from_ncbi.py
           - validate_DNA_sequence.py
             
      * ORF_finder_lib/ (Erin Nicole Decocker) — Core ORF detection logic
           - init.py 
           - Find_orfs.py 
           - Orfs_metadata.py
    
      * analysis_lib/ (Amanda Yaworsky)  — ORF analysis and repeat detection
           - init.py 
           - find_repeated_orfs.py
           - Calculate_similarity_scores.py 

      * statistics_lib/ (Who ever gets done the soonest) — Summary statistics
           - init.py 
           - Calculate_orf_stats.py
           - Write_stats_to_file.py


# Installation Instructions
1. Clone the repository:
```bash
git clone https://github.com/TahmidA139/Project-Programming-2-Group-6.git
```
2. Ensure Python 3 is installed (≥3 recommended).

# License
This project is licensed under the GNU GPL v2.1. Chosen for open collaboration, ease of edits, and public use.

# Authors:
Amanda Yaworsky
- Student ID:801489950 
- Email: ayaworsk@charlotte.edu

Erin Nicole Decocker
- Student ID:801442694
- Email: edecocke@charlotte.edu

Tahmid Anwar
- Email: tanwar@charlotte.edu
- Student ID: 801501080

## Project-Programming-2-Group-6
Github Repo Link: https://github.com/TahmidA139/Project-Programming-2-Group-6.git

## 20260226_input_validate.py_function_execution
We created a function called fetch_fasta_from_ncbi which will output raw DNA sequence string from the NCBI nucleotide accession number. Using try/except and inside it using Entrez.efetch, using the accession number, we collect the fasta file of that accession number. Then using SeqIO.read(), we take only the raw sequence as string ignoring the headers of the fasta file. The cleaned DNA sequence will be passed directly to ORF_finder.py through main.py

