# Title
Identifying Repeating Genes in a DNA Sequence Using Open Reading Frames (ORFs).

## Objective 
The project will develop a program that automate ORF detection and analysis in DNA sequences. It helps researchers identify potential protein-coding regions, repeated genes, and conserved sequences efficiently, while generating other ORF statistics for further study.

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
      * Input_validate.py: (Tahmid Anwar)
           - init.py 
           - fetch_fasta_from_ncbi.py
           - validate_DNA_sequence.py
             
      * ORF_finder_libarry: (Erin Nicole Decocker)
           - init.py 
           - Find_orfs.py 
           - Orfs_metadata.py

      * Statistics_summary_libarary: (Who ever gets done the soonest)
           - init.py 
           - Calculate_orf_stats.py
           - Write_stats_to_file.py

      * Orf_analysis_library.py: (Amanda Yaworsky)
           - init.py 
           - find_repeated_orfs.py
           - Calculate_similarity_scores.py 

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
