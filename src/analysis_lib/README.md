# ORF Analysis Module

This Python module provides basic tools for analyzing Open Reading Frames (ORFs). It focuses on identifying repeated ORFs and calculating simple similarity scores between ORF sequences.

## nicoles commment:
**Note on calculate_similarity_scores:**
- We are dropping the position-by-position similarity score comparison between every ORF pair. The reason is that comparing every ORF to every other ORF scales very poorly. Like it means for a sequence with possibly hundreds of ORFs this becomes thousands of comparisons that are also not biologically meaningful without lotsssss more functions and analysis. So instead we think it would be best to replace similiarity scores   with codon usage analysis (simpler math analysis) and comparative statistics between two sequences (comp stats will be in the stats files not ur stuff), which is more informative and computationally reasonable.

**Functions that need to be in this file:**
- **find_repeated_orfs** — you have this! you do need to add list of ORF dicts instead of raw sequences tho 
- **calculate_orf_stats** — Loops over every ORF and computes per-ORF statistics (GC content, protein length, codon usage). Calls the helper functions below to do the actual math (helpers keep functions less than 40 lines lol)
- **_extract_sequence** — pulls the actual nucleotide sequence of each ORF out of the full DNA string. Needs to handle both forward and reverse strand ORFs separately
-  **_gc_content** — calculates what percentage of a sequence is G or C bases. Called by calculate_orf_stats for each ORF
- **_codon_usage** — counts how many times each codon appears in an ORF sequence. Called by calculate_orf_stats for each ORF
- **_protein_length** — estimates how many amino acids the ORF would produce. Complete ORFs subtract the stop codon, incomplete ORFs use the full length
- **_tally_by_frame** — counts how many ORFs are in each reading frame (e.g. there are six reading frames we looked at: +1, +2, +3, -1, -2, -3). 
- **_length_stats** — computes min, max, and mean ORF length across all ORFs. 
- **_codon_usage_totals** — aggregates codon usage counts across all ORFs into one summary dictionary. 
- **_top_codons** — takes the aggregated codon usage and returns just the top 10 most used codons. 

If you want me to do this part and then you do the stats file which will legit just be to make the two txt summary output files with the comparison let me know. If you want to keep your part let me know and ill start working on the stats file. The stats file will only need to have three functions:

- **write_stats_to_file**  needs tp write a human-readable summary report to a text file that includes dataset-level counts, genomic GC content, longest ORF details, and a per-ORF table with GC content, protein length, and codon usage.

- **write_comparative_report** Write a human-readable side-by-side comparative report for two sequences.

- **write_comparative_csv** Write a codon-usage delta table to a CSV file this will be used for the grapghics part which will be our last module!!

## Features

- Find Repeated ORFs
  - Identifies ORF sequences that appear more than once in a list.
  - Returns a dictionary with ORFs and their counts.

- Calculate Similarity Scores
  - Compares each ORF to every other ORF.
  - Computes a similarity score based on matching characters at the same positions.
  - Uses the length of the shorter sequence for comparison.

## Functions

### find_repeated_orfs(orfs)

Input:
- orfs (list): A list of ORF sequences (strings)

Output:
- dict: ORFs that appear more than once, with their counts

Example:
orfs = ["ATGAAA", "ATGAAA", "ATGCCA", "ATGAAA"]
print(find_repeated_orfs(orfs))

Output:
{'ATGAAA': 3}

---

### calculate_similarity_scores(orfs)

Input:
- orfs (list): A list of ORF sequences (strings)

Output:
- dict: Pairwise similarity scores using index pairs as keys

How it works:
- Compares each pair of sequences
- Counts matching positions
- Divides by the length of the shorter sequence

Example:
orfs = ["ATGAAA", "ATGAAA", "ATGCCA", "ATGAAA"]
print(calculate_similarity_scores(orfs))

Output:
{
 (0, 1): 1.0,
 (0, 2): 0.6666666666666666,
 (0, 3): 1.0,
 (1, 2): 0.6666666666666666,
 (1, 3): 1.0,
 (2, 3): 0.6666666666666666
}

## Notes

- Similarity is position-based and does not perform sequence alignment.
- Scores range from 0.0 (no matches) to 1.0 (identical sequences).
- The function assumes sequences are non-empty; empty sequences may cause errors.

## Requirements

- Python 3.x

## Usage

Run the script directly:

python ORF_analysis.py
