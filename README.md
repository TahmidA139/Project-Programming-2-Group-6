<p align="center">
  <img src="pseudocode/orca_banner.png" alt="ORCA Banner" width="350"/>
</p>

# ORCA — ORF Recognition and Comparative Annotator

Tahmid Anwar's Part
---

## Description

ORCA is a Python-based bioinformatics pipeline that automates Open Reading Frame (ORF) detection and annotation across DNA sequences. It accepts an NCBI accession number or a local FASTA file, fetches or loads the DNA sequence, detects all ORFs across all six reading frames, and computes per-ORF statistics including GC content and protein length. When two sequences are provided, ORCA performs side-by-side comparative analysis including shared ORF detection and codon usage comparison (RSCU).

---

## Features

- Downloads DNA sequences directly from NCBI using an accession number, or loads them from a local FASTA file
- Validates and cleans input sequences (handles IUPAC ambiguity codes, whitespace, and invalid characters)
- Scans all six reading frames (+1, +2, +3, −1, −2, −3) using NumPy vectorization
- Detects canonical (ATG) and non-canonical (GTG, TTG) start codons
- Computes per-ORF statistics: GC content, protein length (number of complete codons), and strand/frame information
- Exports ORF annotations in GFF3 format for use with genome browsers
- **Comparative mode** (`--accession2` / `--fasta2`): side-by-side ORF statistics, shared ORF sequence detection, and an RSCU codon-usage heatmap
- Generates an ORF map visualisation across all six reading frames (PNG)
- Writes all results to the `output/` folder as GFF3 files, CSV files, and plain-text reports

---

## Requirements

### Python Version

- Python **3.10**

### Packages

| Package | Purpose |
|---------|---------|
| `numpy` | Vectorized codon array construction and ORF scanning |
| `matplotlib` | ORF map and RSCU heatmap figure generation |
| `biopython` | NCBI Entrez sequence fetching and FASTA parsing |

Standard library modules used (`argparse`, `os`, `sys`, `re`, `csv`, `datetime`) require no separate installation.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/TahmidA139/ORCA.git
cd ORCA
```

### 2. Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate ORCA
```

The `environment.yml` file was exported with `--from-history` so conda will resolve the correct platform-specific sub-dependencies on any operating system.

### 3. Verify the installation

```bash
python -m src.main --help
```

You should see the full argument list printed to the terminal. If you see an import error, confirm that the `ORCA` environment is active and that all packages installed without errors.

---

## Usage

Run all commands from the **project root directory** with the `ORCA` environment activated.

### Input source rules

- Sequence 1: provide **either** `--accession` **or** `--fasta`, not both.
- Sequence 2: provide **either** `--accession2` **or** `--fasta2`, not both.
- Local FASTA files must contain **exactly one** sequence record. Files with multiple records will cause the pipeline to exit with a descriptive error message.

### Example 1 — Single sequence from NCBI (all defaults)

Uses the default start codon (ATG only) and minimum ORF length (30 nt):

```bash
python -m src.main \
    --accession <ACCESSION_NUMBER> \
    --email     <YOUR_EMAIL>
```

### Example 2 — Single sequence from a local FASTA file with custom settings

Searches for all three start codons and requires ORFs to be at least 150 nt:

```bash
python -m src.main \
    --fasta        <PATH_TO_FASTA> \
    --email        <YOUR_EMAIL> \
    --start-codons ATG GTG TTG \
    --min-length   150
```

### Example 3 — Comparative mode with two NCBI accession numbers

Runs the full comparative pipeline including shared ORF detection and RSCU heatmap:

```bash
python -m src.main \
    --accession  <ACCESSION_NUMBER_1> \
    --accession2 <ACCESSION_NUMBER_2> \
    --email      <YOUR_EMAIL>
```

Please see ENVIRONMENT.md for more examples!

Replace each placeholder with your own values:

| Placeholder | Description |
|-------------|-------------|
| `<ACCESSION_NUMBER>` | NCBI nucleotide accession number, e.g. `NM_012367.1` |
| `<ACCESSION_NUMBER_1>` | NCBI accession for sequence 1 (comparative mode) |
| `<ACCESSION_NUMBER_2>` | NCBI accession for sequence 2 (comparative mode) |
| `<PATH_TO_FASTA>` | Path to a single-sequence FASTA file, e.g. `example_input_files/OR2B6_sequence.fasta` |
| `<PATH_TO_FASTA_1>` | Path to the first FASTA file (comparative mode) |
| `<PATH_TO_FASTA_2>` | Path to the second FASTA file (comparative mode) |
| `<YOUR_EMAIL>` | A valid email address, e.g. `you@example.com` (required by NCBI Entrez) |

---

## Command-Line Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--accession` | See note | — | NCBI accession number for sequence 1. Cannot be used with `--fasta`. |
| `--fasta` | See note | — | Path to a local FASTA file for sequence 1. Must contain exactly one sequence. Cannot be used with `--accession`. |
| `--accession2` | No | — | NCBI accession number for sequence 2. Enables comparative mode. |
| `--fasta2` | No | — | Path to a local FASTA file for sequence 2. Enables comparative mode. Cannot be used with `--accession2`. |
| `--email` | Yes | — | Email address required by NCBI Entrez. Must be supplied even when loading from local files. |
| `--start-codons` | No | `ATG` | One or more start codons: `ATG`, `GTG`, `TTG`. |
| `--min-length` | No | `30` | Minimum ORF length in nucleotides. Must be at least 3. |
| `--outdir` | No | `output/` | Directory for all output files. Created automatically if it does not exist. |

**Note:** Exactly one of `--accession` or `--fasta` must be provided for sequence 1.

---

## Output Format

All output files are written to `output/` by default. A custom directory can be set with `--outdir`.

### Single-sequence mode

| File | Description |
|------|-------------|
| `output/cleaned_sequence_1.fasta` | Input sequence after validation and cleaning. |
| `output/<accession>.gff3` | ORF annotations in GFF3 format (e.g. `NM_012367.1.gff3`). |
| `output/orf_summary.txt` | Human-readable summary of all ORF statistics. |
| `output/orf_map.png` | ORF map across all six reading frames. |

### Comparative mode (additional files)

| File | Description |
|------|-------------|
| `output/comp_cleaned_sequence_1.fasta` | Cleaned sequence 1. |
| `output/comp_cleaned_sequence_2.fasta` | Cleaned sequence 2. |
| `output/<accession1>.gff3` | ORF annotations for sequence 1 in GFF3 format. |
| `output/<accession2>.gff3` | ORF annotations for sequence 2 in GFF3 format. |
| `output/orf_comparison_report.txt` | Full comparative report (see format below). |
| `output/codon_usage_comparison.png` | RSCU codon-usage heatmap comparing both sequences. |
| `output/orf_map.png` | Two-panel comparative ORF map. |

### Comparison report format (`orf_comparison_report.txt`)

The report is divided into three sections:

**1. Run parameters** — date generated, accessions analysed, start codons used, and minimum ORF length.

**2. Per-sequence statistics** — for each sequence: total ORF count, average GC content, the longest ORF (ID, length, strand, frame, GC%), and a full per-ORF table:

```
#    ID        Length    GC%   Prot_len  Strand  Frame
────────────────────────────────────────────────────────
0    ORF1         942  43.95       314       +      0
1    ORF2          84  42.86        28       -      0
```

**3. Comparative summary** — total ORFs per sequence, forward/reverse strand breakdown, number of shared ORF sequences, and counts unique to each input.

### Example terminal output

```
[ORCA] Processing sequence 1: NM_001838.4

[VALIDATION] Sequence is valid — 2164 bp ready for analysis.

════════════════════════ ORF Summary — Sequence 1 ════════════════════════
  Total ORFs found  : 27
  Forward strand (+): 14
  Reverse strand (-): 13
  Canonical (ATG)   : 27
────────────────────────────────────────────────────────────────────────

[ORCA] Processing sequence 2: NM_012367.1

[VALIDATION] Sequence is valid — 942 bp ready for analysis.

════════════════════════ ORF Summary — Sequence 2 ════════════════════════
  Total ORFs found  : 14
  Forward strand (+): 7
  Reverse strand (-): 7
  Canonical (ATG)   : 14
────────────────────────────────────────────────────────────────────────

[INFO] ORF map saved to: output/orf_map.png
[INFO] Comparison report written to: output/orf_comparison_report.txt
```

---

## Project Structure

```
ORCA/
├── README.md
├── ENVIRONMENT.md
├── LICENSE
├── idea_and_task_distribution.txt
├── environment.yml
├── run_test.sh
│
├── example_input_files/
│   ├── homo_sapiens_albumin.fasta
│   ├── IUPAC_ambiguity_test.fasta
│   ├── multiple_sequence_file.fasta
│   ├── mus_musculus_albumin.fasta
│   └── OR2B6_sequence.fasta
│
├── example_output/
│   ├── default_comparative_run/      # Example outputs for comparative mode
│   └── default_single_sequence_run/  # Example outputs for single-sequence mode
│
├── pseudocode/
│   ├── flowcharts.txt
│   ├── pseudocode.txt
│   └── *.png                         # Multiple pipeline flowchart diagrams
│
└── src/
    ├── __init__.py
    ├── main.py                       
    ├── graphics.py                    
    ├── input_validate.py          
    │
    ├── analysis_lib/                  
    │   ├── __init__.py
    │   ├── orf_analysis.py           
    │   └── statistics_summary.py     
    │
    └── orf_finder_lib/          
        ├── __init__.py
        ├── frame_scanner.py          
        └── orf_finder.py              
```

Each subdirectory is a proper Python package (`__init__.py` present) so the pipeline is invoked as a module: `python -m src.main`.

---

## Algorithm Description

### Input validation and cleaning (`input_validate.py`)

1. The raw sequence is uppercased and all whitespace and digits are stripped.
2. Characters outside the full IUPAC nucleotide alphabet (A T G C R Y S W K M B D H V N) are removed.
3. Remaining IUPAC ambiguity codes (anything that is not A, T, G, or C) are replaced with `N`. Replacement rather than deletion preserves the reading frame — deleting bases would shift every downstream codon.
4. Sequences shorter than 6 bp are rejected.

### Six-frame ORF scanning (`frame_scanner.py`, `orf_finder.py`)

1. The forward-strand sequence is converted to a NumPy codon array for each frame offset (0, 1, 2).
2. The same process is repeated on the reverse complement to cover frames −1, −2, −3.
3. Within each frame, all positions matching a requested start codon are identified using a Boolean mask.
4. For each start codon position, the first downstream stop codon (TAA, TAG, TGA) is found. If no stop codon exists, the ORF is discarded (open ORFs are not reported — see Known Limitations).
5. When multiple start codons share the same downstream stop codon, only the longest ORF (earliest start codon) is retained.
6. ORFs shorter than `--min-length` are discarded.
7. Reverse-complement coordinates are converted back to the forward-strand reference frame before storage.

### Per-ORF statistics (`orf_analysis.py`)

For each ORF the following are computed:

- **GC content** — counted on the coding region only (stop codon excluded): `(G + C) / length × 100`.
- **Protein length** — number of complete codons in the coding region: `len(coding_seq) // 3`.
- **Nucleotide sequence** — extracted 5'→3' from the forward-strand sequence; reverse-complement strand ORFs are reverse-complemented before storage.
- **Codon usage** — non-overlapping triplet counts across the coding region; codons containing `N` are skipped.

### Metadata and calculated fields

Each ORF record stored in the output CSV and GFF3 file contains the following metadata fields:

| Field | Description |
|-------|-------------|
| `orf_id` | Human-readable label (e.g. `ORF1`, `GTG_ORF1`). Canonical ORFs are numbered sequentially; non-canonical ORFs are prefixed with their start codon. |
| `strand` | `+` for the forward strand, `-` for the reverse complement strand. |
| `frame` | Reading frame offset: 0, 1, or 2. Combined with `strand` this uniquely identifies one of the six frames. |
| `start` | 0-based inclusive start coordinate on the forward strand. |
| `end` | 0-based exclusive end coordinate on the forward strand (GFF3 export converts these to 1-based fully-closed coordinates automatically). |
| `length_nt` | Total ORF length in nucleotides, including the stop codon. |
| `start_codon` | The three-letter start codon found at the ORF's 5' end (e.g. `ATG`, `GTG`, `TTG`). |
| `gc_content` | GC percentage of the coding region (stop codon excluded), as a float rounded to two decimal places. |
| `protein_length` | Number of translated codons (amino acids) in the coding region, computed as `len(coding_seq) // 3`. |
| `status` | Always `complete` in the current release; reserved for future partial-ORF support. |

### RSCU — Relative Synonymous Codon Usage

RSCU measures how often a particular codon is used relative to what would be expected if all synonymous codons (those encoding the same amino acid) were used equally. It is calculated as:

```
RSCU = observed_count / (total_codons_for_amino_acid / number_of_synonymous_codons)
```

- An RSCU of **1.0** means a codon is used exactly as often as expected under equal usage.
- An RSCU **> 1.0** means the codon is used more often than expected (preferred).
- An RSCU **< 1.0** means the codon is used less often than expected (avoided).

In ORCA, RSCU values are computed across all ORF sequences for each sequence in a comparative run and displayed as a heatmap (`codon_usage_comparison.png`). Columns correspond to codons (grouped by amino acid) and rows correspond to the two input sequences, making biases and differences in translational preference immediately visible.

---

## References

- Cock, P. J. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://biopython.org/
- Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. https://numpy.org/
- Sayers, E. W., et al. (2022). Database resources of the National Center for Biotechnology Information. *Nucleic Acids Research*, 50(D1), D20–D26. https://www.ncbi.nlm.nih.gov/books/NBK25499/
- Brent, M. R., & Shi, L. (2024). ORF annotation and the challenge of small proteins. *BMC Genomics*, 25, 1016. https://pmc.ncbi.nlm.nih.gov/articles/PMC11521203/
- Sharp, P. M., & Li, W.-H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. (RSCU concept)
- NCBI ORF Finder: https://www.ncbi.nlm.nih.gov/orffinder/
- IUPAC nucleotide code reference: https://www.bioinformatics.org/sms/iupac.html
- GFF3 specification: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- FASTA format reference: https://en.wikipedia.org/wiki/FASTA_format
- Open Reading Frame overview: https://en.wikipedia.org/wiki/Open_reading_frame
- Python argparse documentation: https://docs.python.org/3/library/argparse.html
- Python typing module documentation: https://docs.python.org/3/library/typing.html
- PEP 8 — Style Guide for Python Code: https://peps.python.org/pep-0008/

---

## Known Limitations

**Open ORFs are not reported.** ORCA only reports complete ORFs, those that contain both a start codon and an in-frame downstream stop codon within the input sequence. ORFs that begin near the end of a sequence and have no stop codon before the sequence terminates are silently discarded. This is most relevant when the input is a genomic fragment rather than a complete coding sequence. If your sequence is likely to contain partial genes at either end, be aware that those ORFs will not appear in the output.

**Only the longest ORF per stop codon is reported.** When multiple start codons share the same downstream stop codon in a given reading frame, ORCA reports only the longest one (i.e. the one with the earliest start codon). Shorter nested ORFs that begin downstream of that first start codon but terminate at the same stop are not included. This follows the standard longest-ORF convention used by most prokaryotic ORF finders, but may result in missed internal start sites in eukaryotic sequences where alternative translation initiation is biologically relevant.

---

## AI Usage Statement

**Tahmid Anwar, Amanda Yaworsky, and Nicole Decocker** acknowledge the use of **AI (Claude, Anthropic)** as a supporting tool during development of the ORCA pipeline. Below is a description of how it was used along with representative example prompts. All AI-generated suggestions were reviewed, tested, and edited by the team before being accepted into the codebase.

**1. Code review and duplicate detection**

We used Claude to identify duplicate or redundant functions spread across multiple files in the project, helping us consolidate logic and avoid inconsistencies between modules.

Example prompts:
- *"Here are two Python files from our project. Are there any functions that are duplicated or doing the same thing across both files?"*
- *"We have a helper function in `input_validate.py` and a similar one in `orf_analysis.py` — are these redundant and can they be merged?"*

**2. Type hints and documentation**

We used Claude to add and improve type hints across our functions and to clean up docstrings so they were consistent in style and complete.

Example prompts:
- *"Add proper Python type hints to all the functions in this file."*
- *"Rewrite the docstrings in this module to follow a consistent format with Input, Output, and How it works sections."*

**3. Debugging and troubleshooting**

We used Claude to diagnose errors and unexpected behavior encountered when running the pipeline, sharing terminal output and code to identify the root cause.

Example prompts:
- *"Here is the error message and the function producing it. What is going wrong and how do I fix it?"*
- *"Our pipeline runs without crashing but the output CSV is empty. Here is the relevant code — what could cause this?"*

**4. Visualization development**

We used Claude to help design and debug the codon usage comparison plot, and to help create an outline of the flowcharts including fixing layout issues and ensuring the output PNG rendered correctly.

Example prompts:
- *"Here is our codon usage comparison function and the plot it produces. The bars are overlapping and the legend is cut off — how do we fix this?"*
- *"Write a matplotlib function that takes two DNA sequences and plots a side-by-side codon usage bar chart, saving it as a PNG."*

---

## License

This project is licensed under the **GNU Lesser General Public License v2.1 (LGPL-2.1)** which was chosen for open collaboration, ease of contribution, and public use. See the `LICENSE` file in the repository root for the full license text.

---

## Authors

**Amanda Yaworsky** — Student ID: 801489950 · GitHub: [amandayaworsky](https://github.com/amandayaworsky) · ayaworsk@charlotte.edu

**Erin Nicole Decocker** — Student ID: 801442694 · GitHub: [edecocke-uncc](https://github.com/edecocke-uncc) · edecocke@charlotte.edu

**Tahmid Anwar** — Student ID: 801501080 · GitHub: [TahmidA139](https://github.com/TahmidA139) · tanwar@charlotte.edu
