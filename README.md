<p align="center">
  <img src="images/orca_banner.png" alt="ORCA Banner" width="350"/>
</p>

# ORCA (ORF Recognition and Comparative Annotator)

## Objective
ORCA is a python-based bioinformatics pipeline that automates Open Reading Frame(ORF) detection and analysis in DNA sequences. It helps researchers identify potential 
protein-coding regions, repeated genes, and conserved sequences efficiently, while generating summary statistics and visualizations for further study.

The pipeline accepts an NCBI accession number as input, fetches the DNA sequence, validates it, detects all ORFs in the sequence, analyzes repeats and 
similarity, and produces statistical summaries and plots.

## Features
- Downloads DNA sequences in FASTA format directly from NCBI using an accession number
- Validates and cleans sequences before analysis
- Scans all six reading frames (+1, +2, +3, −1, −2, −3) using NumPy vectorization
- Detects canonical (ATG) and non-canonical (GTG, TTG) start codons
- Identifies nested ORFs
- Computes per-ORF statistics: GC content, codon usage, and protein length
- **Comparative mode** (`--accession2`): side-by-side ORF structure comparison, codon usage differences, and other analysis between two sequences
- Writes all results to the `output/` folder as .CSV and plain-text reports

## Project Structure

```
ORCA/
├── README.md
├── LICENSE
├── environment.yml
├── src/
│   ├── __init__.py
│   ├── main.py                        
│   │
│   ├── input_lib/                     # Tahmid Anwar
│   │   ├── __init__.py
│   │   └── input_validate.py         
│   │
│   ├── orf_finder_lib/                # Erin Nicole Decocker
│   │   ├── __init__.py
│   │   └── frame_scanner.py   
│   │   └── orf_finder.py
│   │   └── output_writer.py                
│   │
│   ├── analysis_lib/                  # Amanda Yaworsky
│   │   ├── __init__.py
│   │   └── orf_analysis.py            
│   │   └── statistics_summary.py                                    
│   │
│   │
│   ├── graphics_lib/                  # Erin Nicole Decocker 
│   │   ├── __init__.py
│   │   └── graphics.py                                        
│
└── examples/                          # Tahmid Anwar
    ├── example_output.fasta
    └── example_run.txt
```

### File Descriptions
### 0. main.py
The entry point of the entire project. It coordinates all modules by calling them
in the correct order and passing outputs from one module as inputs to the next.
You never need to run `input_validate.py`, `ORF_finder.py`, `ORF_analysis.py`, or
`statistics_summary.py` individually — running `main.py` executes the full pipeline
from start to finish in one command.

**What main.py does step by step:**
```
1. Parses command line arguments (--accession, --email, --min-length, --ignore-nested)
2. Calls input_validate.py  → fetches and validates the DNA sequence from NCBI
3. Calls ORF_finder.py      → detects all ORFs across all 6 reading frames
4. Calls ORF_analysis.py    → finds repeated ORFs and calculates similarity scores
5. Calls statistics_summary.py → generates summary statistics and plots
6. Displays final results in the terminal
```

**Input:** Command line arguments provided by the user

**Output:** Coordinates and triggers all output files from every module

**How to run:**
```bash

# Full run with all options
python -m src.main --accession NM_001301717 --email your_email@example.com --min-length 75 --ignore-nested
```
### Full command with all options explained
```bash
python                        # run Python
    -m src.main               # run main.py inside the src/ folder
    --accession NM_001301717  # NCBI accession number of the DNA sequence to analyse
    --email your_email@example.com  # your email address required by NCBI Entrez API
    --min-length 75           # only report ORFs that are at least 75 base pairs long
    --ignore-nested           # skip ORFs that are fully contained inside a longer ORF
```
**How main.py connects all modules:**
```
main.py
  │
  ├── input_validate.run()        accession number
  │       │                    ──────────────────→  cleaned DNA sequence string
  │       │
  ├── ORF_finder.find_all_orfs()  cleaned DNA sequence
  │       │                    ──────────────────→  list of ORF dictionaries + orfs.csv
  │       │
  ├── ORF_analysis.analyse_orfs() list of ORF dictionaries
  │       │                    ──────────────────→  repeat results + orf_analysis.csv
  │       │
  └── statistics_summary.summarise() ORFs + analysis results
                               ──────────────────→  summary .txt + plots
```

---

### 1. input_validate.py
Fetches a DNA sequence from NCBI using the accession number, validates it, removes
any invalid characters, and writes a cleaned FASTA file to disk.

**Input:** NCBI accession number (e.g. `NM_001301717`)

**Output:** `output/cleaned_sequence.fasta`

Example output:
```
[INFO] Querying NCBI for accession: 'NM_001301717' ...
[INFO] Fetched 'NM_001301717.2' — 2191 bp
[VALIDATION] No invalid characters detected.
[VALIDATION] Sequence is valid — 2191 bp ready for analysis.
[INFO] Cleaned FASTA written to: output/cleaned_sequence.fasta
```

---

### 2. ORF_finder.py  
Detects all ORFs across all 6 reading frames of the DNA sequence (3 forward + 3
reverse complement). An ORF begins at a start codon (ATG) and ends at the next
in-frame stop codon (TAA, TAG, or TGA).

**Input:** Cleaned DNA sequence string from `input_validate.py`

**Output:** Terminal summary + `output/orfs.csv`

Example terminal output:
```
[ORCA] Processing sequence 1: NM_001301717
[VALIDATION] Sequence is valid — 2191 bp ready for analysis.

----------  ORF Summary — Sequence 1 ----------
  Total ORFs found            : 9
  Forward strand (+)          : 3
  Reverse strand (-)          : 6
  Canonical   (ATG)           : 9
  Nested ORFs detected        : 0
----------------------------------------------
```

Example CSV output (`output/orfs.csv`):
```
NM_001301717
orf_id  strand  start_codon  frame  start  end   length_nt  sequence (5'->3')
ORF1    +       ATG          0      99     1218  1119       ATGAAAAGCGTGCTGGTGGTG...
ORF2    +       ATG          0      1287   1476  189        ATGACTCAGGACATCCCCCCG...
ORF3    -       ATG          1      1842   1935  93         ATGTCATCCCCACTCTGGAGC...
ORF4    +       ATG          2      1583   2183  600        ATGAACCTTCTGGCCTCCCAC...
ORF5    -       ATG          2      1013   1178  165        ATGGAGGAGCGCCGGATGTGC...
ORF6    -       ATG          2      914    998   84         ATGTTGAGTTGCTTACTGAGC...
ORF7    -       ATG          2      806    899   93         ATGAAGACCACGACCACAGCG...
ORF8    -       ATG          2      572    797   225        ATGGCCAGCAGGGGGACCAGA...
ORF9    -       ATG          2      200    284   84         ATGATGGAGTACATGATAGGG...
```

CSV column descriptions:

| Column           | Description                                               |
|------------------|-----------------------------------------------------------|
| `orf_id`         | Unique identifier for each ORF (ORF1, ORF2, ...)         |
| `strand`         | `+` for forward strand, `-` for reverse complement strand |
| `start_codon`    | Start codon detected (always ATG for canonical ORFs)      |
| `frame`          | Reading frame offset (0, 1, or 2)                         |
| `start`          | Start position in the original sequence (0-based)         |
| `end`            | End position in the original sequence (0-based)           |
| `length_nt`      | Length of the ORF in nucleotides (base pairs)             |
| `sequence (5'->3')` | ORF nucleotide sequence read in the 5' to 3' direction |

---

### 3. ORF_analysis.py
Identifies repeated ORFs and calculates pairwise similarity scores between all
detected ORFs. Repeated ORFs are sequences that appear more than once in the genome,
which may indicate gene duplication or conserved functional regions.

**Input:** `output/orfs.csv` from `ORF_finder.py`

**Output:** `output/orf_analysis.csv` + `output/plots/similarity_heatmap.png`

Results include:
- Which ORFs are repeated (identical or near-identical sequences)
- Pairwise similarity scores between every pair of ORFs (0.0 to 1.0)
- A heatmap visualization of the similarity matrix

---

### 4. statistics_summary.py
Generates a comprehensive statistical summary of all ORFs found and produces
publication-quality plots for visualization.

**Input:** `output/orfs.csv` and `output/orf_analysis.csv`

**Output:** `output/statistics_summary.txt` + plots in `output/plots/`

Summary statistics include:
- Total number of ORFs detected
- Longest ORF lengths
- codon useage throughout the sequences 
- ORF counts per reading frame
- Forward vs reverse strand distribution

Plots generated:
- `orf_length_dist.png` — histogram of ORF length distribution
- `strand_distribution.png` — bar chart of forward vs reverse strand ORF counts
- `similarity_heatmap.png` — heatmap of pairwise ORF similarity scores

---

## Example Full Run

```bash
# Step 1 — activate environment
conda activate orf_project_env

# Step 2 — run the full pipeline
python -m src.main --accession NM_001301717 --email your_email@example.com --min-length 75 --ignore-nested

# Step 3 — check outputs
ls output/
```

Expected files in `output/` after a successful run:
```
cleaned_sequence.fasta
orfs.csv
orf_analysis.csv
orf_summary.txt
plots/
    orf_map.png
    codon_comparison.png
    similarity_heatmap.png
```

---

## Notes

- Always activate `ORCA` before running any script.
- The NCBI Entrez API requires a valid email — pass it with `--email`.
- Accession numbers can be found at: https://www.ncbi.nlm.nih.gov/nucleotide/
- The `output/` folder is created automatically on the first run.
- Use `--min-length` to filter out very short ORFs that may not be biologically meaningful.
- Use `--ignore-nested` to exclude ORFs that are contained within longer ORFs.

## Installation

### Dependency Requirements
- Python 3.10
- numpy
- Biopython

### Setup

1. Clone the repository:
```bash
git clone https://github.com/TahmidA139/ORCA.git
```

2. Go into the project folder:
```bash
cd ORCA
```

3. Create the environment:
```bash
conda env create -f environment.yml
```

4. Activate the environment:
```bash
conda activate ORCA
```

## Usage

### Command-Line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--accession` | Yes | — | NCBI accession number for sequence 1 |
| `--accession2` | No | — | NCBI accession number for sequence 2 (enables comparative mode) |
| `--email` | Yes | — | Your email address (required by NCBI Entrez) |
| `--start-codons` | No | `ATG` | One or more start codons: `ATG`, `GTG`, `TTG` |
| `--min-length` | No | `30` | Minimum ORF length in nucleotides |
| `--ignore-nested` | No | `False` | Exclude ORFs nested inside another ORF in the same frame |
| `--output` | No | `output/orfs.csv` | Path for the primary ORF output CSV |

### Usage Examples

Single sequence — all defaults (ATG only, min 30 nt, nested ORFs included):
```bash
python -m src.main --accession NM_001301717 --email you@example.com
```

Single sequence — all three start codons, minimum 150 nt, no nested ORFs:
```bash
python -m src.main --accession NM_001301717 --email you@example.com \
    --start-codons ATG GTG TTG --min-length 150 --ignore-nested
```

Comparative mode — two accessions side by side:
```bash
python -m src.main --accession NM_001301717 --accession2 NM_001301718.2 --email you@example.com
```

## References
- NEEDDDD TO ADD STUFF HERE AS WELL
- NCBI Entrez API: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- Biopython: https://biopython.org/
- NumPy: https://numpy.org/

## License
This project is licensed under the GNU GPL v2.1. Chosen for open collaboration, ease of edits, and public use.

## Authors

**Amanda Yaworsky**
- Student ID: 801489950
- Email: ayaworsk@charlotte.edu

**Erin Nicole Decocker**
- Student ID: 801442694
- Email: edecocke@charlotte.edu

**Tahmid Anwar**
- Student ID: 801501080
- Email: tanwar@charlotte.edu
