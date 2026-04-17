<p align="center">
  <img src="images/orca_banner.png" alt="ORCA Banner" width="350"/>
</p>

# ORCA (ORF Recognition and Comparative Annotator)

## Objective
ORCA is a python-based bioinformatics pipeline that automates Open Reading Frame (ORF) detection and analysis in DNA sequences. It helps researchers identify potential protein-coding regions while generating summary statistics and visualizations for further study.

The pipeline accepts an NCBI accession number and email as inputs, fetches the DNA sequence, detects all ORFs across all six reading frames, and infers the most likely genetic code, with comparative analysis available when two sequences are provided.

## Features
- Downloads DNA sequences in FASTA format directly from NCBI using an accession number
- Validates and cleans sequences before analysis
- Scans all six reading frames (+1, +2, +3, −1, −2, −3) using NumPy vectorization
- Detects canonical (ATG) and non-canonical (GTG, TTG) start codons  (non-canotical is not working great at the moment :) )
- Computes per-ORF statistics: GC content, codon usage, and protein length
- infers the most likely genetic code of a sequence (still working on this currently)
- **Comparative mode** (`--accession2`): side-by-side ORF structure comparison, codon usage differences, and other analysis between two sequences
- Writes all results to the `output/` folder as .CSV and plain-text reports

## Project Structure
```
ORCA/
├── README.md                         # Tahmid Anwar
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
│   ├── orf_finder_lib/                # Nicole Decocker
│   │   ├── __init__.py
│   │   └── frame_scanner.py   
│   │   └── orf_finder.py
│   │   └── output_writer.py                
│   │
│   ├── analysis_lib/                  # Amanda Yaworsky
│   │   ├── __init__.py
│   │   └── orf_analysis.py            
│   │   └── statistics_summary.py                                    
│   │   └── genetic_inference.py       # Tahmid Anwar
│   │
│   ├── graphics_lib/                  # Nicole Decocker 
│   │   ├── __init__.py
│   │   └── graphics.py                                        
│
└── examples/                          # Amanda Yaworsky
    ├── example_output.fasta
    └── example_run.txt
```
Flowchart: 
<p align="center">
  <img src="images/flowchart.png" alt="flowchart" width="350"/>
</p>

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
| `--output` | No | `output/orfs.csv` | Path for the primary ORF output CSV |

### Usage Examples

Single sequence — all defaults (ATG only, min 30 nt, nested ORFs included):
```bash
python -m src.main --accession NM_001301717 --email you@example.com
```

Single sequence — all three start codons and minimum 150 nt:
```bash
python -m src.main --accession NM_001301717 --email you@example.com \
    --start-codons ATG GTG TTG --min-length 150 
```

Comparative mode — two accessions side by side:
```bash
python -m src.main --accession NM_001301717 --accession2 NM_001301718.2 --email you@example.com
```

### File Descriptions
### 0. main.py
The entry point of the entire project. It coordinates all modules by calling them
in the correct order and passing outputs from one module as inputs to the next.

### 1. input_validate.py
Fetches a DNA sequence from NCBI using the accession number, validates it, removes
any invalid characters, and writes a cleaned FASTA file to disk.


### 2. ORF_finder.py  
Detects all ORFs across all 6 reading frames of the DNA sequence (3 forward + 3
reverse complement). An ORF begins at a start codon (ATG) and ends at the next
in-frame stop codon (TAA, TAG, or TGA).


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



### 4. statistics_summary.py
Generates a comprehensive statistical summary of all ORFs found and produces
publication-quality plots for visualization.


## Example Full Run


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
```

---

## Notes

- Always activate `ORCA` before running any script.
- The NCBI Entrez API requires a valid email — pass it with `--email`.
- Accession numbers can be found at: https://www.ncbi.nlm.nih.gov/nucleotide/
- The `output/` folder is created automatically on the first run.
- Use `--min-length` to filter out very short ORFs that may not be biologically meaningful.
- Use `--ignore-nested` to exclude ORFs that are contained within longer ORFs.
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
