<p align="center">
  <img src="images/orca_banner.png" alt="ORCA Banner" width="350"/>
</p>

# ORCA (ORF Recognition and Comparative Annotator)

## Objective
ORCA is a command-line bioinformatics pipeline that automates ORF detection and analysis in DNA sequences. It helps researchers identify potential protein-coding regions, compare ORF structure across two transcripts and between species, and generates detailed statistics for further study. All from a single NCBI accession number (or two)!

## Features
- Downloads DNA sequences in FASTA format directly from NCBI using an accession number
- Validates and cleans sequences before analysis
- Scans all six reading frames (+1, +2, +3, −1, −2, −3) using NumPy vectorization
- Detects canonical (ATG) and non-canonical (GTG, TTG) start codons
- Separates complete vs incomplete ORFs and identifies nested ORFs
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
│   │   └── orf_finder.py              
│   │
│   ├── analysis_lib/                  # Amanda Yaworsky
│   │   ├── __init__.py
│   │   └── orf_analysis.py            
│   │                                  
│   │
│   └── statistics_lib/                # Erin Nicole Decocker 
│   │   ├── __init__.py
│   |   └── statistics_summary.py      
│   │
│   │
│   ├── graphics_lib/                  # Tahmid Anwar
│   │   ├── __init__.py
│   │   └── graphics.py                                        
│
└── examples/                          # Amanda Yaworsky
    ├── example_output.fasta
    └── example_run.txt
```

### File Descriptions
need to write 

## Output Files

| File | Mode | Description |
|------|------|-------------|
need to Wrute lol

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
python main.py --accession NM_001301717 --email you@example.com
```

Single sequence — all three start codons, minimum 60 nt, no nested ORFs:
```bash
python main.py --accession NM_001301717 --email you@example.com \
    --start-codons ATG GTG TTG --min-length 60 --ignore-nested
```

Comparative mode — two accessions side by side:
```bash
python main.py --accession NM_001301717 --accession2 NM_001256799 \
    --email you@example.com
```

## Algorithm Description

### Input Validation (`input_validate.py`)
need to write 

### ORF Detection (`orf_finder.py`)
need to write 

### ORF Analysis (`orf_analysis.py`)
need to write 

### Statistics Output (`statistics_summary.py`)
need to write 

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
