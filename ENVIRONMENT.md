# Setting up the environment
Nicole Decocker's part

## About `environment.yml`

The file `environment.yml` was created by running:

```bash
conda activate ORCA
conda env export --from-history > environment.yml
```

The `--from-history` flag tells conda to record only the packages that were explicitly installed (rather than every sub-dependency resolved during installation). This makes the file portable across operating systems: conda will resolve the correct platform-specific dependencies when the environment is recreated on a different machine.

The environment was originally created with:

```bash
conda create -n ORCA python=3.10 matplotlib numpy -y
```

## Creating the environment from `environment.yml`

From the project's root directory, run:

```bash
conda env create -f environment.yml
conda activate ORCA
```

This will install Python 3.10 along with matplotlib and numpy. The standard library modules used by the project (argparse, csv, os, sys) do not require separate installation.

## Running the program

With the environment activated, run the pipeline from the project's root directory.

**Single-sequence mode** (NCBI accession):

```bash
python -m src.main \
    --accession <ACCESSION_NUMBER> \
    --email     <YOUR_EMAIL>
```

**Single-sequence mode** (local FASTA file):

```bash
python -m src.main \
    --fasta <PATH_TO_FASTA> \
    --email <YOUR_EMAIL>
```

**Comparative mode** (two local FASTA files):

```bash
python -m src.main \
    --fasta  <PATH_TO_FASTA_1> \
    --fasta2 <PATH_TO_FASTA_2> \
    --email  <YOUR_EMAIL>
```

**Comparative mode** (two NCBI accessions):
```bash
python -m src.main \
    --accession  <ACCESSION_NUMBER_1> \
    --accession2 <ACCESSION_NUMBER_2> \
    --email      <YOUR_EMAIL>
```

**Comparative mode** (one NCBI accession + one local FASTA file):
```bash
python -m src.main \
    --accession <ACCESSION_NUMBER> \
    --fasta2    <PATH_TO_FASTA_2> \
    --email     <YOUR_EMAIL>
```
Replace each placeholder with your own values:

| Placeholder | Description |
|-------------|-------------|
| `<ACCESSION_NUMBER>` | NCBI nucleotide accession number, e.g. `NM_012367.1` |
| `<PATH_TO_FASTA>` | Path to a single-sequence FASTA file on your machine, e.g. `example_input_files/OR2B6_sequence.fasta` |
| `<PATH_TO_FASTA_1>` | Path to the first FASTA file (comparative mode) |
| `<PATH_TO_FASTA_2>` | Path to the second FASTA file (comparative mode) |
| `<YOUR_EMAIL>` | A valid email address, e.g. `you@example.com` (required by NCBI Entrez) |
| `<N>` | A number that will vary depending on your sequence and settings |

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--accession` | See note | NCBI accession number for sequence 1. Cannot be used with `--fasta`. |
| `--fasta` | See note | Path to a local FASTA file for sequence 1. Must contain exactly one sequence. Cannot be used with `--accession`. |
| `--accession2` | No | NCBI accession number for sequence 2. Enables comparative mode. |
| `--fasta2` | No | Path to a local FASTA file for sequence 2. Enables comparative mode. Cannot be used with `--accession2`. |
| `--email` | Yes | Email address required by NCBI Entrez. Not used when loading from local files, but must still be supplied to avoid an interactive prompt. |
| `--min-length` | No | Minimum ORF length in nucleotides (default: 30). Must be at least 3. |
| `--start-codons` | No | One or more start codons to search for (default: `ATG`). Non-canonical alternatives: `GTG`, `TTG`. |
| `--outdir` | No | Directory for all output files (default: `output/`). Created automatically if it does not exist. |

**Note:** exactly one of `--accession` or `--fasta` must be provided for sequence 1.

## Expected output

All output files are written to the `output/` directory by default. A different directory can be specified with `--outdir`.

### Single-sequence mode

| File | Description |
|------|-------------|
| `output/cleaned_sequence_1.fasta` | The input sequence after validation and cleaning. |
| `output/<accession>.gff3` | ORF annotations in GFF3 format (filename derived from the accession, e.g. `NM_012367.1.gff3`). |
| `output/orf_map.png` | ORF map showing all ORFs across all six reading frames. |
| `output/orf_summary.txt` | Human-readable summary of ORF statistics. |

### Comparative mode (additional files)

| File | Description |
|------|-------------|
| `output/comp_cleaned_sequence_1.fasta` | Cleaned sequence 1. |
| `output/comp_cleaned_sequence_2.fasta` | Cleaned sequence 2. |
| `output/<accession1>.gff3` | ORF annotations for sequence 1 in GFF3 format. |
| `output/<accession2>.gff3` | ORF annotations for sequence 2 in GFF3 format. |
| `output/orf_map.png` | Side-by-side comparative ORF map for both sequences. |
| `output/orf_comparison_report.txt` | Human-readable side-by-side ORF statistics and codon usage for both sequences. |
| `output/codon_usage_comparison.png` | RSCU codon-usage heatmap comparing the two sequences. |

When run with a single FASTA file using default settings, the terminal output should look like this:

```
[ORCA] Processing sequence 1: (local file) <PATH_TO_FASTA>
[VALIDATION] Sequence is valid — <N> bp ready for analysis.

════════════════════════ ORF Summary — Sequence 1 ════════════════════════
  Total ORFs found  : <N>
  Forward strand (+): <N>
  Reverse strand (-): <N>
  Canonical (ATG)   : <N>
────────────────────────────────────────────────────────────────────────
```
