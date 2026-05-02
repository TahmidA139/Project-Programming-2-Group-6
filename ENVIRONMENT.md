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
    --accession NC_001422 \
    --email     <YOUR_EMAIL>
```

**Single-sequence mode** (local FASTA file):

```bash
python -m src.main \
    --fasta example_input_files/OR2B6_sequence.fasta \
    --email <YOUR_EMAIL>
```

**Comparative mode** (two local FASTA files):

```bash
python -m src.main \
    --fasta  example_input_files/homo_sapiens_albumin.fasta \
    --fasta2 example_input_files/mus_musculus_albumin.fasta \
    --email  <YOUR_EMAIL>
```

**Comparative mode** (two NCBI accessions):
```bash
python -m src.main \
    --accession  NC_001416 \
    --accession2 NC_001604 \
    --email      <YOUR_EMAIL>
```

**Comparative mode** (one NCBI accession + one local FASTA file):
```bash
python -m src.main \
    --accession NM_000477.7 \
    --fasta2    example_input_files/mus_musculus_albumin.fasta \
    --email     <YOUR_EMAIL>
```
Replace each placeholder `<YOUR_EMAIL>` with a valid email address, e.g. `you@example.com` (required by NCBI Entrez) 

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

#### Notes on `--min-length` and `--start-codons` and RSCU (Relative Synonymous Codon Usage)

The default minimum ORF length of **30 nt** (10 amino acids) is a common threshold in prokaryotic and eukaryotic ORF annotation. ORFs shorter than this are statistically likely to appear by chance in any random sequence and are usually not biologically meaningful. Raise this value if you want to focus only on longer, more likely functional coding sequences; lower it (minimum 3, one codon) if you are working with very short sequences or want to capture all possible reading frames.

`ATG` is the universal canonical start codon and the default. Non-canonical start codons (`GTG`, `TTG`) are occasionally used in bacteria and bacteriophages but are rare in eukaryotes. Add them with `--start-codons ATG GTG TTG` only if you have reason to believe your sequence uses them.

RSCU (Relative Synonymous Codon Usage) measures how often each codon is used relative to what you would expect if all codons for the same amino acid were used equally, revealing organism-specific translation biases. For the RSCU heatmap to reflect genuine codon bias rather than sampling noise, the coding sequence should contain at least ~200 codons (~600 bp of ORF sequence); below that threshold, many codons will simply not appear in the data and the heatmap will show apparent bias that is an artefact of small sample size rather than a true biological signal.

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

When run in comparative mode with two sequences, the terminal output should look like this:

```
[ORCA] Processing sequence 1: (local file) <PATH_TO_FASTA_1>
[VALIDATION] Sequence is valid — <N> bp ready for analysis.
[ORCA] Processing sequence 2: (local file) <PATH_TO_FASTA_2>
[VALIDATION] Sequence is valid — <N> bp ready for analysis.

════════════════════════ ORF COMPARISON REPORT ═════════════════════════
  Run Parameters
  ──────────────
  Generated    : <DATETIME>
  Sequence 1   : <ACCESSION_1>
  Sequence 2   : <ACCESSION_2>
  Start codons : ATG
  Min length   : 30 nt

────────────────────────────────────────────────────────────────────────

═══════════════════════ Sequence 1 — <ACCESSION_1> ═══════════════════════

  Dataset Statistics
  ──────────────────
  Total ORFs         : <N>
  Average GC Content : <N>%
  ...

═══════════════════════ Sequence 2 — <ACCESSION_2> ═══════════════════════
  ...
────────────────────────────────────────────────────────────────────────
```

The full comparison report is also written to `output/orf_comparison_report.txt`. All `<N>` values are placeholders that will vary depending on your sequences and settings.

## Cleanup

To deactivate the environment when you are done:

```bash
conda deactivate
```

To remove the environment entirely (e.g. to reinstall from scratch):

```bash
conda env remove -n ORCA
```

## Troubleshooting

Below are all error messages you may encounter, what causes them, and how to fix them.

---

### Email errors

**`[ERROR] Invalid email: email cannot be empty.`**  
The `--email` flag was provided but left blank. Supply a valid email address, e.g. `--email you@example.com`.

**`[ERROR] Invalid email: '<value>' is not a valid email address.`**  
The value passed to `--email` does not match the expected `username@domain.tld` format. Check for typos — common mistakes include missing the `@` sign or the domain suffix.

---

### Input / FASTA errors

**`[ERROR] FASTA file not found: '<path>'`**  
The path passed to `--fasta` or `--fasta2` does not point to an existing file. Check that the path is correct and that you are running the command from the project's root directory. Example: `--fasta example_input_files/OR2B6_sequence.fasta`, not `--fasta OR2B6_sequence.fasta`.

**`[ERROR] Could not read FASTA file '<path>': <details>`**  
The file exists but could not be parsed. The file may be corrupted, empty, or in the wrong format. Ensure it is a plain-text FASTA file (starts with a `>` header line followed by nucleotide sequence).

**`[ERROR] No sequences found in '<path>'.`**  
The file was opened successfully but contains no FASTA records. The file may be empty or missing the `>` header line.

**`[ERROR] '<path>' contains <N> sequences. ORCA requires a single-sequence FASTA file.`**  
The FASTA file contains more than one sequence record. ORCA processes one sequence at a time. Extract the sequence you want into its own file and pass that file instead.

---

### NCBI fetch errors

**`[ERROR] Could not fetch sequence from NCBI: <details>`**  
ORCA was unable to retrieve the sequence for the given accession number. Common causes are: an incorrect or retired accession number, no internet connection, or NCBI being temporarily unavailable. Check the accession against the [NCBI Nucleotide database](https://www.ncbi.nlm.nih.gov/nucleotide/) and try again. NCBI also enforces rate limits on unauthenticated requests — if you are running many queries, wait a moment and retry.

---

### Sequence validation errors

**`[VALIDATION] Sequence is empty after initial cleaning.`**  
After stripping whitespace and digits the sequence contained no characters. The FASTA file likely has a header but no sequence data below it.

**`[VALIDATION] IUPAC ambiguity codes detected (<codes>); replacing with 'N' to preserve reading frame.`**  
This is an informational warning, not a fatal error. The sequence contains IUPAC ambiguity codes (e.g. `R`, `Y`, `S`). These are replaced with `N` so that the reading frame is preserved. Codons containing `N` are skipped during ORF scanning and RSCU calculation. If a large proportion of your sequence is ambiguous, ORF detection results may be incomplete.

**`[VALIDATION] Cleaned sequence too short (<N> bp). Minimum required: 6 bp.`**  
After cleaning, the sequence is shorter than 6 bp, which is too short to contain any ORF (even the smallest possible ORF is one codon + stop = 6 bp). Check that the input file contains a real nucleotide sequence.

**`[ERROR] Sequence '<accession>' failed validation. Aborting.`**  
Validation returned false for the named sequence (due to one of the above cleaning failures). Correct the input sequence and re-run.

---

### ORF-finding errors

**`[ERROR] No ORFs found for sequence 1. Cannot produce output. Try lowering --min-length or adding non-canonical start codons.`**  
No ORFs meeting the current criteria were found in sequence 1. Try reducing `--min-length` (e.g. `--min-length 15`) or adding non-canonical start codons (e.g. `--start-codons ATG GTG TTG`). If the sequence is very short or heavily ambiguous (`N`-rich), it may genuinely contain no detectable ORFs.

**`[ERROR] No ORFs found for sequence 2. Cannot produce comparative output. Try lowering --min-length or adding non-canonical start codons.`**  
Same as above, but for the second sequence in a comparative run.

**`[ERROR] Pipeline failed: could not retrieve a valid sequence.`**  
The input/validation stage returned no usable sequence. A more specific error will have been printed above this line — scroll up to find and address it first.

**`[ERROR] Pipeline failed for sequence 2.`**  
The pipeline failed during processing of the second sequence in a comparative run. A more specific error will appear earlier in the terminal output.

---

### Argument errors

**`[ERROR] --min-length must be at least 3 (one codon).`**  
The value passed to `--min-length` is less than 3. The minimum meaningful ORF is a single codon (3 nucleotides). Use `--min-length 3` or higher.

**`[ERROR] Unrecognised start codon(s): <codon>. Allowed values are: ATG, GTG, TTG`**  
One or more values passed to `--start-codons` are not recognised. Only `ATG`, `GTG`, and `TTG` are supported.
