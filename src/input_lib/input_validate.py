#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Tahmid Anwar's Part

## Title: Input Validation Script
## Input:
##   - accession  (str): NCBI nucleotide accession number (sequence 1).
##   - accession2 (str): NCBI nucleotide accession number (sequence 2, optional).
##   - email      (str): Valid email address required by NCBI Entrez API.
##   - fasta_file (str): Path to a local single-sequence FASTA file (alternative to accession).
##   - fasta_file2(str): Path to a local single-sequence FASTA file for sequence 2 (optional).
## Output:
##   - Raw DNA sequence string(s) retrieved from NCBI or loaded from disk,
##     or None if the fetch/load fails.
## How it works:
##   Accession number(s) are inputted into the function. The function retrieves the
##   FASTA file(s) from NCBI using the accession number(s), removes the header, and
##   returns the raw sequence string(s) as output, which will be used by ORF_finder.py.
##   If a local FASTA file path is provided instead, the sequence is loaded from disk.
##   If two sequences are provided, both sequences are fetched/loaded, validated,
##   and returned for comparative analysis.

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
import sys

VALID_START_CODONS = {"ATG"}


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 0 — Validate email address
# ─────────────────────────────────────────────────────────────────────────────

def validate_email(email: str) -> bool:
    ## Input:
    ##   - email (str): Email address string provided by the user.
    ## Output:
    ##   - bool: True if the email looks valid, False otherwise.
    ## How it works:
    ##   Uses a regex pattern to check the email has the correct format:
    ##   something @ something . something
    ##   e.g. tahmid@gmail.com or student@charlotte.edu
    ##   This does NOT check if the email actually exists — only that the
    ##   format is correct, which is all NCBI requires.

    # Regex pattern breakdown:
    # [^@\s]+   — one or more characters that are not @ or whitespace (the username)
    # @         — the @ symbol
    # [^@\s]+   — one or more characters that are not @ or whitespace (the domain)
    # \.        — a literal dot
    # [^@\s]+   — one or more characters after the dot (e.g. com, edu, org)
    pattern = r"^[^@\s]+@[^@\s]+\.[^@\s]+$"

    if not email or not email.strip():
        print("[ERROR] Invalid email: email cannot be empty.")
        return False

    if not re.match(pattern, email.strip()):
        print(f"[ERROR] Invalid email: '{email}' is not a valid email address.")
        print("        Expected format: username@domain.com")
        print("        Example:         tahmid@gmail.com")
        return False

    return True


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 1 — Fetch FASTA from NCBI
# ─────────────────────────────────────────────────────────────────────────────

def fetch_fasta_from_ncbi(accession: str, db: str = "nucleotide") -> str | None:
    ## Input:
    ##   - accession (str): NCBI nucleotide accession number
    ##                      e.g. 'NM_001301717' or 'NC_000913.3'
    ##   - db (str):        Entrez database to query (default: 'nucleotide')
    ## Output:
    ##   - sequence (str):  Raw DNA sequence string retrieved from NCBI,
    ##                      or None if the fetch fails.

    # try/except is used here because many things can go wrong such as bad
    # accession number, no internet, NCBI is down, rate limit hit, etc.
    # Instead of crashing the whole program, any error is caught and None
    # is returned so the caller can handle the failure gracefully.
    try:
        # Step 1 & 2 — query NCBI and parse the FASTA response
        # Entrez.efetch is the actual API call to NCBI:
        # db: which NCBI database to search (e.g. 'nucleotide')
        # id: the accession number to look up
        # rettype: the file format to return ('fasta')
        # retmode: how the data is encoded ('text' = plain readable text)
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")

        # 'handle' is a file-like object — SeqIO reads and parses it,
        # automatically separating the header from the sequence
        record = SeqIO.read(handle, "fasta")

        # Always close the handle — leaving it open can cause memory
        # leaks or hit NCBI's rate limits
        handle.close()

        # Step 3 — extract the raw sequence string
        # record.seq is a Biopython Seq object (not a plain Python string).
        # We wrap it in str() to convert it to a regular string, which
        # takes only the sequence (header already removed by SeqIO)
        sequence = str(record.seq)
        return sequence

    except Exception as e:
        # 'Exception as e' catches any type of error and stores its message in 'e'.
        # We print the error message so the user knows what went wrong,
        # then return None to signal that no sequence was retrieved.
        print(f"[ERROR] Could not fetch sequence from NCBI: {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 1b — Load a single sequence from a local FASTA file
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta_from_file(filepath: str) -> tuple[str, str] | tuple[None, None]:
    ## Input:
    ##   - filepath (str): Path to a local FASTA file on disk.
    ## Output:
    ##   - (record_id, sequence) on success, (None, None) on failure.
    ## How it works:
    ##   Opens the file, counts how many sequences are inside it, and
    ##   immediately exits with an error if more than one sequence is found.
    ##   ORCA is designed to work with one sequence per run; multi-sequence
    ##   FASTA files are ambiguous and not supported as input.

    # Step 1 — make sure the file actually exists before trying to open it
    if not os.path.isfile(filepath):
        print(f"[ERROR] FASTA file not found: '{filepath}'")
        print(f"        Please check the path and try again.")
        return None, None

    # Step 2 — parse all records from the file
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        print(f"[ERROR] Could not read FASTA file '{filepath}': {e}")
        return None, None

    # Step 3 — reject empty files
    if len(records) == 0:
        print(f"[ERROR] No sequences found in '{filepath}'.")
        print(f"        The file may be empty or not in valid FASTA format.")
        return None, None

    # Step 4 — HARD STOP: multi-sequence FASTA files are not allowed.
    # A multi-sequence FASTA contains more than one '>' header line.
    # ORCA processes exactly one sequence at a time. If the user supplies
    # a file with multiple sequences it is unclear which one to use, so
    # the pipeline exits immediately with a clear explanation.
    if len(records) > 1:
        print(
            f"\n[ERROR] '{filepath}' contains {len(records)} sequences.\n"
            f"        ORCA requires a single-sequence FASTA file.\n"
            f"\n"
            f"        A valid single-sequence FASTA looks like this:\n"
            f"          >NM_001301717.1 Homo sapiens ...\n"
            f"          ATGCGATCGATCGATCG...\n"
            f"\n"
            f"        Your file has {len(records)} entries starting with:\n"
        )
        # Show the IDs of the first few sequences so the user knows what's in it
        for rec in records[:5]:
            print(f"          > {rec.id}")
        if len(records) > 5:
            print(f"          ... and {len(records) - 5} more.")
        print(
            f"\n"
            f"        Please extract a single sequence into its own FASTA file\n"
            f"        and re-run ORCA with that file.\n"
        )
        sys.exit(1)

    # Step 5 — single record confirmed, return it
    record = records[0]
    sequence = str(record.seq)
    return record.id, sequence


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 2 — Validate and clean the DNA sequence
# ─────────────────────────────────────────────────────────────────────────────

def validate_dna_sequence(sequence: str) -> tuple[bool, str]:
    ## Input:
    ##   - sequence (str): DNA sequence string (any case, may contain noise).
    ## Output:
    ##   - is_valid  (bool): True if the sequence is usable after cleaning.
    ##   - clean_seq (str):  Uppercase DNA string with invalid characters removed.

    # Step 1 — normalise: uppercase, remove whitespace and digits
    sequence = re.sub(r"[\s\d]", "", sequence).upper()

    if not sequence:
        print("[VALIDATION] Sequence is empty after initial cleaning.")
        return False, ""

    # Step 2 — identify characters outside the IUPAC nucleotide alphabet
    # Accepted bases: A T G C + ambiguity codes R Y S W K M B D H V N
    invalid_chars = set(re.findall(r"[^ATGCRYSWKMBDHVN]", sequence))

    # Step 3 — flag and remove invalid characters
    if invalid_chars:
        sequence = re.sub(r"[^ATGCRYSWKMBDHVN]", "", sequence)

    # Step 4 — check for IUPAC ambiguity codes (valid but not plain A/T/G/C)
    ambiguous_chars = set(re.findall(r"[^ATGC]", sequence))
    if ambiguous_chars:
        print(f"[VALIDATION] IUPAC ambiguity codes detected ({', '.join(sorted(ambiguous_chars))}), sequence cleaned.")

    # Step 5 — verify the cleaned sequence is long enough for ORF analysis
    if len(sequence) < 6:
        print(f"[VALIDATION] Cleaned sequence too short ({len(sequence)} bp). "
              "Minimum required: 6 bp.")
        return False, sequence

    print(f"[VALIDATION] Sequence is valid — {len(sequence)} bp ready for analysis.")
    return True, sequence


# ─────────────────────────────────────────────────────────────────────────────
# HELPER — Write a single cleaned sequence to a FASTA file
# ─────────────────────────────────────────────────────────────────────────────

def write_cleaned_fasta(
    clean_seq: str,
    accession: str,
    output_path: str = "output/cleaned_sequence.fasta",
) -> None:
    ## Input:
    ##   - clean_seq   (str): Validated, cleaned DNA string.
    ##   - accession   (str): Accession number used as the FASTA record ID.
    ##   - output_path (str): Destination file path.
    ## Output:
    ##   - Writes a cleaned FASTA file to disk.

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    record = SeqRecord(
        Seq(clean_seq),
        id=accession,
        description="cleaned sequence | ready for ORF analysis",
    )
    with open(output_path, "w") as fh:
        SeqIO.write(record, fh, "fasta")

# ─────────────────────────────────────────────────────────────────────────────
# HELPER — Write both cleaned sequences into one combined FASTA file
# ─────────────────────────────────────────────────────────────────────────────

def write_combined_cleaned_fasta(
    clean_seq1: str,
    accession1: str,
    clean_seq2: str,
    accession2: str,
    output_path: str = "output/cleaned_sequences.fasta",
) -> None:
    ## Input:
    ##   - clean_seq1  (str): Validated, cleaned DNA string for sequence 1.
    ##   - accession1  (str): Accession / ID label for sequence 1.
    ##   - clean_seq2  (str): Validated, cleaned DNA string for sequence 2.
    ##   - accession2  (str): Accession / ID label for sequence 2.
    ##   - output_path (str): Destination file path.
    ## Output:
    ##   - Writes both cleaned sequences into a single two-record FASTA file.
    ## Why:
    ##   In comparative mode the user may want a single file containing both
    ##   cleaned sequences for downstream tools.  This file uses the name
    ##   'cleaned_sequences.fasta' (plural) to distinguish it from the
    ##   single-sequence 'cleaned_sequence.fasta' produced in single mode.

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    records = [
        SeqRecord(
            Seq(clean_seq1),
            id=accession1,
            description="cleaned sequence 1 | ready for ORF analysis",
        ),
        SeqRecord(
            Seq(clean_seq2),
            id=accession2,
            description="cleaned sequence 2 | ready for ORF analysis",
        ),
    ]
    with open(output_path, "w") as fh:
        SeqIO.write(records, fh, "fasta")


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 3 — Validate start codons
# ─────────────────────────────────────────────────────────────────────────────

def validate_start_codons(requested: list) -> list:
    """
    Upper-case and validate the user-supplied start codons.
    Exits with a helpful message if any unrecognised codon is given.

    Input:
        requested (list): List of start codon strings e.g. ['ATG']
    Output:
        list: Uppercased, validated list of start codons.
    """
    upper = [c.upper() for c in requested]

    unknown = [c for c in upper if c not in VALID_START_CODONS]
    if unknown:
        print(
            f"[ERROR] Unrecognised start codon(s): {', '.join(unknown)}\n"
            f"        Allowed value is: ATG"
        )
        sys.exit(1)
    return upper


# ─────────────────────────────────────────────────────────────────────────────
# HELPER — Fetch/load and validate a single sequence
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_and_validate_one(
    accession: str,
    output_fasta: str,
    fasta_file: str | None = None,
) -> tuple[str, str] | tuple[None, None]:
    ## Input:
    ##   - accession   (str):       NCBI accession number (used when fasta_file is None).
    ##   - output_fasta(str):       Path for the cleaned FASTA output file.
    ##   - fasta_file  (str|None):  Path to a local FASTA file.
    ##                              When provided, the local file is used instead of NCBI.
    ## Output:
    ##   - (accession, clean_seq) on success, (None, None) on failure.
    ## Note:
    ##   This is an internal helper used by run() to avoid repeating
    ##   the same fetch-validate-write steps for each sequence.

    if fasta_file is not None:
        # ── Local file path: load from disk ───────────────────────────────────
        # load_fasta_from_file() will sys.exit(1) if the file has multiple
        # sequences, so if we reach the next line the file is guaranteed to
        # contain exactly one record.
        record_id, raw_sequence = load_fasta_from_file(fasta_file)
        if raw_sequence is None:
            return None, None
        # Use the sequence ID embedded in the FASTA header as the label,
        # overriding whatever placeholder accession string was passed in.
        accession = record_id
    else:
        # ── NCBI path: fetch by accession number ──────────────────────────────
        raw_sequence = fetch_fasta_from_ncbi(accession)
        if raw_sequence is None:
            return None, None

    # Step 2 — validate and clean the sequence (shared for both paths)
    is_valid, clean_seq = validate_dna_sequence(raw_sequence)
    if not is_valid:
        print(f"[ERROR] Sequence '{accession}' failed validation. Aborting.")
        return None, None

    # Step 3 — write the individual cleaned FASTA to disk
    write_cleaned_fasta(clean_seq, accession, output_fasta)

    return accession, clean_seq


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE — called by main.py
# ─────────────────────────────────────────────────────────────────────────────

def run(
    accession:     str,
    email:         str,
    accession2:    str | None = None,
    output_fasta:  str = "output/cleaned_sequence.fasta",
    output_fasta2: str = "output/cleaned_sequence_2.fasta",
    fasta_file:    str | None = None,
    fasta_file2:   str | None = None,
) -> tuple:
    ## Input:
    ##   - accession    (str):       NCBI accession number for sequence 1.
    ##                               Ignored when fasta_file is provided.
    ##   - email        (str):       User email required by NCBI Entrez.
    ##                               Still validated even in local-file mode.
    ##   - accession2   (str|None):  NCBI accession number for sequence 2 (optional).
    ##   - output_fasta (str):       Output path for cleaned sequence 1 FASTA.
    ##   - output_fasta2(str):       Output path for cleaned sequence 2 FASTA.
    ##   - fasta_file   (str|None):  Local FASTA file for sequence 1 (overrides NCBI).
    ##   - fasta_file2  (str|None):  Local FASTA file for sequence 2 (overrides NCBI).
    ## Output:
    ##   Single mode:      (accession,  clean_seq,  None,        None)
    ##   Comparative mode: (accession,  clean_seq,  accession2,  clean_seq2)
    ##   On failure:       (None, None, None, None)

    # ── Step 0 — validate the email before making any NCBI requests ──────────
    # Even when using local FASTA files we validate the email format so that
    # if the user switches to NCBI mode later the email is already confirmed.
    if not validate_email(email):
        return None, None, None, None

    # Email is valid — set it for all Entrez queries
    # NCBI requires this to identify who is making the request
    Entrez.email = email.strip()

    # ── Step 1 — fetch/load and validate sequence 1 ──────────────────────────
    acc1, seq1 = _fetch_and_validate_one(accession, output_fasta, fasta_file)
    if acc1 is None:
        return None, None, None, None

    # ── Step 2 — fetch/load and validate sequence 2 (if provided) ────────────
    if accession2 or fasta_file2:
        acc2, seq2 = _fetch_and_validate_one(
            accession2 or "",
            output_fasta2,
            fasta_file2,
        )
        if acc2 is None:
            return acc1, seq1, None, None

        write_combined_cleaned_fasta(
            clean_seq1=seq1, accession1=acc1,
            clean_seq2=seq2, accession2=acc2,
            output_path="output/cleaned_sequences.fasta",
        )

        return acc1, seq1, acc2, seq2

    # Single sequence mode — return None for the second sequence slots
    return acc1, seq1, None, None

