#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Tahmid Anwar's Part

## Title: Input Validation Script
## Input:
##   - accession  (str): NCBI nucleotide accession number (sequence 1).
##   - accession2 (str): NCBI nucleotide accession number (sequence 2, optional).
##   - email      (str): Valid email address required by NCBI Entrez API.
## Output:
##   - Raw DNA sequence string(s) retrieved from NCBI, or None if the fetch fails.
## How it works:
##   Accession number(s) are inputted into the function. The function retrieves the
##   FASTA file(s) from NCBI using the accession number(s), removes the header, and
##   returns the raw sequence string(s) as output, which will be used by ORF_finder.py.
##   If two accession numbers are provided, both sequences are fetched, validated,
##   and returned for comparative analysis.

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
import sys

VALID_START_CODONS = {"ATG", "GTG", "TTG"}


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
        print(f"[INFO] Querying NCBI for accession: '{accession}' ...")

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
        print(f"[INFO] Fetched '{record.id}' — {len(sequence)} bp")
        return sequence

    except Exception as e:
        # 'Exception as e' catches any type of error and stores its message in 'e'.
        # We print the error message so the user knows what went wrong,
        # then return None to signal that no sequence was retrieved.
        print(f"[ERROR] Could not fetch sequence from NCBI: {e}")
        return None


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
        print(f"[VALIDATION] Invalid characters found and removed: {invalid_chars}")
        sequence = re.sub(r"[^ATGCRYSWKMBDHVN]", "", sequence)
    else:
        print("[VALIDATION] No invalid characters detected.")

    # Step 4 — verify the cleaned sequence is long enough for ORF analysis
    if len(sequence) < 6:
        print(f"[VALIDATION] Cleaned sequence too short ({len(sequence)} bp). "
              "Minimum required: 6 bp.")
        return False, sequence

    print(f"[VALIDATION] Sequence is valid — {len(sequence)} bp ready for analysis.")
    return True, sequence


# ─────────────────────────────────────────────────────────────────────────────
# HELPER — Write cleaned sequence to a FASTA file
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
    print(f"[INFO] Cleaned FASTA written to: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 3 — Validate start codons
# ─────────────────────────────────────────────────────────────────────────────

def validate_start_codons(requested: list) -> list:
    """
    Upper-case and validate the user-supplied start codons.
    Exits with a helpful message if any unrecognised codon is given.

    Input:
        requested (list): List of start codon strings e.g. ['ATG', 'GTG']
    Output:
        list: Uppercased, validated list of start codons.
    """
    upper = [c.upper() for c in requested]

    # Warn if any non-canonical start codons are requested
    noncanonical_requested = [c for c in upper if c in {"GTG", "TTG"}]
    if noncanonical_requested:
        print(
            f"[WARNING] Non-canonical start codon(s) detected: {', '.join(noncanonical_requested)}\n"
            f"          GTG and TTG produce false positives in eukaryotic sequences.\n"
            f"          Consider using ATG only unless working with prokaryotic sequences."
        )

    unknown = [c for c in upper if c not in VALID_START_CODONS]
    if unknown:
        print(
            f"[ERROR] Unrecognised start codon(s): {', '.join(unknown)}\n"
            f"        Allowed values are: {', '.join(sorted(VALID_START_CODONS))}"
        )
        sys.exit(1)
    return upper


# ─────────────────────────────────────────────────────────────────────────────
# HELPER — Fetch and validate a single sequence
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_and_validate_one(
    accession: str,
    output_fasta: str,
) -> tuple[str, str] | tuple[None, None]:
    ## Input:
    ##   - accession   (str): NCBI accession number.
    ##   - output_fasta(str): Path for the cleaned FASTA output file.
    ## Output:
    ##   - (accession, clean_seq) on success, (None, None) on failure.
    ## Note:
    ##   This is an internal helper used by run() to avoid repeating
    ##   the same fetch-validate-write steps for each sequence.

    # Step 1 — fetch raw sequence from NCBI
    raw_sequence = fetch_fasta_from_ncbi(accession)
    if raw_sequence is None:
        return None, None

    # Step 2 — validate and clean the sequence
    is_valid, clean_seq = validate_dna_sequence(raw_sequence)
    if not is_valid:
        print(f"[ERROR] Sequence {accession} failed validation. Aborting.")
        return None, None

    # Step 3 — write cleaned FASTA to disk
    write_cleaned_fasta(clean_seq, accession, output_fasta)

    return accession, clean_seq


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE — called by main.py
# ─────────────────────────────────────────────────────────────────────────────

def run(
    accession:    str,
    email:        str,
    accession2:   str | None = None,
    output_fasta: str = "output/cleaned_sequence.fasta",
    output_fasta2: str = "output/cleaned_sequence_2.fasta",
) -> tuple:
    ## Input:
    ##   - accession    (str):       NCBI accession number for sequence 1.
    ##   - email        (str):       User email required by NCBI Entrez.
    ##   - accession2   (str|None):  NCBI accession number for sequence 2 (optional).
    ##                               If provided, comparative mode is enabled.
    ##   - output_fasta (str):       Output path for cleaned sequence 1 FASTA.
    ##   - output_fasta2(str):       Output path for cleaned sequence 2 FASTA.
    ## Output:
    ##   Single mode:      (accession,  clean_seq,  None,        None)
    ##   Comparative mode: (accession,  clean_seq,  accession2,  clean_seq2)
    ##   On failure:       (None, None, None, None)

    # ── Step 0 — validate the email before making any NCBI requests ──────────
    # This is checked first so the user gets an immediate clear error message
    # if the email is wrong, rather than a confusing NCBI error later.
    if not validate_email(email):
        # validate_email() already printed the specific error message
        return None, None, None, None

    # Email is valid — set it for all Entrez queries
    # NCBI requires this to identify who is making the request
    Entrez.email = email.strip()
    print(f"[INFO] NCBI email set to: {Entrez.email}")

    # ── Step 1 — fetch and validate sequence 1 ───────────────────────────────
    print(f"\n[INFO] Processing sequence 1: {accession}")
    acc1, seq1 = _fetch_and_validate_one(accession, output_fasta)
    if acc1 is None:
        return None, None, None, None

    # ── Step 2 — fetch and validate sequence 2 (if provided) ─────────────────
    # accession2 is only provided when the user passes --accession2 in main.py
    # If not provided, we run in single sequence mode
    if accession2:
        print(f"\n[INFO] Processing sequence 2: {accession2}")
        acc2, seq2 = _fetch_and_validate_one(accession2, output_fasta2)
        if acc2 is None:
            # Sequence 1 succeeded but sequence 2 failed
            # Return sequence 1 results and None for sequence 2
            print("[WARNING] Sequence 2 failed. Continuing with sequence 1 only.")
            return acc1, seq1, None, None
        print(f"\n[INFO] Both sequences ready for comparative analysis.")
        return acc1, seq1, acc2, seq2

    # Single sequence mode — return None for the second sequence slots
    return acc1, seq1, None, None


# ─────────────────────────────────────────────────────────────────────────────
# Run interactively when the module is run directly
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    print("=" * 60)
    print("  ORCA — Input Validation")
    print("=" * 60)
    print("  Fetch and validate DNA sequences from NCBI.")
    print("  Find accession numbers at:")
    print("  https://www.ncbi.nlm.nih.gov/nucleotide/")
    print("=" * 60)

    # ── Step 1 — get email from user ──────────────────────────────────────────
    # Keep asking until a valid email is entered
    # The loop only exits when validate_email() returns True
    print()
    while True:
        email = input("Enter your email address (required by NCBI): ").strip()
        if validate_email(email):
            # Email is valid — break out of the loop and continue
            break
        # If invalid, validate_email() already printed the error message
        # The loop repeats and asks again

    # ── Step 2 — get accession number 1 from user ─────────────────────────────
    # Keep asking until the user types something (not blank)
    print()
    while True:
        accession = input("Enter accession number for sequence 1 "
                          "(e.g. NM_001301717): ").strip()
        if accession:
            break
        print("[ERROR] Accession number cannot be empty. Please try again.")

    # ── Step 3 — ask if user wants a second sequence (comparative mode) ────────
    print()
    second = input("Do you want to enter a second accession number "
                   "for comparative mode? (yes/no): ").strip().lower()

    # Accept 'yes', 'y' as positive answers — anything else means no
    accession2 = None
    if second in ("yes", "y"):
        while True:
            accession2 = input("Enter accession number for sequence 2 "
                               "(e.g. NM_001838.4): ").strip()
            if accession2:
                break
            print("[ERROR] Accession number cannot be empty. Please try again.")

    # ── Step 4 — run the pipeline with the user-provided inputs ───────────────
    print()
    acc1, seq1, acc2, seq2 = run(
        accession=accession,
        email=email,
        accession2=accession2,
    )

    # ── Step 5 — display results ───────────────────────────────────────────────
    print()
    print("=" * 60)
    print("  RESULTS")
    print("=" * 60)

    if seq1:
        print(f"  Sequence 1 : {acc1}")
        print(f"  Length     : {len(seq1)} bp")
        print(f"  First 60bp : {seq1[:60]}")
    else:
        print("  Sequence 1 : FAILED — check accession number and internet connection")

    if seq2:
        print()
        print(f"  Sequence 2 : {acc2}")
        print(f"  Length     : {len(seq2)} bp")
        print(f"  First 60bp : {seq2[:60]}")
    elif accession2:
        print()
        print("  Sequence 2 : FAILED — check accession number and internet connection")

    print("=" * 60)



