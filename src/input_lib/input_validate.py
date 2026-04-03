#!/usr/bin/env python3
# -*- coding: utf-8 -*-


## Tahmid Anwar's Part 


## Title: Input Validation Script
## Input:
##   - accession (str): NCBI nucleotide accession number.
## Output:
##   - Raw DNA sequence string retrieved from NCBI, or None if the fetch fails.
## How it works:
##   Accession number is inputted into the function. The function retrieves the
##   FASTA file from NCBI using that accession number, removes the header, and
##   returns the raw sequence string as output, which will be used by ORF_finder.py.

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
import argparse

VALID_START_CODONS = {"ATG", "GTG", "TTG"}

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


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE — called by main.py
# ─────────────────────────────────────────────────────────────────────────────

def run(
    accession: str,
    email: str,
    output_fasta: str = "output/cleaned_sequence.fasta",
) -> tuple[str, str] | tuple[None, None]:
    ## Input:
    ##   - accession    (str): NCBI accession number.
    ##   - email        (str): User email required by NCBI Entrez.
    ##   - output_fasta (str): Path for the cleaned FASTA output file.
    ## Output:
    ##   - (accession, clean_seq): tuple on success, (None, None) on failure.

    # NCBI requires a valid email address for all Entrez queries.
    # This is not a login — NCBI uses it only to contact you if your
    # script sends too many requests.
    Entrez.email = email

    # Step 1 — fetch raw sequence from NCBI
    raw_sequence = fetch_fasta_from_ncbi(accession)
    if raw_sequence is None:
        return None, None

    # Step 2 — validate and clean the sequence
    is_valid, clean_seq = validate_dna_sequence(raw_sequence)
    if not is_valid:
        print("[ERROR] Sequence failed validation. Aborting pipeline.")
        return None, None

    # Step 3 — write cleaned FASTA to disk
    write_cleaned_fasta(clean_seq, accession, output_fasta)

    # Step 4 — return clean sequence for ORF_finder.py
    return accession, clean_seq

def validate_start_codons(requested: list) -> list:
    """
    Upper-case and validate the user-supplied start codons.
    Exits with a helpful message if any unrecognised codon is given.
    """
    upper   = [c.upper() for c in requested]
     
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


