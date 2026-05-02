#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
input_validate.py

Tahmid Anwar's Part

Purpose:
    Input validation and sequence retrieval for the ORCA pipeline.
    Handles fetching DNA sequences from NCBI via accession number or loading
    them from local FASTA files, then validates and cleans them before
    passing them to the ORF finder.

Public API
----------
validate_email          Validate the format of a user-supplied email address.
fetch_fasta_from_ncbi   Fetch a raw DNA sequence string from NCBI Entrez.
load_fasta_from_file    Load a single sequence from a local FASTA file.
validate_dna_sequence   Clean and validate a raw DNA sequence string.
write_cleaned_fasta     Write one cleaned sequence to a FASTA file.
validate_start_codons   Validate and normalise a list of start codon strings.
run                     Main pipeline entry point called by main.py.
"""

from __future__ import annotations

import os
import re
import sys

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VALID_START_CODONS: set[str] = {"ATG", "GTG", "TTG"}

# Validate email address
def validate_email(email: str) -> bool:
    """
    Validate the format of a user-supplied email address.

    Uses a regex pattern to check that the address follows the shape
    ``username@domain.tld``.  This does not verify that the address
    actually exists — only that the format is correct, which is all
    NCBI Entrez requires.

    Parameters
    ----------
    email : str
        Email address string provided by the user.

    Returns
    -------
    bool
        True if the format is valid, False otherwise.
        Prints a descriptive error message on failure.
    """
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


def fetch_fasta_from_ncbi(accession: str, db: str = "nucleotide") -> str | None:
    """
    Fetch a raw DNA sequence string from NCBI Entrez.

    Parameters
    ----------
    accession : str
        NCBI nucleotide accession number, e.g. ``'NM_001301717'`` or
        ``'NC_000913.3'``.
    db : str, optional
        Entrez database to query.  Defaults to ``'nucleotide'``.

    Returns
    -------
    str or None
        Raw DNA sequence string on success, or ``None`` if the fetch
        fails for any reason (bad accession, no internet, NCBI down,
        rate limit exceeded, etc.).
    """
    try:
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)

    except Exception as e:
        print(f"[ERROR] Could not fetch sequence from NCBI: {e}")
        return None


def load_fasta_from_file(filepath: str) -> tuple[str, str] | tuple[None, None]:
    """
    Load a single DNA sequence from a local FASTA file.

    ORCA processes exactly one sequence at a time.  If the file contains
    more than one record the function raises ``ValueError`` with a
    descriptive message rather than silently choosing one or killing the
    process.  Callers (i.e. CLI entry points) are responsible for catching
    the exception and deciding how to handle it.

    Parameters
    ----------
    filepath : str
        Path to a local FASTA file on disk.

    Returns
    -------
    tuple[str, str]
        ``(record_id, sequence)`` on success.
    tuple[None, None]
        On failure (file not found, unreadable, or empty).

    Raises
    ------
    ValueError
        If the file contains more than one sequence record.
    """
    if not os.path.isfile(filepath):
        print(f"[ERROR] FASTA file not found: '{filepath}'")
        print(f"        Please check the path and try again.")
        return None, None

    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        print(f"[ERROR] Could not read FASTA file '{filepath}': {e}")
        return None, None

    if len(records) == 0:
        print(f"[ERROR] No sequences found in '{filepath}'.")
        print(f"        The file may be empty or not in valid FASTA format.")
        return None, None

    if len(records) > 1:
        # Build a readable preview of the first few record IDs.
        preview = "\n".join(f"          > {rec.id}" for rec in records[:5])
        if len(records) > 5:
            preview += f"\n          ... and {len(records) - 5} more."

        raise ValueError(
            f"'{filepath}' contains {len(records)} sequences.\n"
            f"        ORCA requires a single-sequence FASTA file.\n"
            f"\n"
            f"        A valid single-sequence FASTA looks like this:\n"
            f"          >NM_001301717.1 Homo sapiens ...\n"
            f"          ATGCGATCGATCGATCG...\n"
            f"\n"
            f"        Your file has {len(records)} entries starting with:\n"
            f"{preview}\n"
            f"\n"
            f"        Please extract a single sequence into its own FASTA file\n"
            f"        and re-run ORCA with that file."
        )

    record = records[0]
    return record.id, str(record.seq)


def validate_dna_sequence(sequence: str) -> tuple[bool, str]:
    """
    Validate and clean a raw DNA sequence string.

    Cleaning is applied in the following order:

    1. Uppercase and strip all whitespace and digits.
    2. Remove any characters outside the full IUPAC nucleotide alphabet
       (A T G C R Y S W K M B D H V N).
    3. Replace remaining IUPAC ambiguity codes (anything that is not a
       plain A, T, G, or C) with ``'N'``.  Replacement rather than
       deletion preserves the reading frame — deleting bases would shift
       every codon downstream.  The codon scanner and RSCU calculator
       both naturally skip any codon that contains ``'N'``.
    4. Reject sequences shorter than 6 bp (too short for any ORF).

    Parameters
    ----------
    sequence : str
        Raw DNA sequence string (any case; may contain whitespace or
        non-standard characters).

    Returns
    -------
    tuple[bool, str]
        ``(True, clean_seq)`` if the sequence passes validation,
        ``(False, partial_seq)`` otherwise.
    """
    sequence = re.sub(r"[\s\d]", "", sequence).upper()

    if not sequence:
        print("[VALIDATION] Sequence is empty after initial cleaning.")
        return False, ""

    invalid_chars = set(re.findall(r"[^ATGCRYSWKMBDHVN]", sequence))
    if invalid_chars:
        sequence = re.sub(r"[^ATGCRYSWKMBDHVN]", "", sequence)

    ambiguous_chars = set(re.findall(r"[^ATGC]", sequence))
    if ambiguous_chars:
        print(
            f"[VALIDATION] IUPAC ambiguity codes detected "
            f"({', '.join(sorted(ambiguous_chars))}); replacing with 'N' "
            f"to preserve reading frame."
        )
        sequence = re.sub(r"[^ATGC]", "N", sequence)

    if len(sequence) < 6:
        print(
            f"[VALIDATION] Cleaned sequence too short ({len(sequence)} bp). "
            "Minimum required: 6 bp."
        )
        return False, sequence

    print(f"[VALIDATION] Sequence is valid — {len(sequence)} bp ready for analysis.")
    return True, sequence


def write_cleaned_fasta(
    clean_seq:   str,
    accession:   str,
    output_path: str = "output/cleaned_sequence.fasta",
) -> None:
    """
    Write a single cleaned DNA sequence to a FASTA file.

    Parameters
    ----------
    clean_seq : str
        Validated, cleaned DNA string (ACGT and N only).
    accession : str
        Accession number or record ID used as the FASTA header.
    output_path : str, optional
        Destination file path.  The parent directory is created if it
        does not already exist.  Defaults to
        ``'output/cleaned_sequence.fasta'``; overridden by :func:`run`.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    record = SeqRecord(
        Seq(clean_seq),
        id=accession,
        description="cleaned sequence | ready for ORF analysis",
    )
    with open(output_path, "w") as fh:
        SeqIO.write(record, fh, "fasta")


def validate_start_codons(requested: list[str]) -> list[str]:
    """
    Uppercase and validate a list of user-supplied start codons.

    Parameters
    ----------
    requested : list[str]
        Start codon strings as supplied on the command line,
        e.g. ``['ATG']`` or ``['ATG', 'GTG']``.

    Returns
    -------
    list[str]
        The same list with every codon uppercased.

    Raises
    ------
    ValueError
        If any codon is not one of ``ATG``, ``GTG``, or ``TTG``.
        The caller (CLI entry point) is responsible for catching this
        and handling it appropriately (e.g. printing an error and
        calling ``sys.exit``).
    """
    upper = [c.upper() for c in requested]
    unknown = [c for c in upper if c not in VALID_START_CODONS]
    if unknown:
        raise ValueError(
            f"Unrecognised start codon(s): {', '.join(unknown)}\n"
            f"        Allowed values are: ATG, GTG, TTG"
        )
    return upper


def fetch_and_validate_one(
    accession:    str,
    output_fasta: str,
    fasta_file:   str | None = None,
) -> tuple[str, str] | tuple[None, None]:
    """
    Fetch or load one sequence, validate it, and write the cleaned FASTA.

    Internal helper called by :func:`run` to avoid duplicating the
    fetch-validate-write logic across single and comparative modes.

    Parameters
    ----------
    accession : str
        NCBI accession number.  Ignored when *fasta_file* is provided;
        the record ID from the FASTA header is used instead.
    output_fasta : str
        Destination file path for the cleaned FASTA.
    fasta_file : str or None, optional
        Path to a local FASTA file.  When provided, the sequence is
        loaded from disk rather than fetched from NCBI.

    Returns
    -------
    tuple[str, str]
        ``(accession, clean_seq)`` on success.
    tuple[None, None]
        On failure.
    """
    if fasta_file is not None:
        try:
            record_id, raw_sequence = load_fasta_from_file(fasta_file)
        except ValueError as e:
            # load_fasta_from_file raises ValueError for multi-sequence files.
            # Print the message here so the error is visible, then signal
            # failure to the caller via the (None, None) sentinel.
            print(f"\n[ERROR] {e}\n")
            return None, None

        if raw_sequence is None:
            return None, None
        accession = record_id
    else:
        raw_sequence = fetch_fasta_from_ncbi(accession)
        if raw_sequence is None:
            return None, None

    is_valid, clean_seq = validate_dna_sequence(raw_sequence)
    if not is_valid:
        print(f"[ERROR] Sequence '{accession}' failed validation. Aborting.")
        return None, None

    write_cleaned_fasta(clean_seq, accession, output_fasta)
    return accession, clean_seq

def run(
    accession:   str,
    email:       str,
    outdir:      str = "output",
    fasta_file:  str | None = None,
    comparative: bool = False,
    seq_num:     int  = 1,
) -> tuple[str, str] | tuple[None, None]:
    """
    Main pipeline entry point called by ``main.py``.

    Fetches or loads a single DNA sequence, cleans it, and writes the
    cleaned FASTA file to *outdir*.  The output filename reflects both
    the run mode and the sequence number:

    - Single mode  → ``cleaned_sequence_1.fasta``
    - Comparative  → ``comp_cleaned_sequence_1.fasta`` / ``comp_cleaned_sequence_2.fasta``

    Email validation is the caller's responsibility and should be done
    once before the pipeline starts, not per-sequence inside this function.

    Parameters
    ----------
    accession : str
        NCBI accession number.  Ignored when *fasta_file* is provided.
    email : str
        User email address passed directly to NCBI Entrez.  Must already
        be validated by the caller before this function is invoked.
    outdir : str, optional
        Directory for all output files.  Defaults to ``'output/'``.
    fasta_file : str or None, optional
        Path to a local FASTA file.  Overrides *accession* when provided.
    comparative : bool, optional
        When ``True``, the ``comp_`` prefix is added to the output
        filename to signal that this sequence is part of a comparative
        run.  Defaults to ``False``.
    seq_num : int, optional
        Sequence number (1 or 2) used in the output filename.
        Defaults to ``1``.

    Returns
    -------
    tuple[str, str]
        ``(accession, clean_seq)`` on success.
    tuple[None, None]
        On failure.
    """
    Entrez.email = email.strip()

    if comparative:
        filename = f"comp_cleaned_sequence_{seq_num}.fasta"
    else:
        filename = f"cleaned_sequence_{seq_num}.fasta"

    output_fasta = os.path.join(outdir, filename)

    acc, seq = fetch_and_validate_one(accession, output_fasta, fasta_file)
    if acc is None:
        return None, None

    return acc, seq
