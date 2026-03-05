#!/usr/bin/env python3

"""

Tahmid's portion:

input_validate.py
-----------------
Fetches a DNA sequence from NCBI in FASTA format, validates it (removing
any invalid characters), and the cleaned DNA sequence will be passed directly to ORF_finder.py through main.py

Dependency:
    pip install biopython
"""

#import necessary packages
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# NCBI requires a valid email address for all Entrez queries
Entrez.email = "your_email@example.com"   # replace before running


# ─────────────────────────────────────────────────────────────────────────────
# FUNCTION 1 — Fetch FASTA from NCBI
# ─────────────────────────────────────────────────────────────────────────────

def fetch_fasta_from_ncbi(accession: str, db: str = "nucleotide") -> str | None:
    """
    Objective:
        Download a DNA sequence in FASTA format from NCBI.

    Input:
        accession (str): NCBI nucleotide accession number
                         e.g. 'NM_001301717' or 'NC_000913.3'
        db (str):        Entrez database to query (default: 'nucleotide')

    Output:
        sequence (str):  Raw DNA sequence string retrieved from NCBI,
                         or None if the fetch fails.

    High-Level Steps:
        1. Query NCBI using the accession number via Entrez.efetch
        2. Retrieve and parse the FASTA record with SeqIO
        3. Extract and return the DNA sequence as a plain string
    """
    try:
        print(f"[INFO] Querying NCBI for accession: '{accession}' ...")
        # Step 1 & 2 — query NCBI and parse the FASTA response
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") # Entrez.efetch is the actual API call to NCBI
        # db: which NCBI database to search (e.g. 'nucleotide') id: the accession number to look up, rettype: the file format to return, retmode: how the data is encoded; 'text' means plain readable text,
        
        record = SeqIO.read(handle, "fasta") # 'handle' is a file-like object
        handle.close() # Always close the handle as leaving handles open can cause memory leaks or hit NCBI's rate limits

    # try/except is used here because many things can go wrong such as bad accession number, no internet, NCBI is down, rate limit hit, etc.
    # Instead of crashing the whole program, any error is found and return "None" so the caller can handle the failure gracefully.

    # Step 3 — extract the raw sequence string
        sequence = str(record.seq)
        print(f"[INFO] Fetched '{record.id}' — {len(sequence)} bp")
        return sequence
    # record.seq is a Biopython Seq object (not a plain Python string). We wrap it in str() to convert it to a regular string, which takes only the sequence.

    except Exception as e:
        print(f"[ERROR] Could not fetch sequence from NCBI: {e}")
        return None
    # 'Exception as e' catches any type of error and stores its message in 'e'. We print the error message so the user knows what went wrong,
    # then return None to signal that no sequence was retrieved.


