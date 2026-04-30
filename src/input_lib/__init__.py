"""
input_lib — Input validation and NCBI sequence fetching package.

What this folder contains:
    input_validate.py — fetches DNA sequences from NCBI, validates them
                        against the IUPAC alphabet, removes invalid characters,
                        and writes cleaned FASTA files to disk.

Tells Python that input_lib/ is a package so other modules can import from it:
    from src.input_lib.input_validate import run
"""
