"""
ORCA — ORF Recognition and Comparative Annotator
=================================================
Main package initialisation file for the src/ folder.

What this file does:
    Tells Python that the src/ folder is a package, and exposes
    the most commonly used functions so main.py can import them
    cleanly without needing to know the full subfolder path.

Without this file:
    from src.input_lib.input_validate import run   # would fail
With this file:
    from src import run                             # works
"""

# ── Input validation ──────────────────────────────────────────────────────────
# Import the three main functions from input_lib/input_validate.py
# run()                  — fetches and validates sequences from NCBI
# validate_start_codons()— checks that start codons are valid (ATG/GTG/TTG)
# validate_email()       — checks the user email format before NCBI queries
from src.input_validate import (
    run,
    validate_start_codons,
    validate_email,
)

# ── Graphics ──────────────────────────────────────────────────────────────────
# Import the plot functions from graphics_lib/graphics.py
# plot_orf_map()            — generates the ORF map for a single sequence
# plot_comparative_orf_map()— generates a side-by-side ORF map for two sequences
from src.graphics import (
    plot_orf_map,
    plot_comparative_orf_map,
)

# ── Package metadata ──────────────────────────────────────────────────────────
__version__ = "1.0.0"
__author__  = "Tahmid Anwar, Nicole Decocker, Amanda Yaworsky"

# __all__ controls what gets exported when someone does: from src import *
# Only the functions listed here will be available that way.
__all__ = [
    "run",
    "validate_start_codons",
    "validate_email",
    "plot_orf_map",
    "plot_comparative_orf_map",
]
