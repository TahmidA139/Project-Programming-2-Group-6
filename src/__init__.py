"""ORCA — ORF Recognition and Comparative Analysis."""

from src.input_validate import run, validate_start_codons, validate_email
from src.graphics import plot_orf_map, plot_comparative_orf_map

__version__ = "1.0.0"
__all__ = ["run", "validate_start_codons", "validate_email",
           "plot_orf_map", "plot_comparative_orf_map"]
