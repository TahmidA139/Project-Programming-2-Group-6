"""
graphics_lib — Visualisation package.

What this folder contains:
    graphics.py — generates all plots for the ORCA pipeline using matplotlib.
                  Includes the ORF map (orf_map.png), codon usage charts,
                  strand distribution bar chart, and comparative ORF map
                  when two sequences are analysed side by side.

Tells Python that graphics_lib/ is a package so other modules can import from it:
    from src.graphics_lib.graphics import plot_orf_map
    from src.graphics_lib.graphics import plot_comparative_orf_map
"""
