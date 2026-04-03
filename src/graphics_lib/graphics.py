#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graphics.py

Nicole Decocker's part 
 
Purpose:
    Visual output for the ORCA pipeline.
    Generates a ORF map showing all ORFs across all six reading frames.

    again unitll all is done im not messing with the docstrings sorry
"""

from __future__ import annotations

from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow

CODON_COLORS = {
    "ATG": "#e05c5c",
    "GTG": "#5c8ee0",
    "TTG": "#5cc47a",
}
DEFAULT_COLOR = "#aaaaaa"

FRAME_Y = {
    ("+", 0): 5,
    ("+", 1): 4,
    ("+", 2): 3,
    ("-", 0): 2,
    ("-", 1): 1,
    ("-", 2): 0,
}

FRAME_LABELS = {
    ("+", 0): "+1",
    ("+", 1): "+2",
    ("+", 2): "+3",
    ("-", 0): "-1",
    ("-", 1): "-2",
    ("-", 2): "-3",
}



def _draw_orf_map(
    ax,
    flat_list: list,
    seq_len:   int,
    accession: str,
) -> None:
    """Draw one ORF map onto an existing Axes object."""

    ax.set_xlim(0, seq_len)
    ax.set_ylim(-0.6, 6.4)


    for y in range(6):
        color = "#f5f5f5" if y % 2 == 0 else "#ebebeb"
        ax.axhspan(y - 0.4, y + 0.4, color=color, zorder=0)

    
    ax.axhline(y=2.5, color="#999999", linewidth=1, linestyle="--", zorder=1)


    for orf in flat_list:
        start  = orf["start"]
        end    = orf["end"]
        strand = orf["strand"]
        frame  = orf["frame"]
        codon  = orf["start_codon"]

        if end is None:
            continue

        y     = FRAME_Y[(strand, frame)]
        color = CODON_COLORS.get(codon, DEFAULT_COLOR)
        width = end - start

        # Arrow direction shows strand
        dx = width * 0.08 if strand == "+" else -width * 0.08
        arrow_x = end - dx if strand == "+" else start - dx

        rect = mpatches.FancyArrow(
            x      = start if strand == "+" else end,
            y      = y,
            dx     = width if strand == "+" else -width,
            dy     = 0,
            width  = 0.55,
            length_includes_head = True,
            head_width  = 0.65,
            head_length = max(width * 0.08, 15),
            color  = color,
            alpha  = 0.85,
            zorder = 2,
        )
        ax.add_patch(rect)


    ax.set_yticks(list(range(6)))
    ax.set_yticklabels(
        [FRAME_LABELS[(s, f)] for (s, f) in sorted(FRAME_Y, key=FRAME_Y.get)],
        fontsize=9,
    )


    ax.set_xlabel("Nucleotide position (bp)", fontsize=9)
    ax.set_title(f"{accession}  —  {seq_len} bp", fontsize=10, fontweight="bold")


    ax.text(-seq_len * 0.01, 4, "(+) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90)
    ax.text(-seq_len * 0.01, 1, "(−) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90)


def _make_legend() -> list:
    """Return legend handles for start codon colors."""
    return [
        mpatches.Patch(color=color, label=codon)
        for codon, color in CODON_COLORS.items()
    ]


def plot_orf_map(
    flat_list:   list,
    seq_len:     int,
    accession:   str,
    output_path: str,
) -> None:
    """
    Plot a single-sequence ORF map and save to file.

    """
    fig, ax = plt.subplots(figsize=(14, 4))
    _draw_orf_map(ax, flat_list, seq_len, accession)
    fig.legend(handles=_make_legend(), title="Start codon",
               loc="upper right", fontsize=8, title_fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[INFO] ORF map saved to: {output_path}")


def plot_comparative_orf_map(
    flat1:       list,
    seq_len1:    int,
    acc1:        str,
    flat2:       list,
    seq_len2:    int,
    acc2:        str,
    output_path: str,
) -> None:
    """
    Plot a two-sequence comparative ORF map and save to file.

    The two maps share the same x-axis scale (the longer sequence sets
    the maximum) so ORF positions are visually comparable.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=False)

    _draw_orf_map(ax1, flat1, seq_len1, acc1)
    _draw_orf_map(ax2, flat2, seq_len2, acc2)

    fig.legend(handles=_make_legend(), title="Start codon",
               loc="upper right", fontsize=8, title_fontsize=8)
    fig.suptitle("Comparative ORF Map", fontsize=12, fontweight="bold", y=1.01)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
