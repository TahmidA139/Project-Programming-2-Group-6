#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graphics.py

Nicole Decocker's part

Purpose:
    Visual output for the ORCA pipeline.
    Generates an ORF map showing all ORFs across all six reading
    frames for single and comparative sequence modes.

Public API
----------
draw_orf_map                Draw one ORF panel onto existing Axes.
make_legend                 Return legend Patch handles.
plot_orf_map                Single-sequence ORF map → file.
plot_comparative_orf_map    Two-sequence ORF map → file.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, TypedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.transforms import blended_transform_factory

log = logging.getLogger(__name__)


# Typed ORF record 
class OrfRecord(TypedDict, total=False):
    """Shape of one ORF dict as produced by ``find_orfs()``."""
    strand:      str
    frame:       int
    start:       int
    end:         int | None
    start_codon: str


# Codon colours 
CODON_COLORS: dict[str, str] = {
    "ATG": "#e05c5c",
    "GTG": "#5c8ee0",
    "TTG": "#5cc47a",
}
DEFAULT_COLOR: str = "#aaaaaa"

#  Rendering constants 
MIN_LABEL_FRAC: float = 0.012  # min ORF width as fraction of seq_len to show label
FONT_SCALE: int = 600          # numerator used to scale in-rect label font size

#  Frame layout 
FRAME_Y: dict[tuple[str, int], int] = {
    ("+", 0): 5, ("+", 1): 4, ("+", 2): 3,
    ("-", 0): 2, ("-", 1): 1, ("-", 2): 0,}

FRAME_LABELS: dict[tuple[str, int], str] = {
    ("+", 0): "+1", ("+", 1): "+2", ("+", 2): "+3",
    ("-", 0): "-1", ("-", 1): "-2", ("-", 2): "-3",}

# ORF map helpers
def setup_frame_axes(ax: Axes, seq_len: int, accession: str) -> None:
    """Configure axis limits, row backgrounds, strand divider, and tick labels for one ORF panel."""
    ax.set_xlim(0, seq_len)
    ax.set_ylim(-0.6, 6.4)

    for y in range(6):
        bg = "#f5f5f5" if y % 2 == 0 else "#ebebeb"
        ax.axhspan(y - 0.4, y + 0.4, color=bg, zorder=0)

    ax.axhline(y=2.5, color="#999999", linewidth=1, linestyle="--", zorder=1)

    ax.set_yticks(list(range(6)))
    ax.set_yticklabels(
        [FRAME_LABELS[(s, f)] for (s, f) in sorted(FRAME_Y, key=FRAME_Y.get)],
        fontsize=9,
    )
    ax.set_xlabel("Nucleotide position (bp)", fontsize=9)
    ax.set_title(f"{accession}  —  {seq_len} bp", fontsize=10, fontweight="bold")

    # Blended transform (axes-fraction x, data y) keeps strand labels clear of
    # the frame tick labels regardless of figure width or DPI.
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(-0.03, 4, "(+) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90, transform=trans)
    ax.text(-0.03, 1, "(−) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90, transform=trans)

def draw_single_orf(
    ax: Axes,
    orf: dict[str, Any],
    orf_num: int,
    seq_len: int,
) -> None:
    """Draw one ORF rectangle onto *ax*, with an ordinal label if the rectangle is wide enough."""
    start  = orf["start"]
    end    = orf["end"]
    strand = orf["strand"]
    frame  = orf["frame"]
    codon  = orf.get("start_codon", "")

    y     = FRAME_Y[(strand, frame)]
    color = CODON_COLORS.get(codon, DEFAULT_COLOR)
    width = end - start

    ax.add_patch(mpatches.Rectangle(
        xy=(start, y - 0.275), width=width, height=0.55,
        color=color, alpha=0.85, zorder=2,
    ))

    if width >= seq_len * MIN_LABEL_FRAC:
        fontsize = max(5, min(8, int(width / seq_len * FONT_SCALE)))
        ax.text(
            (start + end) / 2, y, str(orf_num),
            ha="center", va="center", fontsize=fontsize,
            color="white", fontweight="bold", zorder=3,
        )

# Public API
def draw_orf_map(
    ax: Axes,
    flat_list: list[dict[str, Any]],
    seq_len: int,
    accession: str,
) -> None:
    """
    Draw one ORF map onto an existing Axes object.

    Parameters
    ----------
    ax : Axes
        Target axes; must already exist. This function does not create figures.
    flat_list : list[dict[str, Any]]
        Flat ORF list as returned by ``find_orfs()``. Each dict must contain
        ``strand``, ``frame``, ``start``, ``end``, and ``start_codon``.
    seq_len : int
        Total sequence length in bp; sets the x-axis scale so ORF
        positions are shown in true sequence coordinates.
    accession : str
        Accession number or label shown in the subplot title.

    Ensures
    -------
    ORFs are drawn in ascending ``start`` order for stable ordinal numbering.
    Open ORFs (``end=None``) are silently excluded before drawing.

    Returns
    -------
    None
        Draws directly onto *ax*; the caller is responsible for saving
        or displaying the figure.
    """
    display_list = sorted(
        [o for o in flat_list if o.get("end") is not None],
        key=lambda o: o.get("start", 0),
    )
    setup_frame_axes(ax, seq_len, accession)
    for orf_num, orf in enumerate(display_list, start=1):
        draw_single_orf(ax, orf, orf_num, seq_len)


def make_legend() -> list[mpatches.Patch]:
    """
    Build legend handles for start-codon colours.

    Returns
    -------
    list[mpatches.Patch]
        One ``Patch`` per entry in ``CODON_COLORS``, labelled by codon name.
    """
    return [
        mpatches.Patch(color=color, label=codon)
        for codon, color in CODON_COLORS.items()]

def plot_orf_map(
    flat_list: list[dict[str, Any]],
    seq_len: int,
    accession: str,
    output_path: str | Path,
) -> None:
    """
    Plot a single-sequence ORF map and save to *output_path*.

    Parameters
    ----------
    flat_list : list[dict[str, Any]]
        ORF records for the sequence; see :func:`draw_orf_map`.
    seq_len : int
        Sequence length in bp.
    accession : str
        Accession or label for the subplot title.
    output_path : str | Path
        Destination file path. The format is inferred from the extension
        (PNG, PDF, SVG, …).

    Returns
    -------
    None
    """
    fig, ax = plt.subplots(figsize=(14, 4))
    draw_orf_map(ax, flat_list, seq_len, accession)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("ORF map saved to: %s", output_path)

def plot_comparative_orf_map(
    flat1: list[dict[str, Any]],
    seq_len1: int,
    acc1: str,
    flat2: list[dict[str, Any]],
    seq_len2: int,
    acc2: str,
    output_path: str | Path,
) -> None:
    """
    Plot a two-sequence comparative ORF map and save to *output_path*.

    Parameters
    ----------
    flat1 : list[dict[str, Any]]
        ORF records for the first sequence; see :func:`draw_orf_map`.
    seq_len1 : int
        Length of the first sequence in bp.
    acc1 : str
        Accession or label for the first subplot.
    flat2 : list[dict[str, Any]]
        ORF records for the second sequence.
    seq_len2 : int
        Length of the second sequence in bp.
    acc2 : str
        Accession or label for the second subplot.
    output_path : str | Path
        Destination file path.

    Ensures
    -------
    Each subplot uses an independent x-axis (``sharex=False``) so that ORF
    positions reflect true sequence coordinates even when the two sequences
    differ substantially in length.

    Returns
    -------
    None
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=False)
    draw_orf_map(ax1, flat1, seq_len1, acc1)
    draw_orf_map(ax2, flat2, seq_len2, acc2)
    fig.suptitle("Comparative ORF Map", fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Comparative ORF map saved to: %s", output_path)
