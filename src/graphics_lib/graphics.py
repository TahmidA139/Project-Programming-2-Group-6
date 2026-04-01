#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graphics.py

Purpose:
    Visual output for the ORCA pipeline.
    Generates a genome-browser style ORF map showing all ORFs across
    all six reading frames.
"""

from __future__ import annotations

from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np

CODON_COLORS = {
    "ATG": "#c0392b",   # rich red
    "GTG": "#2471a3",   # steel blue
    "TTG": "#1e8449",   # forest green
}
DEFAULT_COLOR = "#7f8c8d"

# Lane background colours — (+) strand cool blue-grey, (−) strand warm rose-grey
LANE_COLORS = {
    "+": "#d6e4f0",   # soft blue
    "-": "#fadbd8",   # soft rose
}
LANE_EDGE_COLOR = "#ffffff"   # white border between lanes

FRAME_ROW = {
    ("+", 0): 5,
    ("+", 1): 4,
    ("+", 2): 3,
    ("-", 0): 2,
    ("-", 1): 1,
    ("-", 2): 0,
}

FRAME_LABEL = {
    ("+", 0): "+1",
    ("+", 1): "+2",
    ("+", 2): "+3",
    ("-", 0): "−1",
    ("-", 1): "−2",
    ("-", 2): "−3",
}

ROW_HEIGHT   = 1.0   # each lane occupies 1 unit on the y axis
BAR_HEIGHT   = 0.55  # rectangle height within the lane
ARROW_FRAC   = 0.06  # arrow head occupies this fraction of ORF width (min 12 px)

def _draw_lane_backgrounds(ax, seq_len: int) -> None:
    """Draw coloured background bands for each reading-frame lane."""
    for (strand, frame), row in FRAME_ROW.items():
        color = LANE_COLORS[strand]
        y_lo  = row - 0.5
        y_hi  = row + 0.5
        ax.axhspan(y_lo, y_hi, facecolor=color, edgecolor=LANE_EDGE_COLOR,
                   linewidth=0.8, zorder=0)


def _draw_strand_divider(ax) -> None:
    """Dashed line between (+) and (−) strand groups."""
    ax.axhline(y=2.5, color="#555555", linewidth=1.2, linestyle="--",
               zorder=1, alpha=0.6)


def _arrow_inside(ax, x0: float, x1: float, row: int,
                  strand: str, color: str) -> None:
    """
    Draw a small directional arrow INSIDE the ORF rectangle.

    The arrow sits near the right end for (+) strand and near the
    left end for (−) strand, pointing in the direction of transcription.
    """
    width  = x1 - x0
    hw     = BAR_HEIGHT * 0.42          # arrow head half-width (y direction)
    hl     = max(width * ARROW_FRAC, 12)  # arrow head length (x direction)
    hl     = min(hl, width * 0.35)       # never bigger than 35 % of ORF

    if strand == "+":
        ax_x  = x1 - hl * 1.1          # tip near right edge
        dx    = hl
    else:
        ax_x  = x0 + hl * 1.1          # tip near left edge
        dx    = -hl

    arrow = mpatches.FancyArrow(
        x      = ax_x - dx,            # tail position
        y      = row,
        dx     = dx,
        dy     = 0,
        width  = hw * 0.6,
        length_includes_head = True,
        head_width  = hw,
        head_length = abs(dx),
        color  = "white",
        alpha  = 0.75,
        zorder = 4,
    )
    ax.add_patch(arrow)


def _draw_orfs(ax, flat_list: list, seq_len: int) -> None:
    """Draw one rectangle + inner arrow + label per ORF."""
    for orf in flat_list:
        start  = orf.get("start")
        end    = orf.get("end")
        strand = orf.get("strand")
        frame  = orf.get("frame")
        codon  = orf.get("start_codon", "ATG")
        label  = orf.get("orf_id", "")

        if end is None or start is None:
            continue
        if (strand, frame) not in FRAME_ROW:
            continue

        row    = FRAME_ROW[(strand, frame)]
        color  = CODON_COLORS.get(codon, DEFAULT_COLOR)
        width  = end - start

        # Rectangle
        rect = mpatches.Rectangle(
            (start, row - BAR_HEIGHT / 2),
            width, BAR_HEIGHT,
            linewidth = 0.6,
            edgecolor = "white",
            facecolor = color,
            alpha     = 0.88,
            zorder    = 2,
        )
        ax.add_patch(rect)

        # Directional arrow inside the rectangle
        _arrow_inside(ax, start, end, row, strand, color)

        # ORF ID label — centred in the rectangle, white text
        if width > seq_len * 0.02:          # only label if wide enough to read
            font_size = max(5, min(8, width / seq_len * 120))
            ax.text(
                (start + end) / 2, row,
                label,
                ha        = "center",
                va        = "center",
                fontsize  = font_size,
                color     = "white",
                fontweight= "bold",
                clip_on   = True,
                zorder    = 5,
                path_effects = [
                    pe.withStroke(linewidth=1.2, foreground=color)
                ],
            )


def _configure_axes(ax, seq_len: int, accession: str) -> None:
    """Set axis limits, tick labels, title, and strand annotations."""
    n_rows = 6
    ax.set_xlim(0, seq_len)
    ax.set_ylim(-0.5, n_rows - 0.5)   # tight — no phantom row above +1

    # y-axis: one tick per row
    row_order = sorted(FRAME_ROW.values())               # [0,1,2,3,4,5]
    labels    = [FRAME_LABEL[k]
                 for k in sorted(FRAME_ROW, key=FRAME_ROW.get)]
    ax.set_yticks(row_order)
    ax.set_yticklabels(labels, fontsize=9, fontfamily="monospace")
    ax.tick_params(axis="y", length=0, pad=4)

    ax.set_xlabel("Nucleotide position (bp)", fontsize=9, labelpad=6)
    ax.set_title(f"{accession}  —  {seq_len:,} bp",
                 fontsize=10, fontweight="bold", pad=8)

    # Strand bracket annotations on the right margin
    ax.annotate(
        "(+) strand", xy=(1.002, (3 + 5.5) / 2 / n_rows),
        xycoords="axes fraction",
        fontsize=7.5, color="#1a5276", va="center",
        rotation=90,
    )
    ax.annotate(
        "(−) strand", xy=(1.002, (0 + 2.5) / 2 / n_rows),
        xycoords="axes fraction",
        fontsize=7.5, color="#922b21", va="center",
        rotation=90,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")


def _draw_orf_map(ax, flat_list: list, seq_len: int, accession: str) -> None:
    """Compose one complete ORF map on an existing Axes."""
    _draw_lane_backgrounds(ax, seq_len)
    _draw_strand_divider(ax)
    _draw_orfs(ax, flat_list, seq_len)
    _configure_axes(ax, seq_len, accession)


def _make_legend() -> list:
    """Legend patches for start codons."""
    return [
        mpatches.Patch(facecolor=color, edgecolor="white",
                       linewidth=0.5, label=codon)
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
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    _draw_orf_map(ax, flat_list, seq_len, accession)

    fig.legend(
        handles       = _make_legend(),
        title         = "Start codon",
        loc           = "upper right",
        fontsize      = 8,
        title_fontsize= 8,
        framealpha    = 0.9,
        edgecolor     = "#cccccc",
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
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

    The two panels use independent x-axis scales so each sequence fills
    the full width — ORF *shapes* are comparable even if lengths differ.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8),
                                   sharex=False,
                                   gridspec_kw={"hspace": 0.55})
    fig.patch.set_facecolor("#fafafa")
    for ax in (ax1, ax2):
        ax.set_facecolor("#fafafa")

    _draw_orf_map(ax1, flat1, seq_len1, acc1)
    _draw_orf_map(ax2, flat2, seq_len2, acc2)

    fig.legend(
        handles       = _make_legend(),
        title         = "Start codon",
        loc           = "upper right",
        fontsize      = 8,
        title_fontsize= 8,
        framealpha    = 0.9,
        edgecolor     = "#cccccc",
    )
    fig.suptitle("Comparative ORF Map", fontsize=12,
                 fontweight="bold", y=1.01)

    plt.savefig(output_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[INFO] Comparative ORF map saved to: {output_path}")
