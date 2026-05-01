#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graphics.py

Nicole Decocker's part

Purpose:
    Visual output for the ORCA pipeline.
    Generates an ORF map showing all ORFs across all six reading
    frames, and (in comparative mode) an RSCU codon-usage heatmap.

Public API
----------
draw_orf_map                Draw one ORF panel onto existing Axes.
make_legend                 Return legend Patch handles.
plot_orf_map                Single-sequence ORF map → file.
plot_comparative_orf_map    Two-sequence ORF map → file.
plot_codon_usage_comparison RSCU heatmap comparison → file.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, TypedDict

import numpy as np
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

# Genetic code for RSCU 
GENETIC_CODE: dict[str, str] = {
    "TTT": "Phe", "TTC": "Phe",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu",
    "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "AGT": "Ser", "AGC": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr",
    "TAA": "Stop", "TAG": "Stop", "TGA": "Stop",
    "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys",
    "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Amino acids with ≥2 synonymous codons — display order for RSCU heatmap
AA_ORDER: list[str] = [
    "Phe", "Leu", "Ile", "Val", "Ser", "Pro", "Thr", "Ala",
    "Tyr", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys", "Arg", "Gly",
]

# RSCU helpers
def compute_rscu(sequence: str) -> dict[str, float]:
    """Return RSCU values for all sense codons in *sequence*; Met, Trp, and stops excluded."""
    seq = sequence.upper()
    counts: dict[str, int] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3 and all(c in "ACGT" for c in codon):
            counts[codon] = counts.get(codon, 0) + 1

    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        if aa not in ("Stop", "Met", "Trp"):
            aa_codons[aa].append(codon)

    rscu: dict[str, float] = {}
    for aa, codons in aa_codons.items():
        total = sum(counts.get(c, 0) for c in codons)
        n = len(codons)
        for codon in codons:
            rscu[codon] = counts.get(codon, 0) / (total / n) if total else 0.0
    return rscu


def build_codon_order() -> tuple[list[str], list[tuple[int, int, str]]]:
    """Return (codon_list, aa_boundaries) in AA_ORDER for use in the RSCU heatmap."""
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        if aa not in ("Stop", "Met", "Trp"):
            aa_codons[aa].append(codon)
    for aa in aa_codons:
        aa_codons[aa].sort()

    codon_list: list[str] = []
    aa_boundaries: list[tuple[int, int, str]] = []
    pos = 0
    for aa in AA_ORDER:
        codons = aa_codons[aa]
        aa_boundaries.append((pos, pos + len(codons), aa))
        codon_list.extend(codons)
        pos += len(codons)
    return codon_list, aa_boundaries


def configure_heatmap_axes(
    ax: Axes,
    ax_top: Axes,
    codon_list: list[str],
    aa_boundaries: list[tuple[int, int, str]],
    acc1: str,
    acc2: str,
) -> None:
    """Apply codon ticks, amino-acid group labels, and sequence divider to the RSCU heatmap axes."""
    ax.axhline(y=0.5, color="black", linewidth=1.0, zorder=10)

    ax.set_yticks([0, 1])
    ax.set_yticklabels([acc1, acc2], fontsize=9, fontweight="bold")
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(range(len(codon_list)))
    ax.set_xticklabels(codon_list, rotation=90, fontsize=6.5, family="monospace")
    ax.tick_params(axis="x", length=2, pad=1)

    ax_top.set_xlim(ax.get_xlim())
    aa_mid   = [(s + e - 1) / 2 for s, e, _ in aa_boundaries]
    aa_names = [name              for _, _, name in aa_boundaries]
    ax_top.set_xticks(aa_mid)
    ax_top.set_xticklabels(aa_names, fontsize=8, fontweight="bold", rotation=45, ha="left")
    ax_top.tick_params(top=False, pad=2)


def add_rscu_colorbar(fig: Figure, im: Any) -> None:
    """Attach a horizontal RSCU colorbar with a dashed reference line at 1.0."""
    cbar_ax = fig.add_axes([0.15, 0.10, 0.70, 0.07])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("RSCU  (1.0 = equal synonymous codon usage)", fontsize=8, labelpad=4)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.axvline(x=1.0, color="black", linewidth=1.5, linestyle="--", zorder=5)
    cbar.ax.text(
        1.0, 1.7, "1.0",
        ha="center", va="bottom", fontsize=7, fontweight="bold",
        color="black", transform=cbar.ax.get_xaxis_transform(),
    )

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

def plot_codon_usage_comparison(
    flat1: list[dict[str, Any]],
    acc1:  str,
    seq1:  str,
    flat2: list[dict[str, Any]],
    acc2:  str,
    seq2:  str,
    output_path: str | Path,
) -> None:
    """
    Plot a comparative RSCU heatmap using only ORF sequences and save to file.

    Parameters
    ----------
    flat1 : list[dict[str, Any]]
        Flat ORF list for sequence 1, as returned by ``find_orfs()``.
        Open ORFs (``end=None``) are excluded automatically.
    acc1 : str
        Accession or label for *seq1*; used as the row label and in the title.
    seq1 : str
        Full forward-strand DNA for sequence 1; used to extract ORF sequences.
    flat2 : list[dict[str, Any]]
        Flat ORF list for sequence 2.
    acc2 : str
        Accession or label for *seq2*.
    seq2 : str
        Full forward-strand DNA for sequence 2.
    output_path : str | Path
        Destination file path (PNG, PDF, SVG, …).

    Ensures
    -------
    RSCU is computed from ORF nucleotides only — not the full background
    sequence — so the heatmap reflects coding-region codon bias.
    Only the 59 sense codons with synonymous alternatives are included.
    Met (ATG), Trp (TGG), and stop codons are excluded.
    The colour scale runs 0–3 RSCU units; a dashed reference line marks 1.0.

    Returns
    -------
    None
    """
    from src.orf_finder_lib.frame_scanner import extract_orf_sequence

    orf_seq1 = "".join(
        extract_orf_sequence(o, seq1)
        for o in flat1 if o.get("end") is not None
    )
    orf_seq2 = "".join(
        extract_orf_sequence(o, seq2)
        for o in flat2 if o.get("end") is not None
    )

    codon_list, aa_boundaries = build_codon_order()
    rscu1 = compute_rscu(orf_seq1)
    rscu2 = compute_rscu(orf_seq2)
    data = np.array([
        [rscu1.get(c, 0.0) for c in codon_list],
        [rscu2.get(c, 0.0) for c in codon_list],
    ])

    fig, ax = plt.subplots(figsize=(22, 3.2))
    fig.subplots_adjust(top=0.72, bottom=0.38, left=0.07, right=0.97)

    im = ax.imshow(data, aspect="auto", cmap="RdYlBu_r", vmin=0.0, vmax=3.0)
    ax_top = ax.twiny()
    configure_heatmap_axes(ax, ax_top, codon_list, aa_boundaries, acc1, acc2)
    add_rscu_colorbar(fig, im)

    # Draw thick black amino-acid dividers AFTER imshow so they sit on top.
    # We use ax.transData for the x position and ax.transAxes for y so the
    # lines can extend from the bottom of the heatmap up through ax_top labels.
    _, aa_boundaries_local = build_codon_order()
    for col_start, _, _ in aa_boundaries_local:
        if col_start > 0:
            x = col_start - 0.5
            # Line through the heatmap (data coords x, axes-fraction y)
            from matplotlib.transforms import blended_transform_factory as btf
            trans = btf(ax.transData, ax.transAxes)
            ax.plot(
                [x, x], [-0.20, 1.35],   # extend slightly below and above the heatmap
                color="black", linewidth=1.5, zorder=20,
                transform=trans, clip_on=False,
            )

    ax.set_title(
        f"Codon Usage Bias  —  {acc1}  vs  {acc2}  (RSCU)",
        fontsize=10, fontweight="bold", pad=30,
    )
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Codon usage comparison saved to: %s", output_path)
