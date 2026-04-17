#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graphics.py

Nicole Decocker's part

Purpose:
    Visual output for the ORCA pipeline.
    Generates an ORF map showing all non-nested ORFs across all six reading
    frames, and (in comparative mode) an RSCU codon-usage heatmap.
"""

from __future__ import annotations

from collections import defaultdict
from typing import List

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Codon colour by start codon (used for rectangles; no legend shown) ──────
CODON_COLORS = {
    "ATG": "#e05c5c",
    "GTG": "#5c8ee0",
    "TTG": "#5cc47a",
}
DEFAULT_COLOR = "#aaaaaa"

# ── Frame layout ─────────────────────────────────────────────────────────────
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

# ── Genetic code for RSCU ────────────────────────────────────────────────────
_GENETIC_CODE: dict[str, str] = {
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

# Amino acids with ≥2 synonymous codons — ordered for display
_AA_ORDER = [
    "Phe", "Leu", "Ile", "Val",
    "Ser", "Pro", "Thr", "Ala",
    "Tyr", "His", "Gln", "Asn",
    "Lys", "Asp", "Glu", "Cys",
    "Arg", "Gly",
]


# ── Internal helpers ─────────────────────────────────────────────────────────

def _filter_nested(flat_list: list) -> list:
    """
    Return a copy of flat_list with nested ORFs removed.

    An ORF is considered nested if its entire [start, end] span falls
    completely within the span of another ORF on the same strand and frame.
    ORFs with end=None (open reading frames) are passed through unchanged.
    """
    # Separate complete ORFs from open ones
    complete = [o for o in flat_list if o.get("end") is not None]
    open_orfs = [o for o in flat_list if o.get("end") is None]

    by_frame: dict = defaultdict(list)
    for orf in complete:
        by_frame[(orf["strand"], orf["frame"])].append(orf)

    non_nested: list = []
    for orfs in by_frame.values():
        for orf in orfs:
            is_nested = any(
                other is not orf
                and other["start"] <= orf["start"]
                and other["end"] >= orf["end"]
                for other in orfs
            )
            if not is_nested:
                non_nested.append(orf)

    return non_nested + open_orfs


def _compute_rscu(sequence: str) -> dict[str, float]:
    """
    Compute Relative Synonymous Codon Usage (RSCU) for a DNA sequence.

    RSCU(codon) = observed_count / (total_AA_count / num_synonymous_codons)

    A value of 1.0 means equal use of all synonymous codons.
    Values >1 indicate preference; values <1 indicate avoidance.
    Met, Trp, and Stop codons are excluded (no synonymous alternatives).
    """
    seq = sequence.upper()
    counts: dict[str, int] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3 and all(c in "ACGT" for c in codon):
            counts[codon] = counts.get(codon, 0) + 1

    # Group synonymous codons per amino acid
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in _GENETIC_CODE.items():
        if aa not in ("Stop", "Met", "Trp"):
            aa_codons[aa].append(codon)

    rscu: dict[str, float] = {}
    for aa, codons in aa_codons.items():
        total = sum(counts.get(c, 0) for c in codons)
        n = len(codons)
        for codon in codons:
            if total == 0:
                rscu[codon] = 0.0
            else:
                rscu[codon] = counts.get(codon, 0) / (total / n)
    return rscu


# ── ORF map ──────────────────────────────────────────────────────────────────

def draw_orf_map(
    ax,
    flat_list: list,
    seq_len:   int,
    accession: str,
) -> None:
    """
    Draw one ORF map onto an existing Axes object.

    Nested ORFs are automatically excluded.  Each ORF is drawn as a
    filled rectangle coloured by start codon; its ordinal number (1, 2, …)
    is printed in white text at the centre of the rectangle.
    """
    # Filter out nested ORFs then sort by start position for stable numbering
    display_list = sorted(
        _filter_nested(flat_list),
        key=lambda o: o.get("start", 0),
    )

    ax.set_xlim(0, seq_len)
    ax.set_ylim(-0.6, 6.4)

    # Alternating row backgrounds
    for y in range(6):
        bg = "#f5f5f5" if y % 2 == 0 else "#ebebeb"
        ax.axhspan(y - 0.4, y + 0.4, color=bg, zorder=0)

    # Strand divider
    ax.axhline(y=2.5, color="#999999", linewidth=1, linestyle="--", zorder=1)

    for orf_num, orf in enumerate(display_list, start=1):
        start  = orf["start"]
        end    = orf.get("end")
        strand = orf["strand"]
        frame  = orf["frame"]
        codon  = orf.get("start_codon", "")

        if end is None:
            continue

        y     = FRAME_Y[(strand, frame)]
        color = CODON_COLORS.get(codon, DEFAULT_COLOR)
        width = end - start

        # Draw rectangle
        rect = mpatches.Rectangle(
            xy     = (start, y - 0.275),
            width  = width,
            height = 0.55,
            color  = color,
            alpha  = 0.85,
            zorder = 2,
        )
        ax.add_patch(rect)

        # ORF number centred in the rectangle (only when wide enough to read)
        min_label_width = seq_len * 0.012
        if width >= min_label_width:
            fontsize = max(5, min(8, int(width / seq_len * 600)))
            ax.text(
                (start + end) / 2,
                y,
                str(orf_num),
                ha        = "center",
                va        = "center",
                fontsize  = fontsize,
                color     = "white",
                fontweight= "bold",
                zorder    = 3,
            )

    # Y-axis: frame labels
    ax.set_yticks(list(range(6)))
    ax.set_yticklabels(
        [FRAME_LABELS[(s, f)] for (s, f) in sorted(FRAME_Y, key=FRAME_Y.get)],
        fontsize=9,
    )

    ax.set_xlabel("Nucleotide position (bp)", fontsize=9)
    ax.set_title(f"{accession}  —  {seq_len} bp", fontsize=10, fontweight="bold")

    # Strand annotations on the left
    ax.text(-seq_len * 0.01, 4, "(+) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90)
    ax.text(-seq_len * 0.01, 1, "(−) strand", va="center", ha="right",
            fontsize=8, color="#444444", rotation=90)


def make_legend() -> list:
    """Return legend handles for start codon colours (kept for external use)."""
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
    """Plot a single-sequence ORF map and save to file."""
    fig, ax = plt.subplots(figsize=(14, 4))
    draw_orf_map(ax, flat_list, seq_len, accession)
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

    The two maps are stacked vertically.  Each uses its own x-axis scale
    so that ORF positions reflect the true sequence coordinates.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=False)

    draw_orf_map(ax1, flat1, seq_len1, acc1)
    draw_orf_map(ax2, flat2, seq_len2, acc2)

    fig.suptitle("Comparative ORF Map", fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Comparative ORF map saved to: {output_path}")


# ── RSCU codon-usage comparison (comparative mode only) ─────────────────────

def plot_codon_usage_comparison(
    seq1:        str,
    acc1:        str,
    seq2:        str,
    acc2:        str,
    output_path: str,
) -> None:
    """
    Plot a comparative RSCU heatmap for two DNA sequences.

    Only meaningful (and only called) in comparative mode.

    The heatmap rows are the two accessions; columns are all 59 sense codons
    (Met and Trp excluded — no synonymous alternatives; stop codons excluded)
    grouped by amino acid.  Cell colour encodes RSCU:

        1.0  = equal use of all synonymous codons for that amino acid
        > 1  = this codon is preferred
        < 1  = this codon is avoided

    This reveals codon bias differences between the two sequences at a glance
    — useful for comparing organisms with different AT/GC content, or for
    spotting expression-relevant bias in genes from the same organism.
    """
    # ── Build ordered codon list ──────────────────────────────────────────
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in _GENETIC_CODE.items():
        if aa not in ("Stop", "Met", "Trp"):
            aa_codons[aa].append(codon)

    for aa in aa_codons:
        aa_codons[aa].sort()

    codon_list: list[str] = []
    aa_boundaries: list[tuple[int, int, str]] = []  # (col_start, col_end, aa_name)
    pos = 0
    for aa in _AA_ORDER:
        codons = aa_codons[aa]
        aa_boundaries.append((pos, pos + len(codons), aa))
        codon_list.extend(codons)
        pos += len(codons)

    # ── Compute RSCU for both sequences ──────────────────────────────────
    rscu1 = _compute_rscu(seq1)
    rscu2 = _compute_rscu(seq2)

    data = np.array([
        [rscu1.get(c, 0.0) for c in codon_list],
        [rscu2.get(c, 0.0) for c in codon_list],
    ])

    # ── Figure layout: extra bottom space for colorbar ────────────────────
    fig, ax = plt.subplots(figsize=(22, 3.2))
    fig.subplots_adjust(top=0.72, bottom=0.38, left=0.07, right=0.97)

    im = ax.imshow(
        data,
        aspect = "auto",
        cmap   = "RdYlBu_r",   # red = high RSCU (preferred), blue = low (avoided)
        vmin   = 0.0,
        vmax   = 3.0,
    )

    # ── Black horizontal line separating the two sequences ────────────────
    ax.axhline(y=0.5, color="black", linewidth=1.0, zorder=4)

    # ── White vertical lines separating amino acid groups ─────────────────
    for col_start, _, _ in aa_boundaries:
        if col_start > 0:
            ax.axvline(x=col_start - 0.5, color="white", linewidth=2.5, zorder=3)

    # ── Y axis: sequence names ────────────────────────────────────────────
    ax.set_yticks([0, 1])
    ax.set_yticklabels([acc1, acc2], fontsize=9, fontweight="bold")
    ax.tick_params(axis="y", length=0)

    # ── Bottom x axis: codon labels ───────────────────────────────────────
    ax.set_xticks(range(len(codon_list)))
    ax.set_xticklabels(codon_list, rotation=90, fontsize=6.5, family="monospace")
    ax.tick_params(axis="x", length=2, pad=1)

    # ── Top x axis: amino acid group labels ───────────────────────────────
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    aa_mid   = [(s + e - 1) / 2 for s, e, _ in aa_boundaries]
    aa_names = [aa               for _, _, aa in aa_boundaries]
    ax_top.set_xticks(aa_mid)
    ax_top.set_xticklabels(aa_names, fontsize=8, fontweight="bold", rotation=45, ha="left")
    ax_top.tick_params(top=False, pad=2)

    # ── Title ─────────────────────────────────────────────────────────────
    ax.set_title(
        f"Codon Usage Bias  —  {acc1}  vs  {acc2}  (RSCU)",
        fontsize=10,
        fontweight="bold",
        pad=30,
    )

    # ── Horizontal colorbar below the heatmap ─────────────────────────────
    cbar_ax = fig.add_axes([0.15, 0.10, 0.70, 0.07])   # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("RSCU  (1.0 = equal synonymous codon usage)", fontsize=8, labelpad=4)
    cbar.ax.tick_params(labelsize=7)

    # Dashed line at RSCU = 1.0 (equal usage baseline); colorbar data range is [0, 3]
    cbar.ax.axvline(x=1.0, color="black", linewidth=1.5, linestyle="--", zorder=5)
    cbar.ax.text(
        1.0, 1.7, "1.0",
        ha="center", va="bottom", fontsize=7, fontweight="bold",
        color="black", transform=cbar.ax.get_xaxis_transform(),
    )

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Codon usage comparison saved to: {output_path}")
