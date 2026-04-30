#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genetic_code_inference.py

Purpose:
    Evaluates which NCBI genetic code best explains the ORFs detected by the
    ORCA pipeline, using log-likelihood scoring of observed codon usage against
    each code's synonymous-group structure.

    This is the standard approach used by gene-prediction tools such as
    GeneMark and AUGUSTUS: every genetic code defines a different amino-acid
    alphabet and therefore a different pattern of synonymous codons.  The code
    whose synonymous-group structure best fits the observed codon frequencies
    inside the detected ORFs is the most likely code in use.

    A secondary signal — the stop-boundary score — checks what fraction of
    ORF-terminal codons are valid stop codons under each code, providing a
    sanity check that is independent of interior codon usage.

Genetic codes evaluated (NCBI table IDs, matching the screenshots):
    1  Standard
    2  Vertebrate Mitochondrial
    3  Yeast Mitochondrial
    4  Mold, Protozoan, and Coelenterate Mitochondrial / Mycoplasma
    6  Ciliate, Dasycladacean and Hexamita Nuclear
    9  Echinoderm and Flatworm Mitochondrial

Public API
----------
infer_genetic_code          Score all candidate codes and return ranked results.
print_inference_report      Print a formatted ranking table to stdout.
write_inference_report      Write the full report to a text file.

Integration
-----------
Call ``calculate_orf_stats()`` from orf_analysis.py *before* calling
``infer_genetic_code()`` so that each ORF dict already has a ``"sequence"``
key.  If that key is absent, sequences are re-extracted automatically from
the forward-strand DNA.
"""

from __future__ import annotations

import math
import os
from collections import Counter
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from Bio.Data import CodonTable

from src.orf_finder_lib.frame_scanner import extract_orf_sequence

# ---------------------------------------------------------------------------
# Candidate genetic codes
# ---------------------------------------------------------------------------

#: NCBI table IDs and human-readable names for all evaluated codes.
CANDIDATE_CODES: Dict[int, str] = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold / Protozoan / Coelenterate Mitochondrial",
    6: "Ciliate / Dasycladacean / Hexamita Nuclear",
    9: "Echinoderm and Flatworm Mitochondrial",
}

# Pseudocount added to every reference frequency before taking log.
# Prevents log(0) for codons absent from the reference table.
_PSEUDOCOUNT: float = 1e-10

# Report formatting width — matches statistics_summary.py
_W: int = 72


# ---------------------------------------------------------------------------
# Formatting helpers (mirrors statistics_summary.py style)
# ---------------------------------------------------------------------------

def _rule(char: str = "─") -> str:
    return char * _W


def _header(title: str, char: str = "═") -> str:
    pad = (_W - len(title) - 2) // 2
    return f"{char * pad} {title} {char * (_W - pad - len(title) - 2)}"


def _subheader(title: str) -> str:
    return f"  {title}\n  {'─' * len(title)}"


# ---------------------------------------------------------------------------
# Codon-table helpers
# ---------------------------------------------------------------------------

def _get_table(table_id: int) -> CodonTable.CodonTable:
    """Return the BioPython unambiguous-DNA CodonTable for *table_id*."""
    return CodonTable.unambiguous_dna_by_id[table_id]


def _build_reference_frequencies(table: CodonTable.CodonTable) -> Dict[str, float]:
    """
    Build a reference codon frequency distribution from a genetic code table.

    Uses an equal-usage-within-synonymous-groups model: each amino acid
    receives equal weight (1 / n_amino_acids), and that weight is split
    evenly among all codons that encode it.  This is a standard null model
    used when organism-specific codon usage data are unavailable.

    Stop codons receive a near-zero probability (_PSEUDOCOUNT) because they
    should not appear inside ORF bodies — the scoring is applied only to the
    interior codons of each ORF, not to its terminal stop.

    Parameters
    ----------
    table : CodonTable.CodonTable
        BioPython codon table object for a specific genetic code.

    Returns
    -------
    Dict[str, float]
        Mapping from three-letter DNA codon to its reference frequency.
        All 64 sense codons have positive frequencies; stop codons have
        frequency == _PSEUDOCOUNT.
    """
    # Group sense codons by amino acid
    aa_groups: Dict[str, List[str]] = {}
    for codon, aa in table.forward_table.items():
        aa_groups.setdefault(aa, []).append(codon)

    n_amino_acids = len(aa_groups)
    frequencies: Dict[str, float] = {}

    for aa, codons in aa_groups.items():
        # Equal share per amino acid, equal share within the synonymous group
        per_codon = 1.0 / (n_amino_acids * len(codons))
        for codon in codons:
            frequencies[codon] = per_codon

    # Stop codons: near-zero — they should not appear in ORF bodies
    for codon in table.stop_codons:
        frequencies[codon] = _PSEUDOCOUNT

    return frequencies


# ---------------------------------------------------------------------------
# Codon extraction
# ---------------------------------------------------------------------------

def _extract_body_codons(
    flat_list: List[Dict[str, Any]],
    dna_sequence: str,
) -> Counter:
    """
    Collect all interior (non-start, non-stop) codons from the detected ORFs.

    The start codon (first 3 nt) and stop codon (last 3 nt) are excluded
    because they are the same under most codes and would dilute the signal
    that distinguishes codes from one another.  Codons containing 'N' are
    skipped to avoid noise from ambiguous bases.

    The ``"sequence"`` key is used when present (set by calculate_orf_stats);
    otherwise the sequence is re-extracted from *dna_sequence*.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Flat ORF list from find_orfs(), optionally enriched by
        calculate_orf_stats().
    dna_sequence : str
        Forward-strand DNA sequence.

    Returns
    -------
    Counter
        Codon → occurrence count for all interior codons across all ORFs.
    """
    counts: Counter = Counter()

    for orf in flat_list:
        seq: str = orf.get("sequence") or extract_orf_sequence(orf, dna_sequence)
        # Need at least start + 1 interior codon + stop = 9 nt
        if len(seq) < 9:
            continue
        body = seq[3:-3]   # strip start and stop codons
        for i in range(0, len(body) - 2, 3):
            codon = body[i : i + 3]
            if len(codon) == 3 and "N" not in codon:
                counts[codon] += 1

    return counts


# ---------------------------------------------------------------------------
# Scoring functions
# ---------------------------------------------------------------------------

def _log_likelihood(
    observed: Counter,
    reference: Dict[str, float],
) -> float:
    """
    Compute the total log-likelihood of the observed codon counts under the
    reference frequency distribution.

    Uses the standard multinomial log-likelihood:
        L = Σ  count(c) * log( P(c | code) + pseudocount )

    The pseudocount inside the log prevents -inf for codons observed in the
    data but absent from the reference table (e.g. codons that are stops
    under the candidate code but happened to appear in the ORF body).

    Parameters
    ----------
    observed : Counter
        Codon → count from _extract_body_codons().
    reference : Dict[str, float]
        Codon → probability from _build_reference_frequencies().

    Returns
    -------
    float
        Total log-likelihood.  More negative == worse fit.
        Returns -inf if no codons were observed.
    """
    if not observed:
        return float("-inf")

    ll = 0.0
    for codon, count in observed.items():
        freq = reference.get(codon, _PSEUDOCOUNT)
        ll += count * math.log(freq + _PSEUDOCOUNT)
    return ll


def _stop_boundary_score(
    flat_list: List[Dict[str, Any]],
    dna_sequence: str,
    table: CodonTable.CodonTable,
) -> float:
    """
    Measure how well a genetic code's stop codons align with ORF boundaries.

    Reads the last three nucleotides of every ORF sequence and checks whether
    that terminal codon is a recognised stop codon under *table*.  Returns the
    fraction of ORFs whose terminal codon is a valid stop under this code.

    A perfect score of 1.0 means every ORF ends with a stop codon that this
    code recognises — strong evidence the code is correct.  A low score means
    many ORFs are "terminated" by codons that are not stops under this code,
    suggesting the code is wrong.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Flat ORF list from find_orfs().
    dna_sequence : str
        Forward-strand DNA sequence.
    table : CodonTable.CodonTable
        BioPython codon table for the candidate code.

    Returns
    -------
    float
        Fraction of ORFs with a valid terminal stop codon, in [0.0, 1.0].
        Returns 0.0 when flat_list is empty.
    """
    if not flat_list:
        return 0.0

    stop_set = set(table.stop_codons)
    valid = 0

    for orf in flat_list:
        seq: str = orf.get("sequence") or extract_orf_sequence(orf, dna_sequence)
        if len(seq) >= 6:
            terminal = seq[-3:]
            if terminal in stop_set:
                valid += 1

    return valid / len(flat_list)


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def infer_genetic_code(
    flat_list: List[Dict[str, Any]],
    dna_sequence: str,
    candidate_ids: Optional[List[int]] = None,
) -> List[Dict[str, Any]]:
    """
    Score each candidate genetic code against the ORFs detected in the
    ORCA pipeline and return a ranked list of results.

    Scoring method
    --------------
    For each candidate code:

    1. Build a reference codon frequency distribution using the equal-usage-
       within-synonymous-groups model (see _build_reference_frequencies).
    2. Extract interior codon counts from all ORF bodies
       (see _extract_body_codons).
    3. Compute the log-likelihood per codon (normalised so sequences of
       different total lengths are directly comparable).
    4. Compute the stop-boundary score (fraction of ORFs with a valid
       terminal stop codon under this code).
    5. Combine: combined_score = 0.8 * ll_per_codon + 0.2 * boundary_score.
       The 0.8 / 0.2 weighting reflects that interior codon usage carries
       ~4× more signal than stop-boundary agreement, because it uses all 61
       sense codons rather than just 3 stop positions.

    The returned list is sorted best-first by combined_score.

    Parameters
    ----------
    flat_list : List[Dict[str, Any]]
        Flat ORF list from find_orfs().  If calculate_orf_stats() has already
        been called, the ``"sequence"`` key is reused; otherwise sequences are
        re-extracted automatically.
    dna_sequence : str
        Forward-strand DNA sequence used to find the ORFs.
    candidate_ids : List[int], optional
        NCBI genetic code table IDs to evaluate.  Defaults to all six codes
        in CANDIDATE_CODES.

    Returns
    -------
    List[Dict[str, Any]]
        Ranked list of result dicts, best code first.  Each dict contains:

        ``code_id``            -- NCBI table ID (int)
        ``code_name``          -- Human-readable name (str)
        ``stop_codons``        -- Stop codons under this code (List[str])
        ``n_body_codons``      -- Total interior codons scored (int)
        ``log_likelihood``     -- Raw total log-likelihood (float)
        ``ll_per_codon``       -- Log-likelihood normalised per codon (float)
        ``stop_boundary_score``-- Fraction of ORFs with valid terminal stop (float)
        ``combined_score``     -- Weighted combination used for ranking (float)

    Notes
    -----
    Returns an empty list and prints a warning if no interior codons can be
    extracted (e.g. all ORFs are shorter than 9 nt).
    """
    if candidate_ids is None:
        candidate_ids = list(CANDIDATE_CODES.keys())

    observed: Counter = _extract_body_codons(flat_list, dna_sequence)
    n_body_codons: int = sum(observed.values())

    if n_body_codons == 0:
        print(
            "[WARNING] genetic_code_inference: no interior codons could be "
            "extracted from the detected ORFs.  All ORFs may be shorter than "
            "9 nt.  Genetic code inference requires longer ORFs."
        )
        return []

    results: List[Dict[str, Any]] = []

    for code_id in candidate_ids:
        table    = _get_table(code_id)
        ref      = _build_reference_frequencies(table)
        ll       = _log_likelihood(observed, ref)
        ll_norm  = ll / n_body_codons
        boundary = _stop_boundary_score(flat_list, dna_sequence, table)

        # Weighted combination: LL dominates; boundary breaks ties
        combined = 0.8 * ll_norm + 0.2 * boundary

        results.append({
            "code_id":             code_id,
            "code_name":           CANDIDATE_CODES.get(code_id, f"Code {code_id}"),
            "stop_codons":         sorted(table.stop_codons),
            "n_body_codons":       n_body_codons,
            "log_likelihood":      ll,
            "ll_per_codon":        ll_norm,
            "stop_boundary_score": boundary,
            "combined_score":      combined,
        })

    results.sort(key=lambda r: r["combined_score"], reverse=True)
    return results


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

def _score_delta(results: List[Dict[str, Any]]) -> Optional[float]:
    """Return the combined-score gap between the top-two candidates."""
    if len(results) < 2:
        return None
    return results[0]["combined_score"] - results[1]["combined_score"]


def _interpret_confidence(delta: Optional[float]) -> str:
    """Translate the score gap into a human-readable confidence label."""
    if delta is None:
        return "N/A (only one code evaluated)"
    if delta > 0.05:
        return "high — best code scores clearly above the next candidate"
    if delta > 0.01:
        return "moderate — best code leads, but alternatives are close"
    return "low — multiple codes are nearly equally supported; more data needed"


def print_inference_report(results: List[Dict[str, Any]]) -> None:
    """
    Print a formatted genetic-code ranking table to stdout.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        Ranked output from infer_genetic_code(), best code first.
    """
    if not results:
        print("[INFO] No genetic code inference results to display.")
        return

    delta = _score_delta(results)

    print(f"\n{_header('GENETIC CODE INFERENCE REPORT')}")
    print(f"  Interior codons scored : {results[0]['n_body_codons']:,}")
    print(f"  Codes evaluated        : {len(results)}")
    print()

    # Column header
    print(
        f"  {'Rank':<5} {'ID':>3}  {'Genetic Code':<44}"
        f"{'LL/codon':>10} {'Boundary':>10} {'Score':>8}"
    )
    print(f"  {_rule('─')}")

    for rank, r in enumerate(results, 1):
        marker = "  ◄ BEST FIT" if rank == 1 else ""
        print(
            f"  {rank:<5} {r['code_id']:>3}  {r['code_name']:<44}"
            f"{r['ll_per_codon']:>10.4f}"
            f"{r['stop_boundary_score']:>10.4f}"
            f"{r['combined_score']:>8.4f}"
            f"{marker}"
        )

    best = results[0]
    print(f"\n{_rule()}")
    print(f"\n  Best fit  : NCBI Genetic Code {best['code_id']} — {best['code_name']}")
    print(f"  Stop codons under this code : {', '.join(best['stop_codons'])}")
    print(f"  Log-likelihood / codon      : {best['ll_per_codon']:.4f}")
    print(f"  Stop-boundary score         : {best['stop_boundary_score']:.4f}  "
          f"({best['stop_boundary_score'] * 100:.1f}% of ORFs end on a valid stop)")
    print(f"  Combined score              : {best['combined_score']:.4f}")
    print(f"  Confidence                  : {_interpret_confidence(delta)}")
    print(f"\n{_rule()}")


def write_inference_report(
    results: List[Dict[str, Any]],
    filename: str = "output/genetic_code_report.txt",
    accession: str = "N/A",
) -> None:
    """
    Write the genetic-code inference report to a text file.

    The layout mirrors the style used by write_stats_to_file() in
    statistics_summary.py: run metadata at the top, then a ranked table,
    then a detailed section for the best-fit code.

    Parameters
    ----------
    results : List[Dict[str, Any]]
        Ranked output from infer_genetic_code(), best code first.
    filename : str
        Destination file path.  Parent directory is created if absent.
    accession : str
        Accession number or label for the sequence being reported on.
    """
    os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)

    if not results:
        with open(filename, "w") as fh:
            fh.write(_header("GENETIC CODE INFERENCE REPORT") + "\n")
            fh.write("  No results — no interior codons could be extracted.\n")
        return

    delta = _score_delta(results)
    best  = results[0]

    with open(filename, "w") as fh:
        fh.write(_header("GENETIC CODE INFERENCE REPORT") + "\n")

        # --- Run metadata ---
        fh.write(_subheader("Run Parameters") + "\n")
        fh.write(f"  Generated              : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"  Sequence               : {accession}\n")
        fh.write(f"  Interior codons scored : {best['n_body_codons']:,}\n")
        fh.write(f"  Codes evaluated        : {len(results)}\n")
        fh.write(f"  Scoring method         : log-likelihood (equal usage within "
                 f"synonymous groups) + stop-boundary fraction\n")
        fh.write(f"  Score formula          : 0.8 × LL/codon + 0.2 × boundary\n\n")

        # --- Ranking table ---
        fh.write(_subheader("Ranked Genetic Codes") + "\n")
        fh.write(
            f"  {'Rank':<5} {'ID':>3}  {'Genetic Code':<44}"
            f"{'LL/codon':>10} {'Boundary':>10} {'Score':>8}\n"
        )
        fh.write(f"  {_rule('─')}\n")

        for rank, r in enumerate(results, 1):
            marker = "  <-- BEST FIT" if rank == 1 else ""
            fh.write(
                f"  {rank:<5} {r['code_id']:>3}  {r['code_name']:<44}"
                f"{r['ll_per_codon']:>10.4f}"
                f"{r['stop_boundary_score']:>10.4f}"
                f"{r['combined_score']:>8.4f}"
                f"{marker}\n"
            )

        # --- Best-fit details ---
        fh.write(f"\n{_rule()}\n")
        fh.write(_subheader("Best-Fit Code Details") + "\n")
        fh.write(f"  NCBI code ID        : {best['code_id']}\n")
        fh.write(f"  Name                : {best['code_name']}\n")
        fh.write(f"  Stop codons         : {', '.join(best['stop_codons'])}\n")
        fh.write(f"  Log-likelihood      : {best['log_likelihood']:.2f}\n")
        fh.write(f"  LL per codon        : {best['ll_per_codon']:.4f}\n")
        fh.write(f"  Stop-boundary score : {best['stop_boundary_score']:.4f}  "
                 f"({best['stop_boundary_score'] * 100:.1f}% of ORFs end on a valid stop)\n")
        fh.write(f"  Combined score      : {best['combined_score']:.4f}\n")
        fh.write(f"  Confidence          : {_interpret_confidence(delta)}\n\n")

        # --- Runner-up comparison ---
        if len(results) >= 2:
            second = results[1]
            fh.write(_subheader("Runner-up Comparison") + "\n")
            fh.write(f"  Runner-up code      : {second['code_id']} — {second['code_name']}\n")
            fh.write(f"  Score gap           : {delta:.4f}\n")
            fh.write(f"  Interpretation      : {_interpret_confidence(delta)}\n\n")

        # --- Score interpretation guide ---
        fh.write(_subheader("Score Interpretation") + "\n")
        fh.write(
            "  LL/codon  : log-likelihood per interior codon.  Higher (less negative)\n"
            "              means the observed codon usage fits the code's synonymous-\n"
            "              group structure better.\n\n"
            "  Boundary  : fraction of ORFs whose terminal codon is a stop codon\n"
            "              under this code.  1.0 = all ORFs terminate correctly.\n\n"
            "  Score gap > 0.05 → high confidence.\n"
            "  Score gap 0.01–0.05 → moderate confidence.\n"
            "  Score gap < 0.01 → low confidence; sequence may be too short or\n"
            "              the codon usage may not be informative enough to\n"
            "              discriminate between codes.\n"
        )
