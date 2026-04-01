#!/usr/bin/env python3
"""
orf_analysis.py
Amandas portion 

This module finds repeated ORFs and calculates simple similarity scores. (Not anymore because we took out similirity scores)
"""
#!/usr/bin/env python3
"""
ORF_analysis.py

Handles ORF extraction and per-ORF statistics.
"""

# -------------------------------
# REQUIRED FUNCTIONS
# -------------------------------

def find_repeated_orfs(orfs):
    """
    Input: list of ORF dicts
    Output: dict of repeated sequences
    """
    counts = {}

    for orf in orfs:
        seq = orf["sequence"]
        counts[seq] = counts.get(seq, 0) + 1

    return {seq: count for seq, count in counts.items() if count > 1}


def extract_sequence(dna, start, end, strand="+"):
    """
    Extract sequence from DNA.
    Handles reverse strand.
    """
    seq = dna[start:end]

    if strand == "-":
        complement = str.maketrans("ATGC", "TACG")
        seq = seq.translate(complement)[::-1]

    return seq


def gc_content(sequence):
    """
    % GC content
    """
    if not sequence:
        return 0

    gc = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100


def protein_length(sequence):
    """
    Number of amino acids
    """
    return len(sequence) // 3


def calculate_orf_stats(orfs, dna_sequence):
    """
    Adds sequence, GC content, and protein length to each ORF
    """
    for orf in orfs:
        start = orf["start"]
        end = orf["end"]
        strand = orf.get("strand", "+")

        seq = extract_sequence(dna_sequence, start, end, strand)

        orf["sequence"] = seq
        orf["gc_content"] = gc_content(seq)
        orf["protein_length"] = protein_length(seq)

    return orfs

import csv
# -------------------------------
# Helper: Codon Usage
# -------------------------------

def codon_usage(sequence):
    counts = {}

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        counts[codon] = counts.get(codon, 0) + 1

    return counts


# -------------------------------
# 1. WRITE SUMMARY FILE
# -------------------------------

def write_stats_to_file(orfs, filename="orf_summary.txt"):
    total_orfs = len(orfs)

    avg_gc = sum(o["gc_content"] for o in orfs) / total_orfs if total_orfs else 0
    longest = max(orfs, key=lambda x: len(x["sequence"]), default=None)

    with open(filename, "w") as f:
        f.write("=== ORF SUMMARY REPORT ===\n\n")

        # Dataset stats
        f.write(f"Total ORFs: {total_orfs}\n")
        f.write(f"Average GC Content: {avg_gc:.2f}%\n")

        if longest:
            f.write("\nLongest ORF:\n")
            f.write(f"Length: {len(longest['sequence'])}\n")
            f.write(f"GC Content: {longest['gc_content']:.2f}%\n")

        # Table
        f.write("\n--- Per ORF Stats ---\n")
        f.write(f"{'Index':<6}{'Length':<10}{'GC%':<10}{'Protein Len':<12}\n")

        for i, orf in enumerate(orfs):
            f.write(f"{i:<6}{len(orf['sequence']):<10}{orf['gc_content']:<10.2f}{orf['protein_length']:<12}\n")

        # Codon usage
        f.write("\n--- Codon Usage ---\n")
        total_codons = {}

        for orf in orfs:
            usage = codon_usage(orf["sequence"])
            for codon, count in usage.items():
                total_codons[codon] = total_codons.get(codon, 0) + count

        for codon, count in sorted(total_codons.items()):
            f.write(f"{codon}: {count}\n")


# -------------------------------
# 2. COMPARATIVE REPORT (TEXT)
# -------------------------------

def write_comparative_report(orfs1, orfs2, filename="comparison.txt"):
    def avg_gc(orfs):
        return sum(o["gc_content"] for o in orfs) / len(orfs) if orfs else 0

    with open(filename, "w") as f:
        f.write("=== COMPARATIVE REPORT ===\n\n")

        f.write("Dataset 1:\n")
        f.write(f"ORFs: {len(orfs1)}\n")
        f.write(f"Avg GC: {avg_gc(orfs1):.2f}%\n\n")

        f.write("Dataset 2:\n")
        f.write(f"ORFs: {len(orfs2)}\n")
        f.write(f"Avg GC: {avg_gc(orfs2):.2f}%\n\n")


# -------------------------------
# 3. COMPARATIVE CSV (CODON DELTA)
# -------------------------------

def write_comparative_csv(orfs1, orfs2, filename="codon_comparison.csv"):
    def total_codons(orfs):
        totals = {}
        for orf in orfs:
            usage = codon_usage(orf["sequence"])
            for codon, count in usage.items():
                totals[codon] = totals.get(codon, 0) + count
        return totals

    codons1 = total_codons(orfs1)
    codons2 = total_codons(orfs2)

    all_codons = set(codons1.keys()).union(codons2.keys())

    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Codon", "Seq1", "Seq2", "Delta"])

        for codon in sorted(all_codons):
            c1 = codons1.get(codon, 0)
            c2 = codons2.get(codon, 0)
            writer.writerow([codon, c1, c2, c1 - c2])
            
"""
#Possible pipline when the other parts are completed :)
from orf_finder import find_all_orfs
from ORF_analysis import calculate_orf_stats, find_repeated_orfs
from orf_stats import write_stats_to_file, write_comparative_report, write_comparative_csv


def main():
    dna = "ATGAAATAGATGCCCTAAATGAAATGA"

    # Find ORFs
    orfs = find_all_orfs(dna)

    # Add stats
    orfs = calculate_orf_stats(orfs, dna)

    # Repeats
    repeats = find_repeated_orfs(orfs)

    print("Repeated ORFs:", repeats)

    # Write outputs
    write_stats_to_file(orfs)
    write_comparative_report(orfs, orfs)
    write_comparative_csv(orfs, orfs)


if __name__ == "__main__":
    main()
"""

