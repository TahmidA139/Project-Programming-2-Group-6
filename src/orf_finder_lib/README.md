# orf_finder_lib — Module Reference

Detects Open Reading Frames (ORFs) across all six reading frames of a DNA sequence. Supports non-canonical start codons, nested ORF detection, and reverse-complement scanning.

Library setup:
```
src/orf_finder_lib/
├── frame_scanner.py   # Low-level sequence utilities and per-frame scanning

Forward DNA string
       │
       ├──► scan frames 0,1,2 on + strand ──────────────┐
       │                                                  │
       └──► reverse_complement ──► scan frames 0,1,2    ──┤
                                   on - strand            │
                                                          ▼
                                              collect all ORF dicts
                                                          │
                                                    _mark_nested()
                                                          │
                                                   final ORF list

├── orf_finder.py      # orchestrates all six frames and builds output
└── output_writer.py   # Summary printing and CSV export
```

---

### Outputs

The outputs for this module are a terminal output and a csv file. The CSV file defaults to output/orfs.csv if --output option is not used.

**For example the command:**
```
python -m src.main --accession NM_001301717 --email edecocke@charlotte.edu --min-length 75 --ignore-nested
```
**Would output this to the terminal:**
```
[ORCA] Processing sequence 1: NM_001301717
[VALIDATION] Sequence is valid — 2191 bp ready for analysis.

----------  ORF Summary — Sequence 1 ----------
  Total ORFs found            : 9
  Forward strand (+)          : 3
  Reverse strand (-)          : 6
  Canonical   (ATG)           : 9
  Nested ORFs detected        : 0
----------------------------------------------
```
**CSV file**

Each sequence block starts with the accession on its own row, followed by column headers and one row per ORF. In comparative mode the two blocks are separated by two blank rows. The `sequence (5'->3')` column always reads in the 5′→3′ direction regardless of strand.

```
NM_001301717
orf_id  strand  start_codon  frame  start  end   length_nt  sequence (5'->3')
ORF1    +       ATG          0      99     1218  1119       ATGAAAAGCGTGCTGGTGGTG...
ORF2    +       ATG          0      1287   1476  189        ATGACTCAGGACATCCCCCCG...
ORF3    -       ATG          1      1842   1935  93         ATGTCATCCCCACTCTGGAGC...
ORF4    +       ATG          2      1583   2183  600        ATGAACCTTCTGGCCTCCCAC...
ORF5    -       ATG          2      1013   1178  165        ATGGAGGAGCGCCGGATGTGC...
ORF6    -       ATG          2      914    998   84         ATGTTGAGTTGCTTACTGAGC...
ORF7    -       ATG          2      806    899   93         ATGAAGACCACGACCACAGCG...
ORF8    -       ATG          2      572    797   225        ATGGCCAGCAGGGGGACCAGA...
ORF9    -       ATG          2      200    284   84         ATGATGGAGTACATGATAGGG...
```
