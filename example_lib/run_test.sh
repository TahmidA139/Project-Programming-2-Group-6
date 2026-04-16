#!/usr/bin/env bash
# =============================================================================
# run_test.sh – ORCA quick-start test
# work done by amanda 
# Tests the two most common use cases using two CCR5 mRNA variants from NCBI.
# Sequences are fetched automatically — no local files needed.
#
# Usage:
#   bash example_lib/run_test.sh youremail@gmail.com
#
# Requirements:
#   conda activate ORCA        (see README for setup instructions)
#   Internet connection        (sequences are fetched live from NCBI)
# =============================================================================

EMAIL="$1"

if [[ -z "$EMAIL" ]]; then
    echo "Please provide your email: bash run_test.sh your@email.com"
    exit 1
fi

# --------------------------------------------------------------------------- #
# Test 1: Single sequence — default settings
# Fetches CCR5 variant 1, scans all 6 reading frames, ATG start codons only,
# minimum ORF length 30 nt, nested ORFs included.
# --------------------------------------------------------------------------- #
echo "Running Test 1: single sequence, default settings..."

python -m src.main \
    --accession NM_001838.4 \
    --email "$EMAIL"

# --------------------------------------------------------------------------- #
# Test 2: Comparative mode — two CCR5 variants side by side
# Runs the full pipeline on both sequences and generates a side-by-side
# comparison of ORF structure and codon usage.
# --------------------------------------------------------------------------- #
echo "Running Test 2: comparative mode..."

python -m src.main \
    --accession  NM_001838.4 \
    --accession2 NM_001301717.2 \
    --email "$EMAIL"

echo ""
echo "Done. Check the output/ folder for results."
