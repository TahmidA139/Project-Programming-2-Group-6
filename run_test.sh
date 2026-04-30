#!/usr/bin/env bash

# ==============================================================================
# run_test.sh
# Test script for ORCA: ORF Recognition and Comparative Analysis
#
# This script verifies that the project runs correctly by executing src/main.py
# with the provided example FASTA files and checking that all expected output
# files are produced. It also performs basic content checks on the results.
#
# Usage:
#   bash run_test.sh
#
# The script assumes that:
#   - The conda environment "ORCA" is activated with the required dependencies
#     (see README for setup instructions).
#   - The example input files are in example_input_files/.
#   - src/main.py is in the src/ directory under the project root.
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

# Directory where this script lives (also the project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input FASTA files (all live in example_input_files/)
OR2B6_FASTA="${SCRIPT_DIR}/example_input_files/OR2B6_sequence.fasta"
IUPAC_FASTA="${SCRIPT_DIR}/example_input_files/IUPAC_ambiguity_test.fasta"
MULTI_FASTA="${SCRIPT_DIR}/example_input_files/multiple_sequence_file.fasta"

TEST_OUTDIR="${SCRIPT_DIR}/test_output"

# Dummy email — required by the --email flag to suppress the interactive
# prompt; NCBI is never contacted when --fasta is used.
DUMMY_EMAIL="test@example.com"

# Counter for passed and failed checks
PASSED=0
FAILED=0


# ---- Helper functions --------------------------------------------------------

print_header() {
    # Print a section header to make the output easier to read.
    # Arguments:
    #   $1  The header text.
    echo ""
    echo "========================================"
    echo "  $1"
    echo "========================================"
    echo ""
}

check_file_exists() {
    # Verify that a file exists and is not empty.
    # Arguments:
    #   $1  Path to the file.
    #   $2  A short description of what the file should contain.
    if [[ -s "$1" ]]; then
        echo "  PASS: $2 exists and is not empty."
        PASSED=$((PASSED + 1))
    else
        echo "  FAIL: $2 is missing or empty ($1)."
        FAILED=$((FAILED + 1))
    fi
}

check_string_in_file() {
    # Verify that a file contains a specific string.
    # Arguments:
    #   $1  Path to the file.
    #   $2  The string to search for.
    #   $3  A short description of the check.
    if grep -q "$2" "$1" 2>/dev/null; then
        echo "  PASS: $3"
        PASSED=$((PASSED + 1))
    else
        echo "  FAIL: $3 (string '$2' not found in $1)."
        FAILED=$((FAILED + 1))
    fi
}

check_csv_rows() {
    # Verify that a CSV file has at least a given number of data rows
    # (excluding the header).
    # Arguments:
    #   $1  Path to the CSV file.
    #   $2  Minimum number of data rows expected.
    #   $3  A short description of the check.
    local row_count
    # 'tail -n +2' skips the header; 'wc -l' counts remaining lines
    row_count=$(tail -n +2 "$1" 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$row_count" -ge "$2" ]]; then
        echo "  PASS: $3 (found ${row_count} rows)."
        PASSED=$((PASSED + 1))
    else
        echo "  FAIL: $3 (expected at least $2 rows, found ${row_count})."
        FAILED=$((FAILED + 1))
    fi
}

check_exit_nonzero() {
    # Verify that a command exited with a non-zero code (i.e. the expected
    # failure occurred).
    # Arguments:
    #   $1  The exit code to check.
    #   $2  A short description of the check.
    if [[ "$1" -ne 0 ]]; then
        echo "  PASS: $2 (exited with code ${1}, as expected)."
        PASSED=$((PASSED + 1))
    else
        echo "  FAIL: $2 (expected non-zero exit code, but got 0)."
        FAILED=$((FAILED + 1))
    fi
}


# ---- Pre-flight checks -------------------------------------------------------

preflight_checks() {
    # Make sure the input files and source module exist before we try to run
    # anything. Runs silently on success; prints errors and exits on failure.

    local all_ok=true

    if [[ ! -f "${SCRIPT_DIR}/src/main.py" ]]; then
        echo "  ERROR: src/main.py not found in ${SCRIPT_DIR}."
        all_ok=false
    fi

    if [[ ! -f "$OR2B6_FASTA" ]]; then
        echo "  ERROR: Input file not found: ${OR2B6_FASTA}"
        all_ok=false
    fi

    if [[ ! -f "$IUPAC_FASTA" ]]; then
        echo "  ERROR: Input file not found: ${IUPAC_FASTA}"
        all_ok=false
    fi

    if [[ ! -f "$MULTI_FASTA" ]]; then
        echo "  ERROR: Input file not found: ${MULTI_FASTA}"
        all_ok=false
    fi

    if ! python -c "import matplotlib, numpy" 2>/dev/null; then
        echo "  ERROR: Python cannot import matplotlib or numpy."
        echo "         Make sure the ORCA conda environment is activated."
        echo "         See README for setup instructions."
        all_ok=false
    fi

    if [[ "$all_ok" == false ]]; then
        echo ""
        echo "Pre-flight checks failed. Exiting."
        exit 1
    fi
}


# ---- Run the program ---------------------------------------------------------

run_tests() {
    # Execute src/main.py for each test case and preserve the outputs.
    # All program output is redirected to /dev/null — only the check results
    # printed by verify_output_files and verify_output_content are shown.
    # Each test passes --outdir directly to ORCA so results land in the
    # correct test_output/testN/ subdirectory without any copying.
    # Test 3 also saves its stdout and stderr as example_stdout.txt and
    # example_stderr.txt in test_output/ for use as reference files.

    # Start with a clean test_output/ directory
    rm -rf "$TEST_OUTDIR"
    mkdir -p "$TEST_OUTDIR"

    # ---------------------------------------------------------------------- #
    # Test 1: Single sequence — OR2B6 (NM_012367.1), default settings
    # Verifies the basic single-sequence code path end-to-end.
    # ---------------------------------------------------------------------- #
    python -m src.main \
        --fasta  "$OR2B6_FASTA" \
        --email  "$DUMMY_EMAIL" \
        --outdir "${TEST_OUTDIR}/test1" \
        > /dev/null 2>&1

    # ---------------------------------------------------------------------- #
    # Test 2: IUPAC ambiguity — CCR7 (NM_001838.4)
    # The sequence contains non-standard IUPAC bases (V, S, D, X, Z, Q, P, W).
    # Verifies that the validation/cleaning step handles ambiguity codes
    # gracefully and the pipeline still produces output.
    # ---------------------------------------------------------------------- #
    python -m src.main \
        --fasta  "$IUPAC_FASTA" \
        --email  "$DUMMY_EMAIL" \
        --outdir "${TEST_OUTDIR}/test2" \
        > /dev/null 2>&1

    # ---------------------------------------------------------------------- #
    # Test 3: Comparative mode — OR2B6 vs CCR7
    # Runs the full comparative pipeline on two local FASTA files and checks
    # that all comparative output files are produced: two per-sequence cleaned
    # FASTAs, a combined ORF CSV, the ORF map, and the comparison report.
    # ---------------------------------------------------------------------- #
    # stdout and stderr from this run are saved as reference files so readers
    # can see what a correct comparative run looks like without running it.
    python -m src.main \
        --fasta  "$OR2B6_FASTA" \
        --fasta2 "$IUPAC_FASTA" \
        --email  "$DUMMY_EMAIL" \
        --outdir "${TEST_OUTDIR}/test3" \
        > "${TEST_OUTDIR}/example_stdout.txt" \
        2> "${TEST_OUTDIR}/example_stderr.txt"

    # ---------------------------------------------------------------------- #
    # Test 4: Expected failure — multi-sequence FASTA
    # main.py requires exactly one sequence per local file and must exit with
    # a non-zero code when given multiple_sequence_file.fasta (17 sequences).
    # This verifies that the input-validation error path works correctly.
    # ---------------------------------------------------------------------- #
    python -m src.main \
        --fasta  "$MULTI_FASTA" \
        --email  "$DUMMY_EMAIL" \
        > /dev/null 2>&1
    local exit4=$?
    check_exit_nonzero $exit4 \
        "Pipeline correctly rejects a FASTA file containing multiple sequences."
}


# ---- Verify output files -----------------------------------------------------

verify_output_files() {
    # Check that every expected output file was created and is non-empty.

    # ---- Test 1: single-sequence outputs --------------------------------- #
    # In single-sequence mode the cleaned FASTA is written as
    # cleaned_sequence_1.fasta (seq_num=1, comparative=False).
    check_file_exists "${TEST_OUTDIR}/test1/orfs.csv"                  "Test 1 ORF CSV"
    check_file_exists "${TEST_OUTDIR}/test1/cleaned_sequence_1.fasta"  "Test 1 cleaned FASTA"
    check_file_exists "${TEST_OUTDIR}/test1/orf_map.png"               "Test 1 ORF map image"
    check_file_exists "${TEST_OUTDIR}/test1/orf_summary.txt"           "Test 1 ORF summary"

    # ---- Test 2: IUPAC ambiguity — same output set as single-sequence mode #
    check_file_exists "${TEST_OUTDIR}/test2/orfs.csv"                  "Test 2 ORF CSV"
    check_file_exists "${TEST_OUTDIR}/test2/cleaned_sequence_1.fasta"  "Test 2 cleaned FASTA"
    check_file_exists "${TEST_OUTDIR}/test2/orf_map.png"               "Test 2 ORF map image"
    check_file_exists "${TEST_OUTDIR}/test2/orf_summary.txt"           "Test 2 ORF summary"

    # ---- Test 3: comparative mode ---------------------------------------- #
    # In comparative mode each sequence gets its own cleaned FASTA prefixed
    # with "comp_" and numbered by position (seq_num 1 and 2).
    # Both sequences share a single combined orfs.csv and a single orf_map.png.
    # The comparison report is written as orf_comparison_report.txt.
    check_file_exists "${TEST_OUTDIR}/test3/orfs.csv"                       "Test 3 combined ORF CSV"
    check_file_exists "${TEST_OUTDIR}/test3/comp_cleaned_sequence_1.fasta"  "Test 3 sequence-1 cleaned FASTA"
    check_file_exists "${TEST_OUTDIR}/test3/comp_cleaned_sequence_2.fasta"  "Test 3 sequence-2 cleaned FASTA"
    check_file_exists "${TEST_OUTDIR}/test3/orf_map.png"                    "Test 3 comparative ORF map"
    check_file_exists "${TEST_OUTDIR}/test3/orf_comparison_report.txt"      "Test 3 comparison report"
}


# ---- Verify output content ---------------------------------------------------

verify_output_content() {
    # Perform basic sanity checks on the content of the text and CSV outputs.

    # ---- Test 1: OR2B6 (accession NM_012367.1 is read from the FASTA header) #

    check_string_in_file "${TEST_OUTDIR}/test1/orfs.csv" \
        "NM_012367" \
        "Test 1 ORF CSV references the OR2B6 accession (NM_012367.1)."

    check_csv_rows "${TEST_OUTDIR}/test1/orfs.csv" \
        1 \
        "Test 1 ORF CSV contains at least 1 ORF."

    check_string_in_file "${TEST_OUTDIR}/test1/cleaned_sequence_1.fasta" \
        "NM_012367" \
        "Test 1 cleaned FASTA header contains the OR2B6 accession."

    check_string_in_file "${TEST_OUTDIR}/test1/orf_summary.txt" \
        "NM_012367" \
        "Test 1 ORF summary references the OR2B6 accession."

    # ---- Test 2: CCR7 with IUPAC ambiguity codes (accession NM_001838.4) --- #

    check_string_in_file "${TEST_OUTDIR}/test2/orfs.csv" \
        "NM_001838" \
        "Test 2 ORF CSV references the CCR7 accession (NM_001838.4)."

    check_csv_rows "${TEST_OUTDIR}/test2/orfs.csv" \
        1 \
        "Test 2 ORF CSV contains at least 1 ORF."

    check_string_in_file "${TEST_OUTDIR}/test2/cleaned_sequence_1.fasta" \
        "NM_001838" \
        "Test 2 cleaned FASTA header contains the CCR7 accession."

    check_string_in_file "${TEST_OUTDIR}/test2/orf_summary.txt" \
        "NM_001838" \
        "Test 2 ORF summary references the CCR7 accession."

    # ---- Test 3: comparative mode ---------------------------------------- #
    # Both accessions appear in the single combined orfs.csv and in their
    # respective cleaned FASTA files.

    check_string_in_file "${TEST_OUTDIR}/test3/orfs.csv" \
        "NM_012367" \
        "Test 3 combined ORF CSV contains the OR2B6 accession (NM_012367.1)."

    check_string_in_file "${TEST_OUTDIR}/test3/orfs.csv" \
        "NM_001838" \
        "Test 3 combined ORF CSV contains the CCR7 accession (NM_001838.4)."

    check_csv_rows "${TEST_OUTDIR}/test3/orfs.csv" \
        2 \
        "Test 3 combined ORF CSV contains ORFs from both sequences."

    check_string_in_file "${TEST_OUTDIR}/test3/comp_cleaned_sequence_1.fasta" \
        "NM_012367" \
        "Test 3 sequence-1 cleaned FASTA contains the OR2B6 accession."

    check_string_in_file "${TEST_OUTDIR}/test3/comp_cleaned_sequence_2.fasta" \
        "NM_001838" \
        "Test 3 sequence-2 cleaned FASTA contains the CCR7 accession."

    check_string_in_file "${TEST_OUTDIR}/test3/orf_comparison_report.txt" \
        "NM_012367" \
        "Test 3 comparison report references the OR2B6 accession."

    check_string_in_file "${TEST_OUTDIR}/test3/orf_comparison_report.txt" \
        "NM_001838" \
        "Test 3 comparison report references the CCR7 accession."
}


# ---- Summary -----------------------------------------------------------------

print_summary() {
    # Print a final tally of passed and failed checks.

    print_header "Test Summary"

    local total=$((PASSED + FAILED))
    echo "  Passed: ${PASSED} / ${total}"
    echo "  Failed: ${FAILED} / ${total}"
    echo ""

    if [[ $FAILED -eq 0 ]]; then
        echo "  All checks passed."
    else
        echo "  Some checks failed. Review the output above for details."
    fi
}


# ---- Main entry point --------------------------------------------------------

main() {
    # Orchestrate all test steps in order.

    preflight_checks
    run_tests
    print_header "Checking output files and output content"
    verify_output_files
    verify_output_content
    print_summary

    # Exit with a non-zero code if any check failed, so CI systems or the
    # user can detect the failure programmatically
    if [[ $FAILED -gt 0 ]]; then
        exit 1
    fi
}

# Run main only when the script is executed directly (not sourced)
main "$@"
