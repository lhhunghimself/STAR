#!/bin/bash
# Stage 4.2 - Subset regression loop
# Purpose: optional quick regression using the subset FASTQs from Stage 0

set -euo pipefail

# Configuration
OUTPUT_DIR="$(dirname "$0")/../output/subset_regression"
SCRIPTS_DIR="$(dirname "$0")"

echo "=== STARsolo Unsorted CB/UB Testing - Stage 4.2: Subset Regression Loop ==="
echo "Timestamp: $(date)"
echo "Output directory: $OUTPUT_DIR"
echo

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to run a script and capture results
run_script() {
    local script="$1"
    local description="$2"
    local expected_exit_code="${3:-0}"
    
    echo "--- Running $script: $description ---"
    
    local start_time=$(date +%s)
    local exit_code=0
    local log_file="$OUTPUT_DIR/$(basename "$script" .sh).log"
    
    if "$SCRIPTS_DIR/$script" --subset >"$log_file" 2>&1; then
        exit_code=0
    else
        exit_code=$?
    fi
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    
    echo "  Exit code: $exit_code (expected: $expected_exit_code)"
    echo "  Runtime: $runtime seconds"
    echo "  Log: $log_file"
    
    if [[ $exit_code -eq $expected_exit_code ]]; then
        echo "  Result: PASS"
        return 0
    else
        echo "  Result: FAIL"
        echo "  Last 10 lines of output:"
        tail -10 "$log_file" | sed 's/^/    /'
        return 1
    fi
}

# Function to run comparison scripts
run_comparison() {
    local script="$1"
    local description="$2"
    
    echo "--- Running $script: $description ---"
    
    local start_time=$(date +%s)
    local exit_code=0
    local log_file="$OUTPUT_DIR/$(basename "$script" .py).log"
    
    # Comparison scripts need to point to subset outputs
    if "$SCRIPTS_DIR/$script" \
        --baseline-dir "$SCRIPTS_DIR/../output/baseline_sorted_subset" \
        --unsorted-dir "$SCRIPTS_DIR/../output/unsorted_twopass_subset" \
        >"$log_file" 2>&1; then
        exit_code=0
    else
        exit_code=$?
    fi
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    
    echo "  Exit code: $exit_code"
    echo "  Runtime: $runtime seconds"
    echo "  Log: $log_file"
    
    if [[ $exit_code -eq 0 ]]; then
        echo "  Result: PASS"
        return 0
    else
        echo "  Result: FAIL"
        echo "  Last 10 lines of output:"
        tail -10 "$log_file" | sed 's/^/    /'
        return 1
    fi
}

# Check prerequisites
echo "Checking prerequisites..."

if [[ ! -f "$SCRIPTS_DIR/../output/manifests/fastq_subset.list" ]]; then
    echo "ERROR: Subset manifest not found. Run 01_smoke_subset_fastqs.sh first" >&2
    exit 1
fi

echo "✓ Subset manifest found"

# Track results
declare -a passed_tests=()
declare -a failed_tests=()
overall_start_time=$(date +%s)

echo
echo "Starting subset regression loop..."
echo "This runs the full pipeline on subset data for quick validation"
echo

# Stage 0: Fixtures (already done, just verify)
echo "=== Stage 0: Fixtures ==="
if [[ -f "$SCRIPTS_DIR/../output/manifests/fastq_subset.list" ]]; then
    echo "✓ Subset fixtures ready"
    passed_tests+=("Stage 0: Fixtures")
else
    echo "✗ Subset fixtures missing"
    failed_tests+=("Stage 0: Fixtures")
fi

# Stage 1: Sorted baseline with subset
echo
echo "=== Stage 1: Sorted Baseline (Subset) ==="
if run_script "10_run_sorted_baseline.sh" "Generate sorted baseline with subset data"; then
    passed_tests+=("Stage 1: Sorted baseline")
    
    # Move output to subset-specific directory
    if [[ -d "$SCRIPTS_DIR/../output/baseline_sorted" ]]; then
        mv "$SCRIPTS_DIR/../output/baseline_sorted" "$SCRIPTS_DIR/../output/baseline_sorted_subset"
        echo "  Moved output to baseline_sorted_subset"
    fi
    
    # Collect metrics
    if run_script "11_collect_sorted_metrics.sh" "Collect sorted metrics"; then
        passed_tests+=("Stage 1: Sorted metrics")
    else
        failed_tests+=("Stage 1: Sorted metrics")
    fi
else
    failed_tests+=("Stage 1: Sorted baseline")
fi

# Stage 2: Two-pass unsorted with subset
echo
echo "=== Stage 2: Two-Pass Unsorted (Subset) ==="
if run_script "20_run_unsorted_two_pass.sh" "Generate two-pass unsorted with subset data"; then
    passed_tests+=("Stage 2: Two-pass unsorted")
    
    # Move output to subset-specific directory
    if [[ -d "$SCRIPTS_DIR/../output/unsorted_twopass" ]]; then
        mv "$SCRIPTS_DIR/../output/unsorted_twopass" "$SCRIPTS_DIR/../output/unsorted_twopass_subset"
        echo "  Moved output to unsorted_twopass_subset"
    fi
    
    # Collect metrics
    if run_script "21_collect_unsorted_metrics.sh" "Collect unsorted metrics"; then
        passed_tests+=("Stage 2: Unsorted metrics")
    else
        failed_tests+=("Stage 2: Unsorted metrics")
    fi
else
    failed_tests+=("Stage 2: Two-pass unsorted")
fi

# Stage 3: Comparisons (only if both Stage 1 and 2 passed)
echo
echo "=== Stage 3: Comparisons (Subset) ==="
if [[ -d "$SCRIPTS_DIR/../output/baseline_sorted_subset" && -d "$SCRIPTS_DIR/../output/unsorted_twopass_subset" ]]; then
    
    # BAM tags comparison
    if run_comparison "30_compare_bam_tags.py" "Compare BAM tags"; then
        passed_tests+=("Stage 3: BAM tags")
    else
        failed_tests+=("Stage 3: BAM tags")
    fi
    
    # Alignment QC comparison
    if run_comparison "31_compare_alignment_qc.py" "Compare alignment QC"; then
        passed_tests+=("Stage 3: Alignment QC")
    else
        failed_tests+=("Stage 3: Alignment QC")
    fi
    
    # Solo matrices comparison
    if "$SCRIPTS_DIR/32_compare_solo_matrices.sh" \
        --baseline-dir "$SCRIPTS_DIR/../output/baseline_sorted_subset" \
        --unsorted-dir "$SCRIPTS_DIR/../output/unsorted_twopass_subset" \
        >"$OUTPUT_DIR/32_compare_solo_matrices.log" 2>&1; then
        passed_tests+=("Stage 3: Solo matrices")
        echo "  Solo matrices comparison: PASS"
    else
        failed_tests+=("Stage 3: Solo matrices")
        echo "  Solo matrices comparison: FAIL"
    fi
    
    # Logs comparison
    if "$SCRIPTS_DIR/33_compare_logs.sh" \
        --baseline-dir "$SCRIPTS_DIR/../output/baseline_sorted_subset" \
        --unsorted-dir "$SCRIPTS_DIR/../output/unsorted_twopass_subset" \
        >"$OUTPUT_DIR/33_compare_logs.log" 2>&1; then
        passed_tests+=("Stage 3: Logs")
        echo "  Logs comparison: PASS"
    else
        failed_tests+=("Stage 3: Logs")
        echo "  Logs comparison: FAIL"
    fi
    
    # Read distributions comparison
    if run_comparison "34_compare_read_distributions.py" "Compare read distributions"; then
        passed_tests+=("Stage 3: Read distributions")
    else
        failed_tests+=("Stage 3: Read distributions")
    fi
    
else
    echo "Skipping comparisons - missing baseline or unsorted outputs"
    failed_tests+=("Stage 3: Comparisons (skipped)")
fi

# Stage 4: Negative tests
echo
echo "=== Stage 4: Negative Tests (Subset) ==="

# Test failure without flag (should fail)
if run_script "40_expect_failure_without_flag.sh" "Expect failure without flag" 1; then
    passed_tests+=("Stage 4: Expect failure")
else
    failed_tests+=("Stage 4: Expect failure")
fi

# Test success without CB/UB (should succeed)
if run_script "41_expect_success_without_cbub.sh" "Expect success without CB/UB"; then
    passed_tests+=("Stage 4: Legacy mode")
else
    failed_tests+=("Stage 4: Legacy mode")
fi

# Calculate total runtime
overall_end_time=$(date +%s)
total_runtime=$((overall_end_time - overall_start_time))

echo
echo "=== Subset Regression Loop Complete ==="
echo "Total runtime: $total_runtime seconds ($(date -d@$total_runtime -u +%H:%M:%S))"
echo
echo "Results summary:"
echo "  Passed tests: ${#passed_tests[@]}"
echo "  Failed tests: ${#failed_tests[@]}"
echo "  Total tests: $((${#passed_tests[@]} + ${#failed_tests[@]}))"

if [[ ${#passed_tests[@]} -gt 0 ]]; then
    echo
    echo "Passed tests:"
    for test in "${passed_tests[@]}"; do
        echo "  ✓ $test"
    done
fi

if [[ ${#failed_tests[@]} -gt 0 ]]; then
    echo
    echo "Failed tests:"
    for test in "${failed_tests[@]}"; do
        echo "  ✗ $test"
    done
fi

# Generate summary report
REPORT_FILE="$OUTPUT_DIR/regression_report.json"
{
    echo "{"
    echo "  \"timestamp\": \"$(date -Iseconds)\","
    echo "  \"total_runtime_seconds\": $total_runtime,"
    echo "  \"total_tests\": $((${#passed_tests[@]} + ${#failed_tests[@]})),"
    echo "  \"passed_tests\": ${#passed_tests[@]},"
    echo "  \"failed_tests\": ${#failed_tests[@]},"
    echo "  \"success_rate\": $(echo "scale=4; ${#passed_tests[@]} * 100 / (${#passed_tests[@]} + ${#failed_tests[@]})" | bc),"
    echo "  \"passed\": ["
    for i in "${!passed_tests[@]}"; do
        echo -n "    \"${passed_tests[$i]}\""
        [[ $i -lt $((${#passed_tests[@]} - 1)) ]] && echo "," || echo
    done
    echo "  ],"
    echo "  \"failed\": ["
    for i in "${!failed_tests[@]}"; do
        echo -n "    \"${failed_tests[$i]}\""
        [[ $i -lt $((${#failed_tests[@]} - 1)) ]] && echo "," || echo
    done
    echo "  ]"
    echo "}"
} > "$REPORT_FILE"

echo
echo "Detailed report saved to: $REPORT_FILE"
echo "Individual test logs in: $OUTPUT_DIR/"

# Final result
if [[ ${#failed_tests[@]} -eq 0 ]]; then
    echo
    echo "=== Stage 4.2 PASSED: All Subset Regression Tests Passed ==="
    echo "The two-pass unsorted mode works correctly on subset data"
    exit 0
else
    echo
    echo "=== Stage 4.2 FAILED: Some Regression Tests Failed ==="
    echo "${#failed_tests[@]} out of $((${#passed_tests[@]} + ${#failed_tests[@]})) tests failed"
    exit 1
fi
