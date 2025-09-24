#!/bin/bash
# Stage 5 - Run all tests
# Purpose: top-level driver to execute the entire matrix of tests

set -euo pipefail

# Configuration
SCRIPTS_DIR="$(dirname "$0")"
OUTPUT_DIR="$SCRIPTS_DIR/../output"
REPORT_FILE="$OUTPUT_DIR/test-report.json"

# Default settings
USE_SUBSET=false
PARALLEL_COMPARISONS=true
KEEP_TEMP_FILES=false
VERBOSE=false
STOP_ON_FAILURE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --subset)
            USE_SUBSET=true
            shift
            ;;
        --full)
            USE_SUBSET=false
            shift
            ;;
        --serial)
            PARALLEL_COMPARISONS=false
            shift
            ;;
        --keep-temp)
            KEEP_TEMP_FILES=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --stop-on-failure)
            STOP_ON_FAILURE=true
            shift
            ;;
        --help)
            cat << EOF
Usage: $0 [OPTIONS]

STARsolo Two-Pass Unsorted CB/UB Testing - Complete Test Suite

OPTIONS:
  --subset              Use subset FASTQs for faster testing (default: full data)
  --full                Use full dataset (default)
  --serial              Run comparisons serially instead of parallel
  --keep-temp           Keep temporary files for debugging
  --verbose             Enable verbose output
  --stop-on-failure     Stop execution on first test failure

DESCRIPTION:
  This script orchestrates the complete testing pipeline for the STARsolo
  two-pass unsorted CB/UB feature. It runs through all stages:
  
  Stage 0: Environment and fixture validation
  Stage 1: Sorted baseline generation and metrics collection
  Stage 2: Two-pass unsorted run and metrics collection
  Stage 3: Cross-run comparisons (BAM tags, QC, matrices, logs, distributions)
  Stage 4: Negative tests and regression checks
  
  The script generates a comprehensive test report with pass/fail status
  for each component and overall timing information.

EXAMPLES:
  $0                    # Run full test suite with complete dataset
  $0 --subset           # Run with subset data for faster testing
  $0 --subset --verbose # Run subset with detailed output
  $0 --stop-on-failure  # Stop at first failure for debugging

EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Use --help for usage information" >&2
            exit 1
            ;;
    esac
done

echo "======================================================================"
echo "STARsolo Two-Pass Unsorted CB/UB Testing - Complete Test Suite"
echo "======================================================================"
echo "Timestamp: $(date)"
echo "Hostname: $(hostname)"
echo "User: $(whoami)"
echo "Working directory: $(pwd)"
echo "Scripts directory: $SCRIPTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo
echo "Configuration:"
echo "  Dataset: $([ "$USE_SUBSET" == "true" ] && echo "Subset (fast)" || echo "Full (complete)")"
echo "  Comparisons: $([ "$PARALLEL_COMPARISONS" == "true" ] && echo "Parallel" || echo "Serial")"
echo "  Keep temp files: $KEEP_TEMP_FILES"
echo "  Verbose output: $VERBOSE"
echo "  Stop on failure: $STOP_ON_FAILURE"
echo "======================================================================"

# Initialize test tracking
declare -a test_results=()
declare -A test_timings=()
declare -A test_logs=()
overall_start_time=$(date +%s)
failed_tests=0
total_tests=0

# Function to run a test and track results
run_test() {
    local test_id="$1"
    local test_name="$2"
    local command="$3"
    local expected_exit_code="${4:-0}"
    
    echo
    echo "--- Test $test_id: $test_name ---"
    
    local start_time=$(date +%s)
    local exit_code=0
    local log_file="$OUTPUT_DIR/logs/$(echo "$test_id" | tr ' ' '_' | tr '[:upper:]' '[:lower:]').log"
    
    # Create logs directory
    mkdir -p "$(dirname "$log_file")"
    
    # Run the command
    if eval "$command" >"$log_file" 2>&1; then
        exit_code=0
    else
        exit_code=$?
    fi
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    
    # Store results
    test_timings["$test_id"]=$runtime
    test_logs["$test_id"]="$log_file"
    total_tests=$((total_tests + 1))
    
    echo "  Runtime: $runtime seconds ($(date -d@$runtime -u +%H:%M:%S))"
    echo "  Exit code: $exit_code (expected: $expected_exit_code)"
    echo "  Log: $log_file"
    
    if [[ $exit_code -eq $expected_exit_code ]]; then
        echo "  Result: PASS ✓"
        test_results+=("$test_id:PASS:$runtime:$exit_code")
        return 0
    else
        echo "  Result: FAIL ✗"
        test_results+=("$test_id:FAIL:$runtime:$exit_code")
        failed_tests=$((failed_tests + 1))
        
        # Show last few lines of output on failure
        echo "  Last 10 lines of output:"
        tail -10 "$log_file" 2>/dev/null | sed 's/^/    /' || echo "    (no output available)"
        
        if [[ "$STOP_ON_FAILURE" == "true" ]]; then
            echo
            echo "STOPPING: --stop-on-failure specified and test failed"
            echo "Check log file: $log_file"
            exit 1
        fi
        
        return 1
    fi
}

# Function to build command with common options
build_command() {
    local script="$1"
    local extra_args="${2:-}"
    
    local cmd="$SCRIPTS_DIR/$script"
    
    if [[ "$USE_SUBSET" == "true" ]]; then
        cmd="$cmd --subset"
    fi
    
    if [[ "$VERBOSE" == "true" && "$script" =~ \.(sh|py)$ ]]; then
        if [[ "$script" =~ \.py$ ]]; then
            cmd="$cmd --verbose"
        elif [[ "$script" =~ \.sh$ ]] && grep -q "verbose" "$SCRIPTS_DIR/$script"; then
            cmd="$cmd --verbose"
        fi
    fi
    
    if [[ "$KEEP_TEMP_FILES" == "true" && "$script" =~ \.(py)$ ]]; then
        cmd="$cmd --keep-temp"
    fi
    
    if [[ -n "$extra_args" ]]; then
        cmd="$cmd $extra_args"
    fi
    
    echo "$cmd"
}

echo
echo "Starting test execution..."

# Stage 0: Environment and fixtures
echo
echo "======================================================================"
echo "STAGE 0: Environment and Fixture Validation"
echo "======================================================================"

run_test "0.1" "Assert fixtures" "$(build_command "00_assert_fixtures.sh")"

if [[ "$USE_SUBSET" == "true" ]]; then
    run_test "0.2" "Create subset FASTQs" "$(build_command "01_smoke_subset_fastqs.sh")"
fi

# Stage 1: Sorted baseline
echo
echo "======================================================================"
echo "STAGE 1: Sorted Baseline Generation"
echo "======================================================================"

run_test "1.1" "Run sorted baseline" "$(build_command "10_run_sorted_baseline.sh")"
run_test "1.2" "Collect sorted metrics" "$(build_command "11_collect_sorted_metrics.sh")"

# Stage 2: Two-pass unsorted
echo
echo "======================================================================"
echo "STAGE 2: Two-Pass Unsorted Generation"  
echo "======================================================================"

# Add keep-tmp flag if requested for debugging
extra_args=""
if [[ "$KEEP_TEMP_FILES" == "true" ]]; then
    extra_args="--keep-tmp"
fi

run_test "2.1" "Run two-pass unsorted" "$(build_command "20_run_unsorted_two_pass.sh" "$extra_args")"
run_test "2.2" "Collect unsorted metrics" "$(build_command "21_collect_unsorted_metrics.sh")"

# Optional: Verify temp trailers if temp files were kept
if [[ "$KEEP_TEMP_FILES" == "true" ]]; then
    run_test "2.3" "Verify temp trailers" "$(build_command "22_verify_tmp_trailers.py")"
fi

# Stage 3: Cross-run comparisons
echo
echo "======================================================================"
echo "STAGE 3: Cross-Run Comparisons"
echo "======================================================================"

if [[ "$PARALLEL_COMPARISONS" == "true" ]]; then
    echo "Running comparisons in parallel..."
    
    # Start all comparison tests in background
    declare -a comparison_pids=()
    
    # BAM tags comparison
    (run_test "3.1" "Compare BAM tags" "$(build_command "30_compare_bam_tags.py")") &
    comparison_pids+=($!)
    
    # Alignment QC comparison  
    (run_test "3.2" "Compare alignment QC" "$(build_command "31_compare_alignment_qc.py")") &
    comparison_pids+=($!)
    
    # Solo matrices comparison
    (run_test "3.3" "Compare Solo matrices" "$(build_command "32_compare_solo_matrices.sh")") &
    comparison_pids+=($!)
    
    # Logs comparison
    (run_test "3.4" "Compare logs" "$(build_command "33_compare_logs.sh")") &
    comparison_pids+=($!)
    
    # Read distributions comparison
    (run_test "3.5" "Compare read distributions" "$(build_command "34_compare_read_distributions.py")") &
    comparison_pids+=($!)
    
    # Wait for all comparisons to complete
    echo "Waiting for parallel comparisons to complete..."
    for pid in "${comparison_pids[@]}"; do
        wait $pid
    done
    echo "All comparisons completed"
    
else
    # Run comparisons serially
    echo "Running comparisons serially..."
    
    run_test "3.1" "Compare BAM tags" "$(build_command "30_compare_bam_tags.py")"
    run_test "3.2" "Compare alignment QC" "$(build_command "31_compare_alignment_qc.py")"
    run_test "3.3" "Compare Solo matrices" "$(build_command "32_compare_solo_matrices.sh")"
    run_test "3.4" "Compare logs" "$(build_command "33_compare_logs.sh")"
    run_test "3.5" "Compare read distributions" "$(build_command "34_compare_read_distributions.py")"
fi

# Stage 4: Negative tests
echo
echo "======================================================================"
echo "STAGE 4: Negative Tests and Regression Checks"
echo "======================================================================"

# These tests have specific expected exit codes
run_test "4.1" "Expect failure without flag" "$(build_command "40_expect_failure_without_flag.sh")" 1
run_test "4.2" "Expect success without CB/UB" "$(build_command "41_expect_success_without_cbub.sh")"

# Optional: Full subset regression loop
if [[ "$USE_SUBSET" == "true" ]]; then
    run_test "4.3" "Subset regression loop" "$(build_command "42_subset_regression_loop.sh")"
fi

# Calculate final statistics
overall_end_time=$(date +%s)
total_runtime=$((overall_end_time - overall_start_time))
passed_tests=$((total_tests - failed_tests))

echo
echo "======================================================================"
echo "TEST EXECUTION COMPLETE"
echo "======================================================================"
echo "Total runtime: $total_runtime seconds ($(date -d@$total_runtime -u +%H:%M:%S))"
echo "Total tests: $total_tests"
echo "Passed: $passed_tests"
echo "Failed: $failed_tests"
echo "Success rate: $(echo "scale=1; $passed_tests * 100 / $total_tests" | bc)%"

# Generate detailed test report
echo
echo "Generating test report..."

{
    echo "{"
    echo "  \"metadata\": {"
    echo "    \"timestamp\": \"$(date -Iseconds)\","
    echo "    \"hostname\": \"$(hostname)\","
    echo "    \"user\": \"$(whoami)\","
    echo "    \"working_directory\": \"$(pwd)\","
    echo "    \"configuration\": {"
    echo "      \"use_subset\": $USE_SUBSET,"
    echo "      \"parallel_comparisons\": $PARALLEL_COMPARISONS,"
    echo "      \"keep_temp_files\": $KEEP_TEMP_FILES,"
    echo "      \"verbose\": $VERBOSE,"
    echo "      \"stop_on_failure\": $STOP_ON_FAILURE"
    echo "    }"
    echo "  },"
    echo "  \"summary\": {"
    echo "    \"total_runtime_seconds\": $total_runtime,"
    echo "    \"total_tests\": $total_tests,"
    echo "    \"passed_tests\": $passed_tests,"
    echo "    \"failed_tests\": $failed_tests,"
    echo "    \"success_rate\": $(echo "scale=4; $passed_tests * 100 / $total_tests" | bc)"
    echo "  },"
    echo "  \"tests\": ["
    
    for i in "${!test_results[@]}"; do
        result="${test_results[$i]}"
        IFS=':' read -r test_id status runtime exit_code <<< "$result"
        
        echo "    {"
        echo "      \"id\": \"$test_id\","
        echo "      \"status\": \"$status\","
        echo "      \"runtime_seconds\": $runtime,"
        echo "      \"exit_code\": $exit_code,"
        echo "      \"log_file\": \"${test_logs[$test_id]}\""
        echo -n "    }"
        
        # Add comma if not last item
        [[ $i -lt $((${#test_results[@]} - 1)) ]] && echo "," || echo
    done
    
    echo "  ]"
    echo "}"
} > "$REPORT_FILE"

echo "✓ Test report saved to: $REPORT_FILE"

# Show summary by stage
echo
echo "Results by stage:"

declare -A stage_stats
for result in "${test_results[@]}"; do
    IFS=':' read -r test_id status runtime exit_code <<< "$result"
    stage=$(echo "$test_id" | cut -d'.' -f1)
    
    if [[ -z "${stage_stats[$stage]+x}" ]]; then
        stage_stats[$stage]="0:0"  # passed:failed
    fi
    
    IFS=':' read -r passed failed <<< "${stage_stats[$stage]}"
    if [[ "$status" == "PASS" ]]; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
    fi
    stage_stats[$stage]="$passed:$failed"
done

for stage in $(printf '%s\n' "${!stage_stats[@]}" | sort -n); do
    IFS=':' read -r passed failed <<< "${stage_stats[$stage]}"
    total=$((passed + failed))
    
    stage_name="Stage $stage"
    case $stage in
        0) stage_name="Stage 0 (Fixtures)" ;;
        1) stage_name="Stage 1 (Sorted Baseline)" ;;
        2) stage_name="Stage 2 (Two-Pass Unsorted)" ;;
        3) stage_name="Stage 3 (Comparisons)" ;;
        4) stage_name="Stage 4 (Negative Tests)" ;;
    esac
    
    echo "  $stage_name: $passed/$total passed $([ $failed -eq 0 ] && echo "✓" || echo "✗")"
done

# Show failed tests in detail
if [[ $failed_tests -gt 0 ]]; then
    echo
    echo "Failed tests:"
    for result in "${test_results[@]}"; do
        IFS=':' read -r test_id status runtime exit_code <<< "$result"
        if [[ "$status" == "FAIL" ]]; then
            echo "  ✗ $test_id (exit code: $exit_code, runtime: ${runtime}s)"
            echo "    Log: ${test_logs[$test_id]}"
        fi
    done
fi

# Cleanup temp files if not requested to keep them
if [[ "$KEEP_TEMP_FILES" == "false" ]]; then
    echo
    echo "Cleaning up temporary files..."
    # Remove any temp directories created during testing
    find "$OUTPUT_DIR" -name "*_temp_*" -type d -exec rm -rf {} + 2>/dev/null || true
    echo "✓ Temporary files cleaned up"
fi

echo
echo "======================================================================"
if [[ $failed_tests -eq 0 ]]; then
    echo "🎉 ALL TESTS PASSED! 🎉"
    echo "The STARsolo two-pass unsorted CB/UB feature is working correctly"
    echo "======================================================================"
    exit 0
else
    echo "❌ SOME TESTS FAILED ❌"
    echo "$failed_tests out of $total_tests tests failed"
    echo "Check the test report and individual log files for details"
    echo "======================================================================"
    exit 1
fi
