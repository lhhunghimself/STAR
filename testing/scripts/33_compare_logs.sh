#!/bin/bash
# Stage 3.3 - Compare logs
# Purpose: diff STAR log artifacts

set -euo pipefail

# Configuration
BASELINE_DIR="$(dirname "$0")/../output/baseline_sorted"
UNSORTED_DIR="$(dirname "$0")/../output/unsorted_twopass"

# Parse command line arguments
VERBOSE=false
IGNORE_TIMESTAMPS=true
while [[ $# -gt 0 ]]; do
    case $1 in
        --baseline-dir)
            BASELINE_DIR="$2"
            shift 2
            ;;
        --unsorted-dir)
            UNSORTED_DIR="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --include-timestamps)
            IGNORE_TIMESTAMPS=false
            shift
            ;;
        --help)
            echo "Usage: $0 [--baseline-dir DIR] [--unsorted-dir DIR] [--verbose] [--include-timestamps]"
            echo "  --baseline-dir       Directory containing sorted baseline outputs"
            echo "  --unsorted-dir       Directory containing unsorted outputs"
            echo "  --verbose            Enable verbose output"
            echo "  --include-timestamps Include timestamps in comparison (default: ignore)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

echo "=== STARsolo Unsorted CB/UB Testing - Stage 3.3: Compare Logs ==="
echo "Timestamp: $(date)"
echo "Baseline directory: $BASELINE_DIR"
echo "Unsorted directory: $UNSORTED_DIR"
echo "Ignore timestamps: $IGNORE_TIMESTAMPS"
echo "Verbose: $VERBOSE"
echo

# Check input directories
if [[ ! -d "$BASELINE_DIR" ]]; then
    echo "ERROR: Baseline directory not found: $BASELINE_DIR" >&2
    exit 1
fi

if [[ ! -d "$UNSORTED_DIR" ]]; then
    echo "ERROR: Unsorted directory not found: $UNSORTED_DIR" >&2
    exit 1
fi

# Define log files to compare
LOG_FILES=(
    "Log.final.out"
    "Log.out"
    "Log.progress.out"
    "Solo.out/Gene/CellReads.stats"
    "Solo.out/Gene/Summary.csv"
    "Solo.out/Gene/UMIperCellSorted.txt"
)

# Optional log files
OPTIONAL_LOG_FILES=(
    "Solo.out/GeneFull/CellReads.stats"
    "Solo.out/GeneFull/Summary.csv"
    "Solo.out/GeneFull/UMIperCellSorted.txt"
    "Solo.out/Velocyto/CellReads.stats"
    "Solo.out/Velocyto/Summary.csv"
)

ALL_LOG_FILES=("${LOG_FILES[@]}" "${OPTIONAL_LOG_FILES[@]}")

# Function to clean log content for comparison
clean_log_content() {
    local file="$1"
    local content
    
    content=$(cat "$file")
    
    if [[ "$IGNORE_TIMESTAMPS" == "true" ]]; then
        # Remove timestamp lines and other time-dependent content
        content=$(echo "$content" | \
            grep -v "started at" | \
            grep -v "finished at" | \
            grep -v "Started job on" | \
            grep -v "Started mapping on" | \
            grep -v "Finished on" | \
            grep -v "Processing speed" | \
            grep -v "Speed of processing" | \
            grep -v "STAR version" | \
            grep -v "Command line" | \
            grep -v "Started STAR run" | \
            sed 's/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]/XX:XX:XX/g' | \
            sed 's/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/XXXX-XX-XX/g' | \
            sed 's/[0-9]\+\.[0-9]\+ seconds/X.X seconds/g' | \
            sed 's/[0-9]\+\.[0-9]\+ minutes/X.X minutes/g')
    fi
    
    echo "$content"
}

# Function to compare log files
compare_log_files() {
    local file="$1"
    local baseline_file="$BASELINE_DIR/$file"
    local unsorted_file="$UNSORTED_DIR/$file"
    
    if [[ ! -f "$baseline_file" && ! -f "$unsorted_file" ]]; then
        if [[ "$VERBOSE" == "true" ]]; then
            echo "  SKIP: $file (not present in either run)"
        fi
        return 0
    fi
    
    if [[ ! -f "$baseline_file" ]]; then
        echo "  FAIL: $file (missing in baseline)"
        return 1
    fi
    
    if [[ ! -f "$unsorted_file" ]]; then
        echo "  FAIL: $file (missing in unsorted)"
        return 1
    fi
    
    # Create temporary cleaned files
    local temp_dir=$(mktemp -d)
    local baseline_clean="$temp_dir/baseline_clean"
    local unsorted_clean="$temp_dir/unsorted_clean"
    
    clean_log_content "$baseline_file" > "$baseline_clean"
    clean_log_content "$unsorted_file" > "$unsorted_clean"
    
    # Compare cleaned content
    if diff -q "$baseline_clean" "$unsorted_clean" >/dev/null; then
        if [[ "$VERBOSE" == "true" ]]; then
            baseline_size=$(stat -c%s "$baseline_file")
            echo "  PASS: $file ($(numfmt --to=iec $baseline_size))"
        else
            echo "  PASS: $file"
        fi
        rm -rf "$temp_dir"
        return 0
    else
        echo "  FAIL: $file (content differs)"
        
        # Show differences
        echo "    Differences (ignoring timestamps):"
        diff "$baseline_clean" "$unsorted_clean" | head -20 | sed 's/^/      /'
        
        if [[ $(diff "$baseline_clean" "$unsorted_clean" | wc -l) -gt 20 ]]; then
            echo "      ... (showing first 20 lines of diff)"
        fi
        
        rm -rf "$temp_dir"
        return 1
    fi
}

# Function to extract and compare key metrics from Log.final.out
compare_final_log_metrics() {
    local baseline_log="$BASELINE_DIR/Log.final.out"
    local unsorted_log="$UNSORTED_DIR/Log.final.out"
    
    if [[ ! -f "$baseline_log" || ! -f "$unsorted_log" ]]; then
        echo "Cannot compare final log metrics (files missing)"
        return 1
    fi
    
    echo "Comparing key metrics from Log.final.out:"
    
    # Extract key metrics
    local metrics=(
        "Number of input reads"
        "Uniquely mapped reads number"
        "Uniquely mapped reads %"
        "Number of reads mapped to multiple loci"
        "% of reads mapped to multiple loci"
        "Number of reads unmapped: too many mismatches"
        "Number of reads unmapped: too short"
        "Number of reads unmapped: other"
    )
    
    local differences=0
    
    for metric in "${metrics[@]}"; do
        baseline_value=$(grep "$metric" "$baseline_log" | awk '{print $NF}' | tr -d '%')
        unsorted_value=$(grep "$metric" "$unsorted_log" | awk '{print $NF}' | tr -d '%')
        
        if [[ "$baseline_value" != "$unsorted_value" ]]; then
            echo "  DIFF: $metric"
            echo "    Baseline: $baseline_value"
            echo "    Unsorted: $unsorted_value"
            differences=$((differences + 1))
        elif [[ "$VERBOSE" == "true" ]]; then
            echo "  MATCH: $metric = $baseline_value"
        fi
    done
    
    if [[ $differences -eq 0 ]]; then
        echo "  ✓ All key metrics match"
        return 0
    else
        echo "  ✗ $differences metric(s) differ"
        return 1
    fi
}

echo "Comparing log files..."
echo

failed_files=()
passed_files=()

for file in "${ALL_LOG_FILES[@]}"; do
    if compare_log_files "$file"; then
        passed_files+=("$file")
    else
        failed_files+=("$file")
    fi
done

echo
echo "Detailed metrics comparison:"
if ! compare_final_log_metrics; then
    failed_files+=("Log.final.out metrics")
fi

echo
echo "Log comparison summary:"
echo "  Total files checked: $((${#passed_files[@]} + ${#failed_files[@]}))"
echo "  Passed: ${#passed_files[@]}"
echo "  Failed: ${#failed_files[@]}"

if [[ ${#failed_files[@]} -gt 0 ]]; then
    echo
    echo "Failed comparisons:"
    for file in "${failed_files[@]}"; do
        echo "  - $file"
    done
fi

# Additional analysis: check for unexpected log files
echo
echo "Checking for unexpected log files..."

find_log_files() {
    local dir="$1"
    find "$dir" -name "*.out" -o -name "*.log" -o -name "*.stats" -o -name "*.csv" -o -name "*.txt" 2>/dev/null | \
        grep -E "(Log\.|Solo\.out/)" | \
        sed "s|^$dir/||" | sort
}

baseline_logs=$(find_log_files "$BASELINE_DIR")
unsorted_logs=$(find_log_files "$UNSORTED_DIR")

# Files only in baseline
only_baseline=$(comm -23 <(echo "$baseline_logs") <(echo "$unsorted_logs"))
if [[ -n "$only_baseline" ]]; then
    echo "Log files only in baseline:"
    echo "$only_baseline" | sed 's/^/  /'
fi

# Files only in unsorted
only_unsorted=$(comm -13 <(echo "$baseline_logs") <(echo "$unsorted_logs"))
if [[ -n "$only_unsorted" ]]; then
    echo "Log files only in unsorted:"
    echo "$only_unsorted" | sed 's/^/  /'
fi

if [[ -z "$only_baseline" && -z "$only_unsorted" ]]; then
    echo "✓ Both runs produced identical log file sets"
fi

echo
if [[ ${#failed_files[@]} -eq 0 ]]; then
    echo "=== Stage 3.3 Complete: Log Comparison PASSED ==="
    echo "All log files and metrics match between runs"
    exit 0
else
    echo "=== Stage 3.3 FAILED: Log Comparison Failed ==="
    echo "${#failed_files[@]} log files or metrics have differences"
    exit 1
fi
