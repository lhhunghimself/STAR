#!/bin/bash
# Stage 3.2 - Compare Solo matrices
# Purpose: guarantee Solo outputs (Gene/GeneFull/filtered/unfiltered) match bit-for-bit

set -euo pipefail

# Configuration
BASELINE_DIR="$(dirname "$0")/../output/baseline_sorted"
UNSORTED_DIR="$(dirname "$0")/../output/unsorted_twopass"

# Parse command line arguments
VERBOSE=false
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
        --help)
            echo "Usage: $0 [--baseline-dir DIR] [--unsorted-dir DIR] [--verbose]"
            echo "  --baseline-dir   Directory containing sorted baseline outputs"
            echo "  --unsorted-dir   Directory containing unsorted outputs"
            echo "  --verbose        Enable verbose output"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

echo "=== STARsolo Unsorted CB/UB Testing - Stage 3.2: Compare Solo Matrices ==="
echo "Timestamp: $(date)"
echo "Baseline directory: $BASELINE_DIR"
echo "Unsorted directory: $UNSORTED_DIR"
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

# Define files to compare
MATRIX_FILES=(
    "Solo.out/Gene/Features.tsv"
    "Solo.out/Gene/Barcodes.tsv"
    "Solo.out/Gene/matrix.mtx"
    "Solo.out/Gene/Filtered/features.tsv"
    "Solo.out/Gene/Filtered/barcodes.tsv"
    "Solo.out/Gene/Filtered/matrix.mtx"
)

# Additional files that might exist
OPTIONAL_FILES=(
    "Solo.out/Gene/Summary.csv"
    "Solo.out/Gene/UMIperCellSorted.txt"
    "Solo.out/Gene/CellReads.stats"
    "Solo.out/GeneFull/Features.tsv"
    "Solo.out/GeneFull/Barcodes.tsv"
    "Solo.out/GeneFull/matrix.mtx"
    "Solo.out/GeneFull/Filtered/features.tsv"
    "Solo.out/GeneFull/Filtered/barcodes.tsv"
    "Solo.out/GeneFull/Filtered/matrix.mtx"
)

ALL_FILES=("${MATRIX_FILES[@]}" "${OPTIONAL_FILES[@]}")

# Function to compare two files
compare_files() {
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
    
    # Compare file sizes first
    baseline_size=$(stat -c%s "$baseline_file")
    unsorted_size=$(stat -c%s "$unsorted_file")
    
    if [[ $baseline_size -ne $unsorted_size ]]; then
        echo "  FAIL: $file (size mismatch: $baseline_size vs $unsorted_size bytes)"
        return 1
    fi
    
    # Compare checksums
    baseline_md5=$(md5sum "$baseline_file" | awk '{print $1}')
    unsorted_md5=$(md5sum "$unsorted_file" | awk '{print $1}')
    
    if [[ "$baseline_md5" != "$unsorted_md5" ]]; then
        echo "  FAIL: $file (checksum mismatch)"
        echo "    Baseline MD5: $baseline_md5"
        echo "    Unsorted MD5: $unsorted_md5"
        
        # Show first few differences if it's a text file
        if file "$baseline_file" | grep -q "text"; then
            echo "    First 10 differences:"
            diff "$baseline_file" "$unsorted_file" | head -10 | sed 's/^/      /'
        fi
        
        return 1
    fi
    
    if [[ "$VERBOSE" == "true" ]]; then
        echo "  PASS: $file ($(numfmt --to=iec $baseline_size), MD5: ${baseline_md5:0:8}...)"
    else
        echo "  PASS: $file"
    fi
    
    return 0
}

# Function to compare matrix dimensions and basic stats
analyze_matrix() {
    local file="$1"
    local baseline_file="$BASELINE_DIR/$file"
    local unsorted_file="$UNSORTED_DIR/$file"
    
    if [[ ! -f "$baseline_file" || ! -f "$unsorted_file" ]]; then
        return 0
    fi
    
    if [[ "$file" == *.mtx ]]; then
        echo "    Matrix analysis for $file:"
        
        # Extract dimensions from MTX files
        baseline_dims=$(grep -v '^%' "$baseline_file" | head -1)
        unsorted_dims=$(grep -v '^%' "$unsorted_file" | head -1)
        
        echo "      Baseline dimensions: $baseline_dims"
        echo "      Unsorted dimensions: $unsorted_dims"
        
        if [[ "$baseline_dims" != "$unsorted_dims" ]]; then
            echo "      WARNING: Dimension mismatch!"
        fi
        
        # Count non-zero entries
        baseline_entries=$(tail -n +4 "$baseline_file" | wc -l)
        unsorted_entries=$(tail -n +4 "$unsorted_file" | wc -l)
        
        echo "      Baseline entries: $baseline_entries"
        echo "      Unsorted entries: $unsorted_entries"
        
    elif [[ "$file" == *Barcodes.tsv || "$file" == *barcodes.tsv ]]; then
        baseline_count=$(wc -l < "$baseline_file")
        unsorted_count=$(wc -l < "$unsorted_file")
        
        echo "    Barcode count analysis for $file:"
        echo "      Baseline barcodes: $baseline_count"
        echo "      Unsorted barcodes: $unsorted_count"
        
    elif [[ "$file" == *Features.tsv || "$file" == *features.tsv ]]; then
        baseline_count=$(wc -l < "$baseline_file")
        unsorted_count=$(wc -l < "$unsorted_file")
        
        echo "    Feature count analysis for $file:"
        echo "      Baseline features: $baseline_count"
        echo "      Unsorted features: $unsorted_count"
    fi
}

echo "Comparing Solo matrix files..."
echo

failed_files=()
passed_files=()

for file in "${ALL_FILES[@]}"; do
    if compare_files "$file"; then
        passed_files+=("$file")
        
        # Detailed analysis for verbose mode
        if [[ "$VERBOSE" == "true" ]]; then
            analyze_matrix "$file"
        fi
    else
        failed_files+=("$file")
    fi
done

echo
echo "Comparison summary:"
echo "  Total files checked: $((${#passed_files[@]} + ${#failed_files[@]}))"
echo "  Passed: ${#passed_files[@]}"
echo "  Failed: ${#failed_files[@]}"

if [[ ${#failed_files[@]} -gt 0 ]]; then
    echo
    echo "Failed files:"
    for file in "${failed_files[@]}"; do
        echo "  - $file"
    done
fi

# Additional validation: check for unexpected files in either directory
echo
echo "Checking for unexpected Solo files..."

find_solo_files() {
    local dir="$1"
    find "$dir/Solo.out" -name "*.tsv" -o -name "*.mtx" -o -name "*.csv" -o -name "*.txt" -o -name "*.stats" 2>/dev/null | \
        sed "s|^$dir/||" | sort
}

if [[ -d "$BASELINE_DIR/Solo.out" && -d "$UNSORTED_DIR/Solo.out" ]]; then
    baseline_files=$(find_solo_files "$BASELINE_DIR")
    unsorted_files=$(find_solo_files "$UNSORTED_DIR")
    
    # Files only in baseline
    only_baseline=$(comm -23 <(echo "$baseline_files") <(echo "$unsorted_files"))
    if [[ -n "$only_baseline" ]]; then
        echo "Files only in baseline:"
        echo "$only_baseline" | sed 's/^/  /'
    fi
    
    # Files only in unsorted
    only_unsorted=$(comm -13 <(echo "$baseline_files") <(echo "$unsorted_files"))
    if [[ -n "$only_unsorted" ]]; then
        echo "Files only in unsorted:"
        echo "$only_unsorted" | sed 's/^/  /'
    fi
    
    if [[ -z "$only_baseline" && -z "$only_unsorted" ]]; then
        echo "✓ Both runs produced identical file sets"
    fi
fi

echo
if [[ ${#failed_files[@]} -eq 0 ]]; then
    echo "=== Stage 3.2 Complete: Solo Matrices Comparison PASSED ==="
    echo "All Solo output files match bit-for-bit between runs"
    exit 0
else
    echo "=== Stage 3.2 FAILED: Solo Matrices Comparison Failed ==="
    echo "${#failed_files[@]} files have differences between runs"
    exit 1
fi
