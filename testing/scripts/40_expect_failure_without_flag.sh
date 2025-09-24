#!/bin/bash
# Stage 4 - Expect failure without flag
# Purpose: confirm we still reject CB/UB on unsorted BAM when the new flag is absent

set -euo pipefail

# Configuration
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
REFERENCE_DIR="/storage/scRNAseq_output/indices-98-32/star"
OUTPUT_DIR="$(dirname "$0")/../output/negative_test_no_flag"
MANIFEST_DIR="$(dirname "$0")/../output/manifests"
STAR_BINARY="$(dirname "$0")/../../bin/Linux_x86_64/STAR"
THREADS=4

# Parse command line arguments
SUBSET=true  # Default to subset for faster negative tests
while [[ $# -gt 0 ]]; do
    case $1 in
        --full-data)
            SUBSET=false
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [--full-data] [--threads N]"
            echo "  --full-data  Use full dataset instead of subset"
            echo "  --threads    Number of threads to use (default: $THREADS)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

echo "=== STARsolo Unsorted CB/UB Testing - Stage 4: Expect Failure Without Flag ==="
echo "Timestamp: $(date)"
echo "Using subset: $SUBSET"
echo "Threads: $THREADS"
echo

# Check prerequisites
if [[ ! -x "$STAR_BINARY" ]]; then
    echo "ERROR: STAR binary not found or not executable: $STAR_BINARY" >&2
    exit 1
fi

# Select appropriate manifest
if [[ "$SUBSET" == "true" ]]; then
    MANIFEST_FILE="$MANIFEST_DIR/fastq_subset.list"
    echo "Using subset manifest: $MANIFEST_FILE"
else
    MANIFEST_FILE="$MANIFEST_DIR/fastq.list"
    echo "Using full manifest: $MANIFEST_FILE"
fi

if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "ERROR: Manifest file not found: $MANIFEST_FILE" >&2
    echo "Run 00_assert_fixtures.sh first to generate manifests" >&2
    exit 1
fi

# Create output directory
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Extract R1 and R2 file lists from manifest
R1_LIST=$(awk 'NR>2 {print $1}' "$MANIFEST_FILE" | paste -sd,)
R2_LIST=$(awk 'NR>2 {print $2}' "$MANIFEST_FILE" | paste -sd,)

if [[ -z "$R1_LIST" || -z "$R2_LIST" ]]; then
    echo "ERROR: Empty file lists from manifest" >&2
    exit 1
fi

echo "R1 files: $(echo "$R1_LIST" | tr ',' '\n' | wc -l) files"
echo "R2 files: $(echo "$R2_LIST" | tr ',' '\n' | wc -l) files"

# STAR arguments that SHOULD FAIL: unsorted BAM with CB/UB but no --soloAddTagsToUnsorted
STAR_ARGS=(
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloCBwhitelist "$WHITELIST"
    --readFilesCommand zcat
    --runThreadN "$THREADS"
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
    --outSAMtype BAM Unsorted
    # NOTE: intentionally omitting --soloAddTagsToUnsorted yes
    --outFileNamePrefix "$OUTPUT_DIR/"
    --genomeDir "$REFERENCE_DIR"
    --readFilesIn "$R2_LIST" "$R1_LIST"
)

echo
echo "Running STAR without --soloAddTagsToUnsorted flag..."
echo "This should FAIL with parameter validation error"
echo "Command: $STAR_BINARY ${STAR_ARGS[*]}"
echo

# Record start time
START_TIME=$(date +%s)

# Run STAR - we EXPECT this to fail
STAR_EXIT_CODE=0
STAR_OUTPUT_FILE="$OUTPUT_DIR/star_output.log"
STAR_ERROR_FILE="$OUTPUT_DIR/star_error.log"

if "$STAR_BINARY" "${STAR_ARGS[@]}" >"$STAR_OUTPUT_FILE" 2>"$STAR_ERROR_FILE"; then
    STAR_EXIT_CODE=0
else
    STAR_EXIT_CODE=$?
fi

# Record end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "STAR completed with exit code: $STAR_EXIT_CODE"
echo "Runtime: $RUNTIME seconds"

# Analyze the output
echo
echo "Analyzing STAR output..."

# Check if STAR failed as expected
if [[ $STAR_EXIT_CODE -eq 0 ]]; then
    echo "UNEXPECTED: STAR succeeded when it should have failed"
    echo "This indicates the parameter validation is not working correctly"
    
    # Show what was produced
    echo
    echo "Files produced:"
    find "$OUTPUT_DIR" -type f | head -10
    
    echo
    echo "=== Stage 4 FAILED: STAR Should Have Failed ===")
    echo "The parameter validation for CB/UB with unsorted BAM is not working"
    exit 1
fi

# Check error message content
echo "Checking error message content..."
if [[ -f "$STAR_ERROR_FILE" ]]; then
    echo "STDERR content:"
    cat "$STAR_ERROR_FILE" | head -20
    
    # Look for expected error patterns
    EXPECTED_PATTERNS=(
        "CB.*UB.*unsorted"
        "soloAddTagsToUnsorted"
        "parameter.*error"
        "invalid.*combination"
    )
    
    found_expected_error=false
    for pattern in "${EXPECTED_PATTERNS[@]}"; do
        if grep -i "$pattern" "$STAR_ERROR_FILE" >/dev/null; then
            echo "✓ Found expected error pattern: $pattern"
            found_expected_error=true
            break
        fi
    done
    
    if [[ "$found_expected_error" == "false" ]]; then
        echo "WARNING: Did not find expected error message pattern"
        echo "STAR failed but with unexpected error message"
    fi
else
    echo "No STDERR file found"
fi

if [[ -f "$STAR_OUTPUT_FILE" ]]; then
    echo
    echo "STDOUT content (first 20 lines):"
    head -20 "$STAR_OUTPUT_FILE"
fi

# Check that no significant output files were created
echo
echo "Checking output files..."
EXPECTED_MISSING_FILES=(
    "Aligned.out.bam"
    "Solo.out/Gene/matrix.mtx"
)

unexpected_files=()
for file in "${EXPECTED_MISSING_FILES[@]}"; do
    if [[ -f "$OUTPUT_DIR/$file" ]]; then
        size=$(stat -c%s "$OUTPUT_DIR/$file")
        unexpected_files+=("$file ($(numfmt --to=iec $size))")
    fi
done

if [[ ${#unexpected_files[@]} -gt 0 ]]; then
    echo "WARNING: Found unexpected output files:"
    for file in "${unexpected_files[@]}"; do
        echo "  $file"
    done
    echo "This suggests STAR ran further than expected before failing"
fi

# Test a variant without CB/UB to ensure basic functionality works
echo
echo "Testing variant without CB/UB tags (should succeed)..."

VARIANT_OUTPUT_DIR="$OUTPUT_DIR/variant_no_cb_ub"
mkdir -p "$VARIANT_OUTPUT_DIR"

VARIANT_STAR_ARGS=(
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloCBwhitelist "$WHITELIST"
    --readFilesCommand zcat
    --runThreadN "$THREADS"
    --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM  # No CB/UB
    --outSAMtype BAM Unsorted
    --outFileNamePrefix "$VARIANT_OUTPUT_DIR/"
    --genomeDir "$REFERENCE_DIR"
    --readFilesIn "$R2_LIST" "$R1_LIST"
)

VARIANT_START_TIME=$(date +%s)
VARIANT_EXIT_CODE=0

if "$STAR_BINARY" "${VARIANT_STAR_ARGS[@]}" >"$VARIANT_OUTPUT_DIR/star_output.log" 2>"$VARIANT_OUTPUT_DIR/star_error.log"; then
    VARIANT_EXIT_CODE=0
else
    VARIANT_EXIT_CODE=$?
fi

VARIANT_END_TIME=$(date +%s)
VARIANT_RUNTIME=$((VARIANT_END_TIME - VARIANT_START_TIME))

echo "Variant run completed with exit code: $VARIANT_EXIT_CODE"
echo "Variant runtime: $VARIANT_RUNTIME seconds"

if [[ $VARIANT_EXIT_CODE -ne 0 ]]; then
    echo "WARNING: Variant run without CB/UB also failed"
    echo "This might indicate a broader issue with the test setup"
    
    if [[ -f "$VARIANT_OUTPUT_DIR/star_error.log" ]]; then
        echo "Variant error:"
        head -10 "$VARIANT_OUTPUT_DIR/star_error.log"
    fi
else
    echo "✓ Variant run without CB/UB succeeded as expected"
    
    # Quick check of variant output
    if [[ -f "$VARIANT_OUTPUT_DIR/Aligned.out.bam" ]]; then
        size=$(stat -c%s "$VARIANT_OUTPUT_DIR/Aligned.out.bam")
        echo "  Produced BAM: $(numfmt --to=iec $size)"
    fi
fi

# Summary
echo
echo "=== Stage 4 Complete: Failure Test Analysis ==="
echo "Primary test (CB/UB without flag):"
echo "  Exit code: $STAR_EXIT_CODE (expected: non-zero)"
echo "  Runtime: $RUNTIME seconds"
echo "  Result: $([ $STAR_EXIT_CODE -ne 0 ] && echo "PASS - Failed as expected" || echo "FAIL - Should have failed")"

echo "Variant test (no CB/UB):"
echo "  Exit code: $VARIANT_EXIT_CODE (expected: 0)"
echo "  Runtime: $VARIANT_RUNTIME seconds"
echo "  Result: $([ $VARIANT_EXIT_CODE -eq 0 ] && echo "PASS - Succeeded as expected" || echo "FAIL - Should have succeeded")"

# Final verdict
if [[ $STAR_EXIT_CODE -ne 0 ]]; then
    if [[ $VARIANT_EXIT_CODE -eq 0 ]]; then
        echo
        echo "=== Stage 4 PASSED: Parameter Validation Working Correctly ==="
        echo "STAR correctly rejects CB/UB tags with unsorted BAM when --soloAddTagsToUnsorted is not specified"
        exit 0
    else
        echo
        echo "=== Stage 4 INCONCLUSIVE: Variant Test Also Failed ==="
        echo "Primary test failed as expected, but variant test suggests broader issues"
        exit 2
    fi
else
    echo
    echo "=== Stage 4 FAILED: Parameter Validation Not Working ==="
    echo "STAR should have rejected CB/UB tags with unsorted BAM"
    exit 1
fi
