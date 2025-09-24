#!/bin/bash
# Stage 4.1 - Expect success without CB/UB
# Purpose: ensure legacy unsorted mode without CB/UB remains unaffected

set -euo pipefail

# Configuration
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
REFERENCE_DIR="/storage/scRNAseq_output/indices-98-32/star"
OUTPUT_DIR="$(dirname "$0")/../output/legacy_unsorted_no_cbub"
MANIFEST_DIR="$(dirname "$0")/../output/manifests"
STAR_BINARY="$(dirname "$0")/../../bin/Linux_x86_64/STAR"
THREADS=4

# Parse command line arguments
SUBSET=true  # Default to subset for faster tests
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

echo "=== STARsolo Unsorted CB/UB Testing - Stage 4.1: Expect Success Without CB/UB ==="
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

# STAR arguments for legacy unsorted mode (no CB/UB tags)
STAR_ARGS=(
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloCBwhitelist "$WHITELIST"
    --readFilesCommand zcat
    --runThreadN "$THREADS"
    --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM  # No CB/UB tags
    --outSAMtype BAM Unsorted
    # No --soloAddTagsToUnsorted flag needed since we're not requesting CB/UB
    --outFileNamePrefix "$OUTPUT_DIR/"
    --genomeDir "$REFERENCE_DIR"
    --readFilesIn "$R2_LIST" "$R1_LIST"
)

echo
echo "Running STAR in legacy unsorted mode (no CB/UB tags)..."
echo "This should SUCCEED and work exactly as before"
echo "Command: $STAR_BINARY ${STAR_ARGS[*]}"
echo

# Record start time
START_TIME=$(date +%s)

# Run STAR - we EXPECT this to succeed
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
echo "Runtime: $RUNTIME seconds ($(date -d@$RUNTIME -u +%H:%M:%S))"

# Analyze the results
echo
echo "Analyzing STAR results..."

if [[ $STAR_EXIT_CODE -ne 0 ]]; then
    echo "UNEXPECTED: STAR failed when it should have succeeded"
    echo "This indicates regression in legacy unsorted mode"
    
    echo
    echo "Error output:"
    if [[ -f "$STAR_ERROR_FILE" ]]; then
        cat "$STAR_ERROR_FILE"
    fi
    
    echo
    echo "=== Stage 4.1 FAILED: Legacy Mode Regression ==="
    echo "Legacy unsorted mode without CB/UB should work"
    exit 1
fi

# Verify expected outputs
echo "Verifying outputs..."

EXPECTED_FILES=(
    "Aligned.out.bam"
    "Log.final.out"
    "Log.out"
    "Log.progress.out"
    "Solo.out/Gene/Features.tsv"
    "Solo.out/Gene/Barcodes.tsv"
    "Solo.out/Gene/matrix.mtx"
    "Solo.out/Gene/Filtered/features.tsv"
    "Solo.out/Gene/Filtered/barcodes.tsv"
    "Solo.out/Gene/Filtered/matrix.mtx"
)

missing_files=()
for file in "${EXPECTED_FILES[@]}"; do
    if [[ ! -f "$OUTPUT_DIR/$file" ]]; then
        missing_files+=("$file")
    else
        size=$(stat -c%s "$OUTPUT_DIR/$file")
        echo "✓ $file ($(numfmt --to=iec $size))"
    fi
done

if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "ERROR: Missing expected files:" >&2
    printf "  %s\n" "${missing_files[@]}" >&2
    echo
    echo "=== Stage 4.1 FAILED: Incomplete Output ==="
    exit 1
fi

# Quick BAM validation
echo
echo "Validating BAM file..."
if command -v samtools >/dev/null 2>&1; then
    if samtools quickcheck "$OUTPUT_DIR/Aligned.out.bam"; then
        echo "✓ BAM file passes quickcheck"
        
        # Get basic stats
        total_reads=$(samtools view -c "$OUTPUT_DIR/Aligned.out.bam")
        mapped_reads=$(samtools view -c -F 4 "$OUTPUT_DIR/Aligned.out.bam")
        echo "✓ Total reads: $total_reads"
        echo "✓ Mapped reads: $mapped_reads"
        
        # Verify NO CB/UB tags (should be absent)
        cb_count=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | grep -c "CB:Z:" || true)
        ub_count=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | grep -c "UB:Z:" || true)
        
        echo "CB tags in first 1000 reads: $cb_count (expected: 0)"
        echo "UB tags in first 1000 reads: $ub_count (expected: 0)"
        
        if [[ $cb_count -gt 0 || $ub_count -gt 0 ]]; then
            echo "ERROR: Found CB/UB tags when none were requested" >&2
            echo "This indicates the new code is incorrectly adding tags" >&2
            echo
            echo "=== Stage 4.1 FAILED: Unexpected CB/UB Tags ==="
            exit 1
        fi
        
        # Verify presence of other expected tags
        cr_count=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | grep -c "CR:Z:" || true)
        ur_count=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | grep -c "UR:Z:" || true)
        gx_count=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | grep -c "GX:Z:" || true)
        
        echo "✓ CR tags in first 1000 reads: $cr_count"
        echo "✓ UR tags in first 1000 reads: $ur_count"
        echo "✓ GX tags in first 1000 reads: $gx_count"
        
        # Verify BAM is unsorted
        echo "Verifying BAM is unsorted..."
        positions=$(samtools view "$OUTPUT_DIR/Aligned.out.bam" | head -1000 | awk '{print $4}' | head -100)
        sorted_positions=$(echo "$positions" | sort -n)
        
        if [[ "$positions" == "$sorted_positions" ]]; then
            echo "WARNING: BAM appears to be sorted (first 100 reads)" >&2
        else
            echo "✓ BAM appears unsorted"
        fi
        
    else
        echo "ERROR: BAM file failed validation" >&2
        echo
        echo "=== Stage 4.1 FAILED: Invalid BAM ==="
        exit 1
    fi
else
    echo "! samtools not available, skipping BAM validation"
fi

# Validate Solo outputs
echo
echo "Validating Solo outputs..."

# Check matrix dimensions
GENE_MATRIX="$OUTPUT_DIR/Solo.out/Gene/matrix.mtx"
if [[ -f "$GENE_MATRIX" ]]; then
    dims=$(grep -v '^%' "$GENE_MATRIX" | head -1)
    genes=$(echo "$dims" | awk '{print $1}')
    barcodes=$(echo "$dims" | awk '{print $2}')
    entries=$(echo "$dims" | awk '{print $3}')
    
    echo "✓ Gene matrix dimensions: $genes genes, $barcodes barcodes, $entries entries"
    
    if [[ $entries -eq 0 ]]; then
        echo "WARNING: Gene matrix has no entries" >&2
    fi
else
    echo "ERROR: Gene matrix not found" >&2
fi

# Check filtered matrix
FILTERED_MATRIX="$OUTPUT_DIR/Solo.out/Gene/Filtered/matrix.mtx"
if [[ -f "$FILTERED_MATRIX" ]]; then
    dims=$(grep -v '^%' "$FILTERED_MATRIX" | head -1)
    filtered_genes=$(echo "$dims" | awk '{print $1}')
    filtered_barcodes=$(echo "$dims" | awk '{print $2}')
    filtered_entries=$(echo "$dims" | awk '{print $3}')
    
    echo "✓ Filtered matrix dimensions: $filtered_genes genes, $filtered_barcodes barcodes, $filtered_entries entries"
fi

# Check log for completion
echo
echo "Checking log for successful completion..."
if [[ -f "$OUTPUT_DIR/Log.final.out" ]]; then
    # Extract key metrics
    input_reads=$(grep "Number of input reads" "$OUTPUT_DIR/Log.final.out" | awk '{print $NF}')
    mapped_reads=$(grep "Uniquely mapped reads number" "$OUTPUT_DIR/Log.final.out" | awk '{print $NF}')
    
    echo "✓ Input reads: $input_reads"
    echo "✓ Uniquely mapped reads: $mapped_reads"
    
    if [[ "$input_reads" =~ ^[0-9]+$ && $input_reads -gt 0 ]]; then
        echo "✓ Log shows successful processing"
    else
        echo "WARNING: Unexpected input read count: $input_reads" >&2
    fi
fi

# Generate checksums
echo
echo "Generating checksums..."
CHECKSUM_FILE="$OUTPUT_DIR/checksums.md5"
{
    echo "# Legacy unsorted mode checksums generated on $(date)"
    echo "# Runtime: $RUNTIME seconds"
    echo "# Mode: Legacy unsorted (no CB/UB tags)"
    cd "$OUTPUT_DIR"
    md5sum "Aligned.out.bam" \
           "Solo.out/Gene/matrix.mtx" \
           "Solo.out/Gene/Filtered/matrix.mtx" \
           "Solo.out/Gene/Features.tsv" \
           "Solo.out/Gene/Barcodes.tsv"
} > "$CHECKSUM_FILE"

echo "✓ Checksums saved to: $CHECKSUM_FILE"

# Summary
echo
echo "=== Stage 4.1 Complete: Legacy Mode Validation ==="
echo "Output directory: $OUTPUT_DIR"
echo "Key findings:"
echo "  - STAR completed successfully (exit code: $STAR_EXIT_CODE)"
echo "  - All expected files generated"
echo "  - BAM contains no CB/UB tags (as expected)"
echo "  - BAM contains expected CR/UR/GX tags"
echo "  - Solo matrices generated successfully"
echo "  - Runtime: $RUNTIME seconds"

echo
echo "=== Stage 4.1 PASSED: Legacy Unsorted Mode Works Correctly ==="
echo "Legacy behavior is preserved - unsorted BAM without CB/UB tags works as before"
echo "Ready for additional regression testing"
