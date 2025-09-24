#!/bin/bash
# Stage 1 - Run sorted baseline
# Purpose: reproduce today's official STARsolo output for comparison

set -euo pipefail

# Configuration
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
REFERENCE_DIR="/storage/scRNAseq_output/indices-98-32/star"
OUTPUT_DIR="$(dirname "$0")/../output/baseline_sorted"
MANIFEST_DIR="$(dirname "$0")/../output/manifests"
STAR_BINARY="$(dirname "$0")/../../bin/Linux_x86_64/STAR"
THREADS=16

# Parse command line arguments
SUBSET=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --subset)
            SUBSET=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [--subset] [--threads N]"
            echo "  --subset    Use subset FASTQs for faster testing"
            echo "  --threads   Number of threads to use (default: $THREADS)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

echo "=== STARsolo Unsorted CB/UB Testing - Stage 1: Sorted Baseline ==="
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

# Default STAR arguments for STARsolo
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
    --outSAMtype BAM SortedByCoordinate
    --outFileNamePrefix "$OUTPUT_DIR/"
    --genomeDir "$REFERENCE_DIR"
    --readFilesIn "$R2_LIST" "$R1_LIST"
)

echo
echo "Running STAR with sorted BAM output..."
echo "Command: $STAR_BINARY ${STAR_ARGS[*]}"
echo "Output directory: $OUTPUT_DIR"
echo

# Record start time
START_TIME=$(date +%s)

# Run STAR
if ! "$STAR_BINARY" "${STAR_ARGS[@]}"; then
    echo "ERROR: STAR failed" >&2
    exit 1
fi

# Record end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo
echo "STAR completed successfully"
echo "Runtime: $RUNTIME seconds ($(date -d@$RUNTIME -u +%H:%M:%S))"

# Verify expected outputs
echo
echo "Verifying outputs..."

EXPECTED_FILES=(
    "Aligned.sortedByCoord.out.bam"
    "Log.final.out"
    "Log.out"
    "Log.progress.out"
    "Solo.out/Gene/raw/features.tsv"
    "Solo.out/Gene/raw/barcodes.tsv"
    "Solo.out/Gene/raw/matrix.mtx"
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
    exit 1
fi

# Generate checksums for key files
echo
echo "Generating checksums..."
CHECKSUM_FILE="$OUTPUT_DIR/checksums.md5"
{
    echo "# Checksums generated on $(date)"
    echo "# Runtime: $RUNTIME seconds"
    cd "$OUTPUT_DIR"
    md5sum "Aligned.sortedByCoord.out.bam" \
           "Solo.out/Gene/matrix.mtx" \
           "Solo.out/Gene/Filtered/matrix.mtx" \
           "Solo.out/Gene/Features.tsv" \
           "Solo.out/Gene/Barcodes.tsv"
} > "$CHECKSUM_FILE"

echo "✓ Checksums saved to: $CHECKSUM_FILE"

# Quick BAM validation
echo
echo "Validating BAM file..."
if command -v samtools >/dev/null 2>&1; then
    if samtools quickcheck "$OUTPUT_DIR/Aligned.sortedByCoord.out.bam"; then
        echo "✓ BAM file passes quickcheck"
        
        # Get basic stats
        total_reads=$(samtools view -c "$OUTPUT_DIR/Aligned.sortedByCoord.out.bam")
        mapped_reads=$(samtools view -c -F 4 "$OUTPUT_DIR/Aligned.sortedByCoord.out.bam")
        echo "✓ Total reads: $total_reads"
        echo "✓ Mapped reads: $mapped_reads"
        
        # Check for CB/UB tags in first few reads
        cb_count=$(samtools view "$OUTPUT_DIR/Aligned.sortedByCoord.out.bam" | head -1000 | grep -c "CB:Z:" || true)
        ub_count=$(samtools view "$OUTPUT_DIR/Aligned.sortedByCoord.out.bam" | head -1000 | grep -c "UB:Z:" || true)
        echo "✓ CB tags in first 1000 reads: $cb_count"
        echo "✓ UB tags in first 1000 reads: $ub_count"
        
        if [[ $cb_count -eq 0 || $ub_count -eq 0 ]]; then
            echo "WARNING: Few or no CB/UB tags found in sample reads" >&2
        fi
    else
        echo "ERROR: BAM file failed validation" >&2
        exit 1
    fi
else
    echo "! samtools not available, skipping BAM validation"
fi

# Summary
echo
echo "=== Stage 1 Complete: Sorted Baseline Generated ==="
echo "Output directory: $OUTPUT_DIR"
echo "Key files:"
echo "  - BAM: Aligned.sortedByCoord.out.bam"
echo "  - Solo matrices: Solo.out/Gene/"
echo "  - Checksums: checksums.md5"
echo "Runtime: $RUNTIME seconds"
echo "Ready for metrics collection (Stage 1.1)"
