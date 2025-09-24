#!/bin/bash
# Stage 0 - Create subset FASTQs
# Purpose: create a tiny subset (first 1M reads per mate) for quicker local dry runs

set -euo pipefail

# Configuration
FASTQ_DIR="/storage/downsampled/SC2300771"
SUBSET_DIR="$FASTQ_DIR/subset"
SUBSET_READS=1000000  # 1M reads per file
OUTPUT_DIR="$(dirname "$0")/../output"
MANIFEST_DIR="$OUTPUT_DIR/manifests"

echo "=== STARsolo Unsorted CB/UB Testing - Stage 0: Create Subset FASTQs ==="
echo "Timestamp: $(date)"
echo "Subset size: $SUBSET_READS reads per file"
echo

# Check if seqtk is available
SEQTK_AVAILABLE=false
if command -v seqtk >/dev/null 2>&1; then
    SEQTK_AVAILABLE=true
    echo "✓ seqtk found, will use for efficient subsampling"
else
    echo "! seqtk not found, will use zcat/head fallback"
fi

# Create subset directory
mkdir -p "$SUBSET_DIR"

# Check if subset already exists and is recent
if [[ -f "$SUBSET_DIR/R1.fastq.gz" && -f "$SUBSET_DIR/R2.fastq.gz" ]]; then
    r1_age=$(stat -c %Y "$SUBSET_DIR/R1.fastq.gz" 2>/dev/null || echo 0)
    r2_age=$(stat -c %Y "$SUBSET_DIR/R2.fastq.gz" 2>/dev/null || echo 0)
    current_time=$(date +%s)
    
    # If files are less than 24 hours old, skip regeneration
    if [[ $((current_time - r1_age)) -lt 86400 && $((current_time - r2_age)) -lt 86400 ]]; then
        echo "✓ Recent subset files found, skipping regeneration"
        echo "  - $SUBSET_DIR/R1.fastq.gz ($(stat -c %s "$SUBSET_DIR/R1.fastq.gz" | numfmt --to=iec) bytes)"
        echo "  - $SUBSET_DIR/R2.fastq.gz ($(stat -c %s "$SUBSET_DIR/R2.fastq.gz" | numfmt --to=iec) bytes)"
        exit 0
    fi
fi

# Load manifest if available, otherwise discover files
if [[ -f "$MANIFEST_DIR/fastq.list" ]]; then
    echo "Using existing manifest: $MANIFEST_DIR/fastq.list"
    R1_FILES=($(awk 'NR>2 {print $1}' "$MANIFEST_DIR/fastq.list"))
    R2_FILES=($(awk 'NR>2 {print $2}' "$MANIFEST_DIR/fastq.list"))
else
    echo "No manifest found, discovering FASTQ files..."
    R1_FILES=($(find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" -type f | sort))
    R2_FILES=($(find "$FASTQ_DIR" -name "*_R2_001.fastq.gz" -type f | sort))
fi

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 files found" >&2
    exit 1
fi

echo "Found ${#R1_FILES[@]} paired FASTQ files"

# Function to create subset using seqtk
create_subset_seqtk() {
    local input_files=("$@")
    local output_file="$1"
    shift
    
    echo "Creating subset using seqtk..."
    
    # Combine all input files and take first N reads
    {
        for file in "${input_files[@]}"; do
            seqtk seq "$file"
        done
    } | seqtk seq -A - | head -n $((SUBSET_READS * 4)) | gzip > "$output_file.tmp"
    
    mv "$output_file.tmp" "$output_file"
}

# Function to create subset using zcat/head fallback
create_subset_fallback() {
    local input_files=("$@")
    local output_file="$1"
    shift
    
    echo "Creating subset using zcat/head fallback..."
    
    # Calculate lines needed (4 lines per read)
    local lines_needed=$((SUBSET_READS * 4))
    
    # Combine all input files and take first N lines
    {
        for file in "${input_files[@]}"; do
            zcat "$file"
        done
    } | head -n $lines_needed | gzip > "$output_file.tmp"
    
    mv "$output_file.tmp" "$output_file"
}

echo
echo "Creating R1 subset..."
if [[ "$SEQTK_AVAILABLE" == "true" ]]; then
    create_subset_seqtk "$SUBSET_DIR/R1.fastq.gz" "${R1_FILES[@]}"
else
    create_subset_fallback "$SUBSET_DIR/R1.fastq.gz" "${R1_FILES[@]}"
fi

echo "Creating R2 subset..."
if [[ "$SEQTK_AVAILABLE" == "true" ]]; then
    create_subset_seqtk "$SUBSET_DIR/R2.fastq.gz" "${R2_FILES[@]}"
else
    create_subset_fallback "$SUBSET_DIR/R2.fastq.gz" "${R2_FILES[@]}"
fi

# Verify subset files
echo
echo "Verifying subset files..."

for file in "$SUBSET_DIR/R1.fastq.gz" "$SUBSET_DIR/R2.fastq.gz"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Failed to create $file" >&2
        exit 1
    fi
    
    # Count reads in subset
    read_count=$(zcat "$file" | wc -l)
    read_count=$((read_count / 4))
    
    size=$(stat -c%s "$file")
    echo "✓ $file: $read_count reads ($(numfmt --to=iec $size))"
    
    # Verify we have reasonable number of reads
    if [[ $read_count -lt $((SUBSET_READS / 2)) ]]; then
        echo "WARNING: Subset has fewer reads than expected ($read_count < $((SUBSET_READS / 2)))" >&2
    fi
done

# Create subset manifest
echo
echo "Creating subset manifest..."
{
    echo "# Subset FASTQ manifest generated on $(date)"
    echo "# Subset size: $SUBSET_READS reads per file"
    echo "# R1_FILE\tR2_FILE"
    echo -e "$SUBSET_DIR/R1.fastq.gz\t$SUBSET_DIR/R2.fastq.gz"
} > "$MANIFEST_DIR/fastq_subset.list"

echo "✓ Subset manifest: $MANIFEST_DIR/fastq_subset.list"

echo
echo "=== Stage 0 Complete: Subset FASTQs created ==="
echo "Subset files available at: $SUBSET_DIR/"
echo "Use --subset flag in downstream scripts to test with these smaller files"
