#!/bin/bash
# Stage 0 - Assert fixtures
# Purpose: fail fast if any required input is missing or unreadable

set -euo pipefail

# Configuration
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
REFERENCE_DIR="/storage/scRNAseq_output/indices-98-32/star"
FASTQ_DIR="/storage/downsampled/SC2300771"
OUTPUT_DIR="$(dirname "$0")/../output"
MANIFEST_DIR="$OUTPUT_DIR/manifests"

echo "=== STARsolo Unsorted CB/UB Testing - Stage 0: Assert Fixtures ==="
echo "Timestamp: $(date)"
echo

# Create output directories if they don't exist
mkdir -p "$MANIFEST_DIR"

# Function to check file/directory existence
check_exists() {
    local path="$1"
    local description="$2"
    
    if [[ ! -e "$path" ]]; then
        echo "ERROR: $description not found at: $path" >&2
        return 1
    fi
    
    if [[ ! -r "$path" ]]; then
        echo "ERROR: $description not readable at: $path" >&2
        return 1
    fi
    
    echo "✓ $description found: $path"
}

# Function to check directory contents
check_dir_contents() {
    local path="$1"
    local description="$2"
    local pattern="$3"
    
    if [[ ! -d "$path" ]]; then
        echo "ERROR: $description directory not found: $path" >&2
        return 1
    fi
    
    local count=$(find "$path" -name "$pattern" -type f 2>/dev/null | wc -l)
    if [[ $count -eq 0 ]]; then
        echo "ERROR: No files matching $pattern found in $description: $path" >&2
        return 1
    fi
    
    echo "✓ $description contains $count files matching $pattern"
}

echo "Checking required fixtures..."

# Check whitelist
check_exists "$WHITELIST" "Whitelist file"

# Check reference directory and key files
check_exists "$REFERENCE_DIR" "Reference directory"
check_exists "$REFERENCE_DIR/genomeParameters.txt" "Genome parameters"
check_exists "$REFERENCE_DIR/Genome" "Genome index"

# Check FASTQ directory and files
check_exists "$FASTQ_DIR" "FASTQ directory"
check_dir_contents "$FASTQ_DIR" "FASTQ directory" "*_R1_001.fastq.gz"
check_dir_contents "$FASTQ_DIR" "FASTQ directory" "*_R2_001.fastq.gz"

echo
echo "Generating FASTQ manifest..."

# Create manifest of R1/R2 pairs
R1_FILES=($(find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" -type f | sort))
R2_FILES=($(find "$FASTQ_DIR" -name "*_R2_001.fastq.gz" -type f | sort))

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 files found" >&2
    exit 1
fi

if [[ ${#R2_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R2 files found" >&2
    exit 1
fi

if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
    echo "ERROR: Mismatch between R1 (${#R1_FILES[@]}) and R2 (${#R2_FILES[@]}) file counts" >&2
    exit 1
fi

# Write manifest files
echo "${R1_FILES[@]}" | tr ' ' '\n' > "$MANIFEST_DIR/R1.list"
echo "${R2_FILES[@]}" | tr ' ' '\n' > "$MANIFEST_DIR/R2.list"

# Create combined manifest
{
    echo "# FASTQ manifest generated on $(date)"
    echo "# R1_FILE\tR2_FILE"
    for i in "${!R1_FILES[@]}"; do
        echo -e "${R1_FILES[$i]}\t${R2_FILES[$i]}"
    done
} > "$MANIFEST_DIR/fastq.list"

echo "✓ Generated manifest with ${#R1_FILES[@]} paired files"
echo "  - R1 list: $MANIFEST_DIR/R1.list"
echo "  - R2 list: $MANIFEST_DIR/R2.list"
echo "  - Combined manifest: $MANIFEST_DIR/fastq.list"

# Verify paired files have matching lane identifiers
echo
echo "Verifying paired file consistency..."
for i in "${!R1_FILES[@]}"; do
    r1_base=$(basename "${R1_FILES[$i]}" _R1_001.fastq.gz)
    r2_base=$(basename "${R2_FILES[$i]}" _R2_001.fastq.gz)
    
    if [[ "$r1_base" != "$r2_base" ]]; then
        echo "ERROR: Mismatched pair at index $i:" >&2
        echo "  R1: ${R1_FILES[$i]}" >&2
        echo "  R2: ${R2_FILES[$i]}" >&2
        exit 1
    fi
done

echo "✓ All ${#R1_FILES[@]} pairs have matching lane identifiers"

# Check file sizes (basic sanity)
echo
echo "Checking file sizes..."
total_size=0
for file in "${R1_FILES[@]}" "${R2_FILES[@]}"; do
    size=$(stat -c%s "$file")
    if [[ $size -lt 1000 ]]; then
        echo "WARNING: File $file is very small ($size bytes)" >&2
    fi
    total_size=$((total_size + size))
done

echo "✓ Total FASTQ size: $(numfmt --to=iec $total_size)"

echo
echo "=== Stage 0 Complete: All fixtures verified ==="
echo "Ready to proceed with testing pipeline"
