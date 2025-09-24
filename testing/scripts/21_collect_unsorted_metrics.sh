#!/bin/bash
# Stage 2.1 - Collect unsorted metrics
# Purpose: gather metrics analogous to Stage 1

set -euo pipefail

# Configuration
INPUT_DIR="$(dirname "$0")/../output/unsorted_twopass"
OUTPUT_FILE="$INPUT_DIR/metrics.json"

echo "=== STARsolo Unsorted CB/UB Testing - Stage 2.1: Collect Unsorted Metrics ==="
echo "Timestamp: $(date)"
echo "Input directory: $INPUT_DIR"
echo

# Check if unsorted run exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Unsorted directory not found: $INPUT_DIR" >&2
    echo "Run 20_run_unsorted_two_pass.sh first" >&2
    exit 1
fi

# Check for required files
REQUIRED_FILES=(
    "Aligned.out.bam"
    "Log.final.out"
    "Solo.out/Gene/Features.tsv"
    "Solo.out/Gene/Barcodes.tsv"
    "Solo.out/Gene/matrix.mtx"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "$INPUT_DIR/$file" ]]; then
        echo "ERROR: Required file not found: $INPUT_DIR/$file" >&2
        exit 1
    fi
done

echo "Collecting metrics..."

# Initialize metrics JSON
cat > "$OUTPUT_FILE" << 'EOF'
{
  "collection_timestamp": "",
  "collection_hostname": "",
  "run_mode": "two_pass_unsorted",
  "bam_metrics": {},
  "star_log_metrics": {},
  "solo_metrics": {},
  "file_sizes": {},
  "checksums": {},
  "temp_file_info": {}
}
EOF

# Helper function to update JSON
update_json() {
    local key="$1"
    local value="$2"
    local temp_file=$(mktemp)
    
    if [[ "$value" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        # Numeric value
        jq --arg k "$key" --argjson v "$value" '. | setpath($k | split("."); $v)' "$OUTPUT_FILE" > "$temp_file"
    else
        # String value
        jq --arg k "$key" --arg v "$value" '. | setpath($k | split("."); $v)' "$OUTPUT_FILE" > "$temp_file"
    fi
    
    mv "$temp_file" "$OUTPUT_FILE"
}

# Collection metadata
update_json "collection_timestamp" "$(date -Iseconds)"
update_json "collection_hostname" "$(hostname)"

# Check for temp file information
TMP_FILE="$INPUT_DIR/Aligned.out.unsorted.solo.tmp"
if [[ -f "$TMP_FILE" ]]; then
    size=$(stat -c%s "$TMP_FILE")
    update_json "temp_file_info.exists" "true"
    update_json "temp_file_info.size" "$size"
    update_json "temp_file_info.path" "$TMP_FILE"
    echo "✓ Temp file found: $(numfmt --to=iec $size)"
else
    update_json "temp_file_info.exists" "false"
    echo "✓ Temp file properly cleaned up"
fi

# BAM metrics using samtools if available
if command -v samtools >/dev/null 2>&1; then
    echo "Collecting BAM metrics..."
    
    BAM_FILE="$INPUT_DIR/Aligned.out.bam"
    
    # Basic read counts
    total_reads=$(samtools view -c "$BAM_FILE")
    mapped_reads=$(samtools view -c -F 4 "$BAM_FILE")
    unmapped_reads=$((total_reads - mapped_reads))
    
    update_json "bam_metrics.total_reads" "$total_reads"
    update_json "bam_metrics.mapped_reads" "$mapped_reads"
    update_json "bam_metrics.unmapped_reads" "$unmapped_reads"
    
    if [[ $total_reads -gt 0 ]]; then
        mapping_rate=$(echo "scale=4; $mapped_reads * 100 / $total_reads" | bc)
        update_json "bam_metrics.mapping_rate_percent" "$mapping_rate"
    fi
    
    # Sample CB/UB tag counts from first 10000 reads
    echo "Sampling CB/UB tags..."
    sample_reads=$(samtools view "$BAM_FILE" | head -10000)
    cb_tags=$(echo "$sample_reads" | grep -c "CB:Z:" || true)
    ub_tags=$(echo "$sample_reads" | grep -c "UB:Z:" || true)
    cr_tags=$(echo "$sample_reads" | grep -c "CR:Z:" || true)
    ur_tags=$(echo "$sample_reads" | grep -c "UR:Z:" || true)
    gx_tags=$(echo "$sample_reads" | grep -c "GX:Z:" || true)
    
    update_json "bam_metrics.sample_cb_tags" "$cb_tags"
    update_json "bam_metrics.sample_ub_tags" "$ub_tags"
    update_json "bam_metrics.sample_cr_tags" "$cr_tags"
    update_json "bam_metrics.sample_ur_tags" "$ur_tags"
    update_json "bam_metrics.sample_gx_tags" "$gx_tags"
    
    # Flag statistics
    primary_alignments=$(samtools view -c -F 256 "$BAM_FILE")
    secondary_alignments=$(samtools view -c -f 256 "$BAM_FILE")
    
    update_json "bam_metrics.primary_alignments" "$primary_alignments"
    update_json "bam_metrics.secondary_alignments" "$secondary_alignments"
    
    # Check if BAM is actually unsorted by sampling positions
    echo "Checking BAM sort order..."
    positions=$(samtools view "$BAM_FILE" | head -1000 | awk '{print $4}' | head -100)
    sorted_positions=$(echo "$positions" | sort -n)
    
    if [[ "$positions" == "$sorted_positions" ]]; then
        update_json "bam_metrics.appears_sorted" "true"
        echo "! BAM appears sorted in sample"
    else
        update_json "bam_metrics.appears_sorted" "false"
        echo "✓ BAM appears unsorted"
    fi
    
    echo "✓ BAM metrics collected"
else
    echo "! samtools not available, skipping BAM metrics"
fi

# STAR log metrics (identical logic to sorted version)
echo "Collecting STAR log metrics..."
LOG_FILE="$INPUT_DIR/Log.final.out"

if [[ -f "$LOG_FILE" ]]; then
    # Extract key metrics from log
    while IFS= read -r line; do
        case "$line" in
            *"Number of input reads"*)
                value=$(echo "$line" | awk '{print $NF}')
                update_json "star_log_metrics.input_reads" "$value"
                ;;
            *"Uniquely mapped reads number"*)
                value=$(echo "$line" | awk '{print $NF}')
                update_json "star_log_metrics.uniquely_mapped_reads" "$value"
                ;;
            *"Uniquely mapped reads %"*)
                value=$(echo "$line" | awk '{print $NF}' | tr -d '%')
                update_json "star_log_metrics.uniquely_mapped_percent" "$value"
                ;;
            *"Number of reads mapped to multiple loci"*)
                value=$(echo "$line" | awk '{print $NF}')
                update_json "star_log_metrics.multi_mapped_reads" "$value"
                ;;
            *"% of reads mapped to multiple loci"*)
                value=$(echo "$line" | awk '{print $NF}' | tr -d '%')
                update_json "star_log_metrics.multi_mapped_percent" "$value"
                ;;
            *"Number of reads unmapped"*)
                if [[ "$line" == *"too short"* ]]; then
                    value=$(echo "$line" | awk '{print $NF}')
                    update_json "star_log_metrics.unmapped_too_short" "$value"
                elif [[ "$line" == *"other"* ]]; then
                    value=$(echo "$line" | awk '{print $NF}')
                    update_json "star_log_metrics.unmapped_other" "$value"
                fi
                ;;
        esac
    done < "$LOG_FILE"
    
    echo "✓ STAR log metrics collected"
fi

# Solo metrics (identical logic to sorted version)
echo "Collecting Solo metrics..."

# Gene matrix dimensions
GENE_MATRIX="$INPUT_DIR/Solo.out/Gene/matrix.mtx"
if [[ -f "$GENE_MATRIX" ]]; then
    # Read matrix dimensions from MTX header
    dims=$(grep -v '^%' "$GENE_MATRIX" | head -1)
    genes=$(echo "$dims" | awk '{print $1}')
    barcodes=$(echo "$dims" | awk '{print $2}')
    entries=$(echo "$dims" | awk '{print $3}')
    
    update_json "solo_metrics.gene_matrix.genes" "$genes"
    update_json "solo_metrics.gene_matrix.barcodes" "$barcodes"
    update_json "solo_metrics.gene_matrix.entries" "$entries"
fi

# Filtered matrix dimensions
FILTERED_MATRIX="$INPUT_DIR/Solo.out/Gene/Filtered/matrix.mtx"
if [[ -f "$FILTERED_MATRIX" ]]; then
    dims=$(grep -v '^%' "$FILTERED_MATRIX" | head -1)
    genes=$(echo "$dims" | awk '{print $1}')
    barcodes=$(echo "$dims" | awk '{print $2}')
    entries=$(echo "$dims" | awk '{print $3}')
    
    update_json "solo_metrics.filtered_matrix.genes" "$genes"
    update_json "solo_metrics.filtered_matrix.barcodes" "$barcodes"
    update_json "solo_metrics.filtered_matrix.entries" "$entries"
fi

# Barcode counts
BARCODES_FILE="$INPUT_DIR/Solo.out/Gene/Barcodes.tsv"
if [[ -f "$BARCODES_FILE" ]]; then
    total_barcodes=$(wc -l < "$BARCODES_FILE")
    update_json "solo_metrics.total_barcodes" "$total_barcodes"
fi

FILTERED_BARCODES_FILE="$INPUT_DIR/Solo.out/Gene/Filtered/barcodes.tsv"
if [[ -f "$FILTERED_BARCODES_FILE" ]]; then
    filtered_barcodes=$(wc -l < "$FILTERED_BARCODES_FILE")
    update_json "solo_metrics.filtered_barcodes" "$filtered_barcodes"
fi

# Feature counts
FEATURES_FILE="$INPUT_DIR/Solo.out/Gene/Features.tsv"
if [[ -f "$FEATURES_FILE" ]]; then
    total_features=$(wc -l < "$FEATURES_FILE")
    update_json "solo_metrics.total_features" "$total_features"
fi

echo "✓ Solo metrics collected"

# File sizes
echo "Collecting file sizes..."
for file in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$INPUT_DIR/$file" ]]; then
        size=$(stat -c%s "$INPUT_DIR/$file")
        # Replace / with . for JSON key
        key=$(echo "$file" | tr '/' '.')
        update_json "file_sizes.$key" "$size"
    fi
done

# Checksums
echo "Collecting checksums..."
CHECKSUM_FILE="$INPUT_DIR/checksums.md5"
if [[ -f "$CHECKSUM_FILE" ]]; then
    while read -r checksum filename; do
        if [[ -n "$checksum" && -n "$filename" && ! "$filename" =~ ^# ]]; then
            # Replace / with . for JSON key
            key=$(echo "$filename" | tr '/' '.')
            update_json "checksums.$key" "$checksum"
        fi
    done < "$CHECKSUM_FILE"
fi

echo "✓ File sizes and checksums collected"

# Validate JSON
if command -v jq >/dev/null 2>&1; then
    if jq empty "$OUTPUT_FILE" 2>/dev/null; then
        echo "✓ Metrics JSON validated"
    else
        echo "ERROR: Invalid JSON generated" >&2
        exit 1
    fi
else
    echo "! jq not available, skipping JSON validation"
fi

# Summary
echo
echo "=== Stage 2.1 Complete: Unsorted Metrics Collected ==="
echo "Metrics file: $OUTPUT_FILE"

# Display key metrics
if command -v jq >/dev/null 2>&1; then
    echo
    echo "Key metrics summary:"
    echo "  Total reads: $(jq -r '.bam_metrics.total_reads // "N/A"' "$OUTPUT_FILE")"
    echo "  Mapped reads: $(jq -r '.bam_metrics.mapped_reads // "N/A"' "$OUTPUT_FILE")"
    echo "  Mapping rate: $(jq -r '.bam_metrics.mapping_rate_percent // "N/A"' "$OUTPUT_FILE")%"
    echo "  CB tags (sample): $(jq -r '.bam_metrics.sample_cb_tags // "N/A"' "$OUTPUT_FILE")"
    echo "  UB tags (sample): $(jq -r '.bam_metrics.sample_ub_tags // "N/A"' "$OUTPUT_FILE")"
    echo "  Total barcodes: $(jq -r '.solo_metrics.total_barcodes // "N/A"' "$OUTPUT_FILE")"
    echo "  Filtered barcodes: $(jq -r '.solo_metrics.filtered_barcodes // "N/A"' "$OUTPUT_FILE")"
    echo "  Temp file preserved: $(jq -r '.temp_file_info.exists // "N/A"' "$OUTPUT_FILE")"
fi

echo "Ready for cross-run comparisons (Stage 3)"
