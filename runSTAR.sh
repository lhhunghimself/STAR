#!/bin/bash

# Run STAR with CB/UB tag injection
# Based on production parameters from SC2300771

set -euo pipefail

# Configuration
NEW_STAR_BINARY="./STAR"                     # Newly built binary in current directory
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/JAX_sequences/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="/mnt/pikachu/forked/Alignments"

NEW_DIR="${OUTPUT_BASE}/${SAMPLE_ID}"
TEMP_DIR="/storage/tmp/${SAMPLE_ID}"

#clean up the temp directory
rm -rf $TEMP_DIR
mkdir -p $TEMP_DIR

# Set ulimit for STAR
ulimit -n 100000 || true

# Helper: clean dir
cleanup_dir() {
    local dir="$1"
    echo "Cleaning up: $dir"
    rm -rf "$dir" "${dir}_STARtmp" 2>/dev/null || true
    mkdir -p "$dir"
}

# Check requirements
echo "=== Checking inputs/binaries ==="
[[ -d "$BASE_DIR" ]] || { echo "ERROR: Input dir not found: $BASE_DIR"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$NEW_STAR_BINARY" ]] || { echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"; exit 1; }

# Build R2/R1 file lists (8 lanes)
R2_FILES=""; R1_FILES=""
for lane in L001 L002 L003 L004 L005 L006 L007 L008; do
    r2_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R2_001.fastq.gz"
    r1_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R1_001.fastq.gz"
    if ls $r2_glob 1>/dev/null 2>&1; then
        [[ -n "$R2_FILES" ]] && R2_FILES+="",
        R2_FILES+="$(ls $r2_glob)"
    fi
    if ls $r1_glob 1>/dev/null 2>&1; then
        [[ -n "$R1_FILES" ]] && R1_FILES+="",
        R1_FILES+="$(ls $r1_glob)"
    fi
done
[[ -n "$R2_FILES" && -n "$R1_FILES" ]] || { echo "ERROR: Could not find R1/R2 FASTQs in $BASE_DIR"; exit 1; }

echo "R2: $R2_FILES"
echo "R1: $R1_FILES"

# Common STAR params
COMMON_PARAMS=(
    --runThreadN $THREADS
    --outTmpDir $TEMP_DIR
    --soloType CB_UMI_Simple
    --soloCBlen 16
    --soloUMIlen 12
    --soloUMIstart 17
    --soloCBstart 1
    --soloBarcodeReadLength 0
    --soloCBwhitelist "$WHITELIST"
    --genomeDir "$GENOME_DIR"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outBAMcompression 6
    --soloMultiMappers Unique
    --alignIntronMax 1
    --alignMatesGapMax 0
    --outFilterMismatchNmax 10
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 16
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 1000
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloCellFilter None
    --clipAdapterType CellRanger4
    --soloFeatures Gene
    --readFilesIn "$R2_FILES" "$R1_FILES"
    --alignEndsType Local
    --readFilesCommand zcat
)

NEW_UNSORTED_TAGS_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloAddTagsToUnsorted yes
)

# Run new STAR with tags
cleanup_dir "$NEW_DIR"

echo "=== Running new STAR with CB/UB tags ==="
$NEW_STAR_BINARY \
    "${NEW_UNSORTED_TAGS_PARAMS[@]}" \
    --outFileNamePrefix "${NEW_DIR}/"

echo "Run finished"
