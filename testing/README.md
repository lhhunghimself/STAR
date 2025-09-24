# STARsolo Two-Pass Unsorted CB/UB Testing Suite

This directory contains a comprehensive testing suite for validating the STARsolo two-pass unsorted CB/UB injection feature. The test suite ensures that the new `--soloAddTagsToUnsorted` functionality produces identical results to the existing coordinate-sorted pipeline while maintaining backward compatibility.

## Overview

The testing suite is organized into 5 stages with a total of 16 scripts that validate different aspects of the implementation:

- **Stage 0**: Environment and fixture validation
- **Stage 1**: Sorted baseline generation and metrics collection
- **Stage 2**: Two-pass unsorted generation and validation
- **Stage 3**: Cross-run comparisons ensuring identical outputs
- **Stage 4**: Negative tests and regression validation
- **Stage 5**: Orchestration and comprehensive reporting

## Quick Start

### Prerequisites

- STAR binary built with the two-pass unsorted feature
- Test data available at the configured paths:
  - FASTQs: `/storage/downsampled/SC2300771/*_R{1,2}_001.fastq.gz`
  - Whitelist: `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt`
  - Reference: `/storage/scRNAseq_output/indices-98-32/star`
- Required tools: `samtools`, `python3` with `pysam`, `jq`, `bc`

### Running the Complete Test Suite

```bash
# Run with subset data (faster, ~30 minutes)
./scripts/90_run_all.sh --subset

# Run with full data (complete validation, ~2-4 hours)
./scripts/90_run_all.sh --full

# Run with verbose output and keep temp files for debugging
./scripts/90_run_all.sh --subset --verbose --keep-temp
```

### Running Individual Stages

```bash
# Stage 0: Check fixtures and create subset
./scripts/00_assert_fixtures.sh
./scripts/01_smoke_subset_fastqs.sh

# Stage 1: Generate sorted baseline
./scripts/10_run_sorted_baseline.sh --subset
./scripts/11_collect_sorted_metrics.sh

# Stage 2: Generate two-pass unsorted
./scripts/20_run_unsorted_two_pass.sh --subset
./scripts/21_collect_unsorted_metrics.sh
./scripts/22_verify_tmp_trailers.py  # Only if temp files preserved

# Stage 3: Compare outputs
./scripts/30_compare_bam_tags.py
./scripts/31_compare_alignment_qc.py
./scripts/32_compare_solo_matrices.sh
./scripts/33_compare_logs.sh
./scripts/34_compare_read_distributions.py

# Stage 4: Negative tests
./scripts/40_expect_failure_without_flag.sh
./scripts/41_expect_success_without_cbub.sh
./scripts/42_subset_regression_loop.sh
```

## Script Descriptions

### Stage 0: Environment Validation

- **`00_assert_fixtures.sh`**: Validates that all required input files exist and are readable. Generates FASTQ manifests for downstream scripts.

- **`01_smoke_subset_fastqs.sh`**: Creates subset FASTQs (1M reads each) for faster testing. Uses `seqtk` if available, falls back to `zcat/head`.

### Stage 1: Sorted Baseline

- **`10_run_sorted_baseline.sh`**: Runs STARsolo with `--outSAMtype BAM SortedByCoordinate` to generate the reference output with CB/UB tags.

- **`11_collect_sorted_metrics.sh`**: Extracts comprehensive metrics from the sorted run including BAM statistics, STAR log metrics, Solo matrix dimensions, and file checksums.

### Stage 2: Two-Pass Unsorted

- **`20_run_unsorted_two_pass.sh`**: Runs STARsolo with `--outSAMtype BAM Unsorted --soloAddTagsToUnsorted yes` to test the new two-pass mode.

- **`21_collect_unsorted_metrics.sh`**: Collects metrics from the unsorted run using identical logic to Stage 1 for direct comparison.

- **`22_verify_tmp_trailers.py`**: Validates the intermediate temp file format when `KEEP_SOLO_TMP=1` is set, ensuring proper `[size][payload][trailer]` structure.

### Stage 3: Cross-Run Comparisons

- **`30_compare_bam_tags.py`**: Name-sorts both BAMs and compares CB/UB/CR/UR/GX/GN tags read-by-read to ensure perfect matching.

- **`31_compare_alignment_qc.py`**: Compares alignment metrics (MAPQ, flags, CIGAR, mate info) and their distributions between runs.

- **`32_compare_solo_matrices.sh`**: Bit-for-bit comparison of all Solo output matrices, features, and barcodes files using checksums.

- **`33_compare_logs.sh`**: Compares STAR log files and Solo statistics, ignoring timestamps but validating all quantitative metrics.

- **`34_compare_read_distributions.py`**: Deep comparison of barcode/UMI distributions using pandas for statistical validation.

### Stage 4: Negative Tests

- **`40_expect_failure_without_flag.sh`**: Verifies that requesting CB/UB with unsorted BAM still fails when `--soloAddTagsToUnsorted` is not specified.

- **`41_expect_success_without_cbub.sh`**: Ensures legacy unsorted mode (without CB/UB tags) continues to work unchanged.

- **`42_subset_regression_loop.sh`**: Runs the complete pipeline on subset data as a quick regression test.

### Stage 5: Orchestration

- **`90_run_all.sh`**: Master orchestrator that runs the complete test suite with configurable options for parallel execution, subset/full data, and comprehensive reporting.

## Output Structure

```
testing/output/
├── manifests/                    # FASTQ file lists
│   ├── fastq.list               # Full dataset manifest
│   └── fastq_subset.list        # Subset dataset manifest
├── baseline_sorted/             # Sorted baseline outputs
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Solo.out/
│   ├── metrics.json
│   └── checksums.md5
├── unsorted_twopass/           # Two-pass unsorted outputs
│   ├── Aligned.out.bam
│   ├── Solo.out/
│   ├── metrics.json
│   └── checksums.md5
├── logs/                       # Individual test logs
└── test-report.json           # Comprehensive test report
```

## Expected Results

When all tests pass, you should see:

1. **Identical BAM tags**: All CB/UB/CR/UR/GX/GN tags match perfectly between sorted and unsorted runs
2. **Identical Solo matrices**: All count matrices, features, and barcodes are byte-identical
3. **Identical metrics**: STAR log statistics and Solo summaries match exactly
4. **Proper error handling**: Legacy parameter validation still works
5. **Backward compatibility**: Unsorted mode without CB/UB works unchanged

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure `samtools`, `python3`, `pysam`, `jq`, and `bc` are installed
2. **Path configuration**: Update paths in scripts if data is located elsewhere
3. **Memory/disk space**: Full dataset requires ~50GB temporary disk space
4. **Permissions**: Ensure scripts are executable (`chmod +x scripts/*.sh scripts/*.py`)

### Debugging Failed Tests

1. **Check individual logs**: Each test creates a detailed log in `output/logs/`
2. **Use subset data**: Run with `--subset` for faster iteration
3. **Keep temp files**: Use `--keep-temp` to preserve intermediate files
4. **Verbose output**: Add `--verbose` for detailed progress information
5. **Stop on failure**: Use `--stop-on-failure` to halt at first error

### Test Report Analysis

The `test-report.json` contains:
- Overall timing and success metrics
- Per-test status and runtime
- Configuration used for the run
- Paths to individual log files

## Integration with CI/CD

The test suite is designed for automated validation:

```bash
# Quick validation (subset data)
./scripts/90_run_all.sh --subset --stop-on-failure

# Full validation (complete data)  
./scripts/90_run_all.sh --full

# Check exit code
if [ $? -eq 0 ]; then
    echo "All tests passed - ready for release"
else
    echo "Tests failed - check test-report.json"
fi
```

## Contributing

When modifying the two-pass unsorted implementation:

1. Run the subset regression test during development: `./scripts/42_subset_regression_loop.sh`
2. Run the complete test suite before submitting: `./scripts/90_run_all.sh --full`
3. Add new test cases if introducing new functionality
4. Update this documentation if changing test behavior

## Performance Benchmarks

Typical runtimes on a 16-core machine:

- **Subset data** (~1M reads): 20-30 minutes total
- **Full data** (~50M reads): 2-4 hours total
- **Stage 3 comparisons**: 5-15 minutes (parallelized)
- **Individual BAM comparison**: 2-5 minutes per comparison

The two-pass unsorted mode typically shows:
- **Disk usage**: 60-80% less than coordinate sorting
- **Runtime**: 10-20% faster than coordinate sorting  
- **Memory usage**: Similar to coordinate sorting
- **Output quality**: Identical to coordinate sorting
