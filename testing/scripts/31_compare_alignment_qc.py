#!/usr/bin/env python3
"""
Stage 3.1 - Compare alignment QC
Purpose: verify MAPQ, mate info, and key SAM flags remain unchanged
"""

import os
import sys
import argparse
import tempfile
import subprocess
from pathlib import Path
from collections import defaultdict, Counter
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare alignment QC metrics between sorted and unsorted runs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--baseline-dir",
        default=Path(__file__).parent / "../output/baseline_sorted",
        type=Path,
        help="Directory containing sorted baseline BAM"
    )
    parser.add_argument(
        "--unsorted-dir",
        default=Path(__file__).parent / "../output/unsorted_twopass",
        type=Path,
        help="Directory containing unsorted BAM"
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=10,
        help="Maximum number of mismatches to report in detail"
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=100000,
        help="Sample N reads for comparison (0 = all reads)"
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary sorted BAM files"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    return parser.parse_args()


def sort_bam_by_name(input_bam, output_bam, threads=4):
    """Sort BAM file by read name using samtools"""
    cmd = [
        "samtools", "sort", "-n",
        "-@", str(threads),
        "-o", str(output_bam),
        str(input_bam)
    ]
    
    print(f"Sorting BAM: {input_bam.name}")
    if subprocess.run(cmd, capture_output=True).returncode != 0:
        raise RuntimeError(f"Failed to sort BAM: {input_bam}")
    
    print(f"✓ Sorted: {output_bam}")


def extract_alignment_info(read):
    """Extract alignment QC information from a read"""
    return {
        'flag': read.flag,
        'reference_id': read.reference_id,
        'reference_start': read.reference_start,
        'mapping_quality': read.mapping_quality,
        'cigar': read.cigarstring,
        'next_reference_id': read.next_reference_id,
        'next_reference_start': read.next_reference_start,
        'template_length': read.template_length,
        'is_paired': read.is_paired,
        'is_proper_pair': read.is_proper_pair,
        'is_unmapped': read.is_unmapped,
        'mate_is_unmapped': read.mate_is_unmapped,
        'is_reverse': read.is_reverse,
        'mate_is_reverse': read.mate_is_reverse,
        'is_read1': read.is_read1,
        'is_read2': read.is_read2,
        'is_secondary': read.is_secondary,
        'is_qcfail': read.is_qcfail,
        'is_duplicate': read.is_duplicate,
        'is_supplementary': read.is_supplementary
    }


def compare_alignment_info(info1, info2, read_name):
    """Compare alignment info between two reads"""
    differences = []
    
    critical_fields = [
        'flag', 'reference_id', 'reference_start', 'mapping_quality',
        'cigar', 'next_reference_id', 'next_reference_start', 'template_length'
    ]
    
    flag_fields = [
        'is_paired', 'is_proper_pair', 'is_unmapped', 'mate_is_unmapped',
        'is_reverse', 'mate_is_reverse', 'is_read1', 'is_read2',
        'is_secondary', 'is_qcfail', 'is_duplicate', 'is_supplementary'
    ]
    
    for field in critical_fields:
        if info1[field] != info2[field]:
            differences.append(f"{field}: {info1[field]} vs {info2[field]}")
    
    for field in flag_fields:
        if info1[field] != info2[field]:
            differences.append(f"{field}: {info1[field]} vs {info2[field]}")
    
    return differences


def analyze_distributions(bam_file, sample_size=100000):
    """Analyze distributions of alignment metrics"""
    print(f"Analyzing alignment distributions: {bam_file.name}")
    
    stats = {
        'mapq_dist': Counter(),
        'flag_dist': Counter(),
        'cigar_dist': Counter(),
        'total_reads': 0,
        'mapped_reads': 0,
        'paired_reads': 0,
        'proper_pairs': 0
    }
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i, read in enumerate(bam):
            if sample_size > 0 and i >= sample_size:
                break
            
            stats['total_reads'] += 1
            stats['mapq_dist'][read.mapping_quality] += 1
            stats['flag_dist'][read.flag] += 1
            
            if read.cigarstring:
                stats['cigar_dist'][read.cigarstring] += 1
            
            if not read.is_unmapped:
                stats['mapped_reads'] += 1
            
            if read.is_paired:
                stats['paired_reads'] += 1
            
            if read.is_proper_pair:
                stats['proper_pairs'] += 1
            
            if i % 50000 == 0 and i > 0:
                print(f"  Analyzed {i:,} reads...")
    
    return stats


def compare_distributions(stats1, stats2, max_diffs=20):
    """Compare alignment statistics distributions"""
    print("Comparing alignment distributions...")
    
    differences = []
    
    # Compare basic counts
    fields = ['total_reads', 'mapped_reads', 'paired_reads', 'proper_pairs']
    for field in fields:
        val1, val2 = stats1[field], stats2[field]
        if val1 != val2:
            diff_pct = 100 * abs(val1 - val2) / max(val1, val2) if max(val1, val2) > 0 else 0
            differences.append(f"{field}: {val1} vs {val2} (diff: {diff_pct:.2f}%)")
    
    # Compare MAPQ distribution
    all_mapq = set(stats1['mapq_dist'].keys()) | set(stats2['mapq_dist'].keys())
    mapq_diffs = []
    for mapq in sorted(all_mapq):
        count1 = stats1['mapq_dist'][mapq]
        count2 = stats2['mapq_dist'][mapq]
        if count1 != count2:
            mapq_diffs.append(f"MAPQ {mapq}: {count1} vs {count2}")
    
    if mapq_diffs:
        differences.append(f"MAPQ distribution differences ({len(mapq_diffs)} values)")
        differences.extend(mapq_diffs[:max_diffs])
        if len(mapq_diffs) > max_diffs:
            differences.append(f"... ({len(mapq_diffs) - max_diffs} more MAPQ differences)")
    
    # Compare FLAG distribution
    all_flags = set(stats1['flag_dist'].keys()) | set(stats2['flag_dist'].keys())
    flag_diffs = []
    for flag in sorted(all_flags):
        count1 = stats1['flag_dist'][flag]
        count2 = stats2['flag_dist'][flag]
        if count1 != count2:
            flag_diffs.append(f"FLAG {flag}: {count1} vs {count2}")
    
    if flag_diffs:
        differences.append(f"FLAG distribution differences ({len(flag_diffs)} values)")
        differences.extend(flag_diffs[:max_diffs])
        if len(flag_diffs) > max_diffs:
            differences.append(f"... ({len(flag_diffs) - max_diffs} more FLAG differences)")
    
    return differences


def compare_alignment_qc(sorted_bam, unsorted_bam, max_mismatches=10, sample_size=100000, verbose=False):
    """Compare alignment QC between two name-sorted BAM files"""
    print(f"Comparing alignment QC:")
    print(f"  Baseline (sorted): {sorted_bam}")
    print(f"  Test (unsorted): {unsorted_bam}")
    print(f"  Sample size: {sample_size if sample_size > 0 else 'all reads'}")
    print()
    
    # First, analyze distributions
    stats1 = analyze_distributions(sorted_bam, sample_size)
    stats2 = analyze_distributions(unsorted_bam, sample_size)
    
    dist_diffs = compare_distributions(stats1, stats2)
    
    # Then do detailed read-by-read comparison
    total_reads = 0
    mismatched_reads = 0
    mismatches_reported = 0
    
    field_stats = defaultdict(int)
    
    print("Performing read-by-read comparison...")
    
    with pysam.AlignmentFile(sorted_bam, "rb") as bam1, \
         pysam.AlignmentFile(unsorted_bam, "rb") as bam2:
        
        iter1 = iter(bam1)
        iter2 = iter(bam2)
        
        while True:
            try:
                read1 = next(iter1)
                read2 = next(iter2)
            except StopIteration:
                break
            
            total_reads += 1
            
            # Sample reads if requested
            if sample_size > 0 and total_reads > sample_size:
                break
            
            # Progress reporting
            if total_reads % 50000 == 0:
                print(f"Compared {total_reads:,} reads...")
            
            # Compare read names
            if read1.query_name != read2.query_name:
                print(f"ERROR: Read name mismatch at position {total_reads}")
                return False
            
            # Compare alignment info
            info1 = extract_alignment_info(read1)
            info2 = extract_alignment_info(read2)
            
            differences = compare_alignment_info(info1, info2, read1.query_name)
            
            if differences:
                mismatched_reads += 1
                
                # Track field-specific differences
                for diff in differences:
                    field = diff.split(':')[0]
                    field_stats[field] += 1
                
                # Report detailed mismatches
                if mismatches_reported < max_mismatches:
                    print(f"ALIGNMENT MISMATCH {mismatched_reads}: {read1.query_name}")
                    for diff in differences:
                        print(f"  {diff}")
                    print()
                    mismatches_reported += 1
                elif mismatches_reported == max_mismatches:
                    print(f"... (suppressing further detailed reports)")
                    mismatches_reported += 1
    
    # Summary
    print(f"Read-by-read comparison complete:")
    print(f"  Total reads compared: {total_reads:,}")
    print(f"  Mismatched reads: {mismatched_reads:,}")
    if total_reads > 0:
        print(f"  Mismatch rate: {100 * mismatched_reads / total_reads:.6f}%")
    print()
    
    # Field-specific statistics
    if field_stats:
        print("Field-specific differences:")
        for field in sorted(field_stats.keys()):
            count = field_stats[field]
            print(f"  {field}: {count:,} differences")
        print()
    
    # Distribution differences
    if dist_diffs:
        print("Distribution differences:")
        for diff in dist_diffs:
            print(f"  {diff}")
        print()
    
    success = mismatched_reads == 0 and len(dist_diffs) == 0
    return success


def main():
    args = parse_args()
    
    print("=== STARsolo Unsorted CB/UB Testing - Stage 3.1: Compare Alignment QC ===")
    print(f"Timestamp: {os.popen('date').read().strip()}")
    print()
    
    # Check input files
    baseline_bam = args.baseline_dir / "Aligned.sortedByCoord.out.bam"
    unsorted_bam = args.unsorted_dir / "Aligned.out.bam"
    
    if not baseline_bam.exists():
        print(f"ERROR: Baseline BAM not found: {baseline_bam}")
        return 1
    
    if not unsorted_bam.exists():
        print(f"ERROR: Unsorted BAM not found: {unsorted_bam}")
        return 1
    
    # Check samtools availability
    if subprocess.run(["samtools", "--version"], capture_output=True).returncode != 0:
        print("ERROR: samtools not available")
        return 1
    
    # Create temporary directory for sorted files
    with tempfile.TemporaryDirectory(prefix="star_qc_compare_") as temp_dir:
        temp_dir = Path(temp_dir)
        
        if args.keep_temp:
            # Use a persistent directory
            temp_dir = Path("testing_temp_qc_sort")
            temp_dir.mkdir(exist_ok=True)
            print(f"Using persistent temp directory: {temp_dir}")
        
        # Sort both BAMs by name
        sorted_baseline = temp_dir / "baseline_sorted_by_name.bam"
        sorted_unsorted = temp_dir / "unsorted_sorted_by_name.bam"
        
        try:
            sort_bam_by_name(baseline_bam, sorted_baseline)
            sort_bam_by_name(unsorted_bam, sorted_unsorted)
            
            # Compare alignment QC
            success = compare_alignment_qc(
                sorted_baseline, 
                sorted_unsorted,
                max_mismatches=args.max_mismatches,
                sample_size=args.sample_size,
                verbose=args.verbose
            )
            
        except Exception as e:
            print(f"ERROR: {e}")
            return 1
        
        finally:
            if not args.keep_temp:
                # Cleanup handled by context manager
                pass
            else:
                print(f"Temporary files preserved in: {temp_dir}")
    
    print()
    if success:
        print("=== Stage 3.1 Complete: Alignment QC Comparison PASSED ===")
        print("All alignment metrics match between sorted and unsorted runs")
        return 0
    else:
        print("=== Stage 3.1 FAILED: Alignment QC Comparison Failed ===")
        print("Differences found in alignment metrics")
        return 1


if __name__ == "__main__":
    sys.exit(main())
