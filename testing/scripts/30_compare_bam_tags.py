#!/usr/bin/env python3
"""
Stage 3 - Compare BAM tags
Purpose: ensure CB/UB (and CR/UR) tags are identical per read
"""

import os
import sys
import argparse
import tempfile
import subprocess
from pathlib import Path
from collections import defaultdict
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare BAM tags between sorted and unsorted runs",
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
        default=0,
        help="Sample only N reads for faster comparison (0 = all reads)"
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
    
    # Index is not needed for name-sorted BAMs
    print(f"✓ Sorted: {output_bam}")


def extract_tags(read):
    """Extract relevant tags from a pysam read"""
    tags = {}
    for tag, value in read.get_tags():
        if tag in ['CB', 'UB', 'CR', 'UR', 'GX', 'GN']:
            tags[tag] = value
    return tags


def compare_reads(read1, read2, read_name):
    """Compare two reads and return differences"""
    differences = []
    
    # Compare flags
    if read1.flag != read2.flag:
        differences.append(f"FLAG: {read1.flag} vs {read2.flag}")
    
    # Compare mapping info (should be identical for same read)
    if read1.reference_id != read2.reference_id:
        differences.append(f"RNAME: {read1.reference_id} vs {read2.reference_id}")
    
    if read1.reference_start != read2.reference_start:
        differences.append(f"POS: {read1.reference_start} vs {read2.reference_start}")
    
    if read1.mapping_quality != read2.mapping_quality:
        differences.append(f"MAPQ: {read1.mapping_quality} vs {read2.mapping_quality}")
    
    if read1.cigarstring != read2.cigarstring:
        differences.append(f"CIGAR: {read1.cigarstring} vs {read2.cigarstring}")
    
    # Compare tags
    tags1 = extract_tags(read1)
    tags2 = extract_tags(read2)
    
    all_tags = set(tags1.keys()) | set(tags2.keys())
    for tag in sorted(all_tags):
        val1 = tags1.get(tag, "MISSING")
        val2 = tags2.get(tag, "MISSING")
        if val1 != val2:
            differences.append(f"{tag}: '{val1}' vs '{val2}'")
    
    return differences


def compare_bam_files(sorted_bam, unsorted_bam, max_mismatches=10, sample_size=0, verbose=False):
    """Compare two name-sorted BAM files"""
    print(f"Comparing BAM files:")
    print(f"  Baseline (sorted): {sorted_bam}")
    print(f"  Test (unsorted): {unsorted_bam}")
    print()
    
    total_reads = 0
    mismatched_reads = 0
    mismatches_reported = 0
    
    # Track tag statistics
    tag_stats = defaultdict(lambda: {'total': 0, 'mismatches': 0})
    
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
            if total_reads % 100000 == 0:
                print(f"Processed {total_reads:,} reads...")
            
            # Compare read names
            if read1.query_name != read2.query_name:
                print(f"ERROR: Read name mismatch at position {total_reads}")
                print(f"  Baseline: {read1.query_name}")
                print(f"  Test: {read2.query_name}")
                return False
            
            # Compare reads
            differences = compare_reads(read1, read2, read1.query_name)
            
            if differences:
                mismatched_reads += 1
                
                # Update tag statistics
                for diff in differences:
                    if ':' in diff and diff.split(':')[0] in ['CB', 'UB', 'CR', 'UR', 'GX', 'GN']:
                        tag = diff.split(':')[0]
                        tag_stats[tag]['mismatches'] += 1
                
                # Report detailed mismatches
                if mismatches_reported < max_mismatches:
                    print(f"MISMATCH {mismatched_reads}: {read1.query_name}")
                    for diff in differences:
                        print(f"  {diff}")
                    print()
                    mismatches_reported += 1
                elif mismatches_reported == max_mismatches:
                    print(f"... (suppressing further detailed reports)")
                    mismatches_reported += 1
            
            # Update tag presence statistics
            tags1 = extract_tags(read1)
            tags2 = extract_tags(read2)
            all_tags = set(tags1.keys()) | set(tags2.keys())
            for tag in all_tags:
                tag_stats[tag]['total'] += 1
        
        # Check if one file has more reads
        try:
            extra_read = next(iter1)
            print(f"ERROR: Baseline BAM has extra reads (starting with {extra_read.query_name})")
            return False
        except StopIteration:
            pass
        
        try:
            extra_read = next(iter2)
            print(f"ERROR: Test BAM has extra reads (starting with {extra_read.query_name})")
            return False
        except StopIteration:
            pass
    
    # Summary
    print(f"Comparison complete:")
    print(f"  Total reads compared: {total_reads:,}")
    print(f"  Mismatched reads: {mismatched_reads:,}")
    print(f"  Mismatch rate: {100 * mismatched_reads / total_reads:.6f}%")
    print()
    
    # Tag statistics
    if tag_stats:
        print("Tag statistics:")
        for tag in sorted(tag_stats.keys()):
            stats = tag_stats[tag]
            if stats['total'] > 0:
                mismatch_rate = 100 * stats['mismatches'] / stats['total']
                print(f"  {tag}: {stats['mismatches']:,} / {stats['total']:,} mismatches ({mismatch_rate:.6f}%)")
        print()
    
    return mismatched_reads == 0


def main():
    args = parse_args()
    
    print("=== STARsolo Unsorted CB/UB Testing - Stage 3: Compare BAM Tags ===")
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
    with tempfile.TemporaryDirectory(prefix="star_bam_compare_") as temp_dir:
        temp_dir = Path(temp_dir)
        
        if args.keep_temp:
            # Use a persistent directory
            temp_dir = Path("testing_temp_bam_sort")
            temp_dir.mkdir(exist_ok=True)
            print(f"Using persistent temp directory: {temp_dir}")
        
        # Sort both BAMs by name
        sorted_baseline = temp_dir / "baseline_sorted_by_name.bam"
        sorted_unsorted = temp_dir / "unsorted_sorted_by_name.bam"
        
        try:
            sort_bam_by_name(baseline_bam, sorted_baseline)
            sort_bam_by_name(unsorted_bam, sorted_unsorted)
            
            # Compare the name-sorted BAMs
            success = compare_bam_files(
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
        print("=== Stage 3 Complete: BAM Tags Comparison PASSED ===")
        print("All CB/UB/CR/UR/GX/GN tags match between sorted and unsorted runs")
        return 0
    else:
        print("=== Stage 3 FAILED: BAM Tags Comparison Failed ===")
        print("Differences found between sorted and unsorted runs")
        return 1


if __name__ == "__main__":
    sys.exit(main())
