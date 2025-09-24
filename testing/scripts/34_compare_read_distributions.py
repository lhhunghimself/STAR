#!/usr/bin/env python3
"""
Stage 3.4 - Compare read distributions
Purpose: deeper sanity on barcode/cell counts
"""

import os
import sys
import argparse
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare read distributions and barcode counts between runs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--baseline-dir",
        default=Path(__file__).parent / "../output/baseline_sorted",
        type=Path,
        help="Directory containing sorted baseline outputs"
    )
    parser.add_argument(
        "--unsorted-dir",
        default=Path(__file__).parent / "../output/unsorted_twopass",
        type=Path,
        help="Directory containing unsorted outputs"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--sample-barcodes",
        type=int,
        default=1000,
        help="Sample N barcodes for detailed comparison (0 = all)"
    )
    return parser.parse_args()


def load_matrix_market(mtx_file, features_file, barcodes_file):
    """Load a Matrix Market format sparse matrix with features and barcodes"""
    try:
        from scipy.io import mmread
        
        print(f"Loading matrix: {mtx_file}")
        matrix = mmread(mtx_file)
        
        # Load features and barcodes
        features = []
        if features_file.exists():
            with open(features_file, 'r') as f:
                for line in f:
                    features.append(line.strip().split('\t'))
        
        barcodes = []
        if barcodes_file.exists():
            with open(barcodes_file, 'r') as f:
                for line in f:
                    barcodes.append(line.strip())
        
        print(f"  Matrix shape: {matrix.shape}")
        print(f"  Features: {len(features)}")
        print(f"  Barcodes: {len(barcodes)}")
        print(f"  Non-zero entries: {matrix.nnz}")
        
        return matrix, features, barcodes
        
    except ImportError:
        print("Warning: scipy not available, skipping matrix analysis")
        return None, None, None


def compare_barcode_files(baseline_file, unsorted_file, verbose=False):
    """Compare barcode files line by line"""
    print(f"Comparing barcode files:")
    print(f"  Baseline: {baseline_file}")
    print(f"  Unsorted: {unsorted_file}")
    
    if not baseline_file.exists():
        print(f"  ERROR: Baseline file missing")
        return False
    
    if not unsorted_file.exists():
        print(f"  ERROR: Unsorted file missing")
        return False
    
    # Read barcodes
    with open(baseline_file, 'r') as f:
        baseline_barcodes = [line.strip() for line in f]
    
    with open(unsorted_file, 'r') as f:
        unsorted_barcodes = [line.strip() for line in f]
    
    print(f"  Baseline count: {len(baseline_barcodes)}")
    print(f"  Unsorted count: {len(unsorted_barcodes)}")
    
    if len(baseline_barcodes) != len(unsorted_barcodes):
        print(f"  FAIL: Different barcode counts")
        return False
    
    # Compare barcode lists
    differences = 0
    for i, (b1, b2) in enumerate(zip(baseline_barcodes, unsorted_barcodes)):
        if b1 != b2:
            differences += 1
            if differences <= 10:  # Show first 10 differences
                print(f"    Diff at line {i+1}: '{b1}' vs '{b2}'")
    
    if differences > 10:
        print(f"    ... ({differences - 10} more differences)")
    
    if differences == 0:
        print(f"  PASS: All barcodes match")
        return True
    else:
        print(f"  FAIL: {differences} barcode differences")
        return False


def compare_feature_files(baseline_file, unsorted_file, verbose=False):
    """Compare feature files line by line"""
    print(f"Comparing feature files:")
    print(f"  Baseline: {baseline_file}")
    print(f"  Unsorted: {unsorted_file}")
    
    if not baseline_file.exists():
        print(f"  ERROR: Baseline file missing")
        return False
    
    if not unsorted_file.exists():
        print(f"  ERROR: Unsorted file missing")
        return False
    
    # Read features
    baseline_features = []
    with open(baseline_file, 'r') as f:
        for line in f:
            baseline_features.append(line.strip().split('\t'))
    
    unsorted_features = []
    with open(unsorted_file, 'r') as f:
        for line in f:
            unsorted_features.append(line.strip().split('\t'))
    
    print(f"  Baseline count: {len(baseline_features)}")
    print(f"  Unsorted count: {len(unsorted_features)}")
    
    if len(baseline_features) != len(unsorted_features):
        print(f"  FAIL: Different feature counts")
        return False
    
    # Compare feature lists
    differences = 0
    for i, (f1, f2) in enumerate(zip(baseline_features, unsorted_features)):
        if f1 != f2:
            differences += 1
            if differences <= 10:  # Show first 10 differences
                print(f"    Diff at line {i+1}: {f1} vs {f2}")
    
    if differences > 10:
        print(f"    ... ({differences - 10} more differences)")
    
    if differences == 0:
        print(f"  PASS: All features match")
        return True
    else:
        print(f"  FAIL: {differences} feature differences")
        return False


def compare_matrices(baseline_dir, unsorted_dir, sample_barcodes=1000, verbose=False):
    """Compare matrix files between runs"""
    print("Comparing matrix files...")
    
    # Define matrix sets to compare
    matrix_sets = [
        {
            'name': 'Gene (unfiltered)',
            'matrix': 'Solo.out/Gene/matrix.mtx',
            'features': 'Solo.out/Gene/Features.tsv',
            'barcodes': 'Solo.out/Gene/Barcodes.tsv'
        },
        {
            'name': 'Gene (filtered)',
            'matrix': 'Solo.out/Gene/Filtered/matrix.mtx',
            'features': 'Solo.out/Gene/Filtered/features.tsv',
            'barcodes': 'Solo.out/Gene/Filtered/barcodes.tsv'
        }
    ]
    
    all_passed = True
    
    for matrix_set in matrix_sets:
        print(f"\n--- {matrix_set['name']} ---")
        
        # File paths
        baseline_matrix = baseline_dir / matrix_set['matrix']
        unsorted_matrix = unsorted_dir / matrix_set['matrix']
        baseline_features = baseline_dir / matrix_set['features']
        unsorted_features = unsorted_dir / matrix_set['features']
        baseline_barcodes = baseline_dir / matrix_set['barcodes']
        unsorted_barcodes = unsorted_dir / matrix_set['barcodes']
        
        # Check if files exist
        if not baseline_matrix.exists():
            print(f"Skipping {matrix_set['name']} (baseline matrix missing)")
            continue
        
        if not unsorted_matrix.exists():
            print(f"FAIL: {matrix_set['name']} (unsorted matrix missing)")
            all_passed = False
            continue
        
        # Compare barcodes
        if not compare_barcode_files(baseline_barcodes, unsorted_barcodes, verbose):
            all_passed = False
        
        # Compare features
        if not compare_feature_files(baseline_features, unsorted_features, verbose):
            all_passed = False
        
        # Compare matrix dimensions and checksums
        print(f"Comparing matrix files:")
        baseline_size = baseline_matrix.stat().st_size
        unsorted_size = unsorted_matrix.stat().st_size
        
        print(f"  Baseline matrix size: {baseline_size:,} bytes")
        print(f"  Unsorted matrix size: {unsorted_size:,} bytes")
        
        if baseline_size != unsorted_size:
            print(f"  FAIL: Matrix sizes differ")
            all_passed = False
            continue
        
        # Compare matrix content
        import subprocess
        baseline_md5 = subprocess.check_output(['md5sum', str(baseline_matrix)]).decode().split()[0]
        unsorted_md5 = subprocess.check_output(['md5sum', str(unsorted_matrix)]).decode().split()[0]
        
        if baseline_md5 != unsorted_md5:
            print(f"  FAIL: Matrix checksums differ")
            print(f"    Baseline: {baseline_md5}")
            print(f"    Unsorted: {unsorted_md5}")
            all_passed = False
        else:
            print(f"  PASS: Matrix files identical (MD5: {baseline_md5[:8]}...)")
        
        # Load and analyze matrices if scipy is available
        baseline_matrix_data, baseline_features_data, baseline_barcodes_data = load_matrix_market(
            baseline_matrix, baseline_features, baseline_barcodes
        )
        
        if baseline_matrix_data is not None:
            unsorted_matrix_data, unsorted_features_data, unsorted_barcodes_data = load_matrix_market(
                unsorted_matrix, unsorted_features, unsorted_barcodes
            )
            
            if unsorted_matrix_data is not None:
                # Compare matrix statistics
                print(f"Matrix statistics comparison:")
                print(f"  Shape match: {baseline_matrix_data.shape == unsorted_matrix_data.shape}")
                print(f"  NNZ match: {baseline_matrix_data.nnz == unsorted_matrix_data.nnz}")
                
                if baseline_matrix_data.shape == unsorted_matrix_data.shape:
                    # Sample some barcodes for detailed comparison
                    if sample_barcodes > 0 and len(baseline_barcodes_data) > sample_barcodes:
                        sample_indices = np.random.choice(
                            len(baseline_barcodes_data), 
                            sample_barcodes, 
                            replace=False
                        )
                        print(f"  Sampling {sample_barcodes} barcodes for detailed comparison...")
                        
                        # Compare UMI counts for sampled barcodes
                        baseline_sample = baseline_matrix_data[:, sample_indices].sum(axis=0).A1
                        unsorted_sample = unsorted_matrix_data[:, sample_indices].sum(axis=0).A1
                        
                        if np.array_equal(baseline_sample, unsorted_sample):
                            print(f"  PASS: Sampled barcode UMI counts match")
                        else:
                            print(f"  FAIL: Sampled barcode UMI counts differ")
                            diff_count = np.sum(baseline_sample != unsorted_sample)
                            print(f"    Differing barcodes: {diff_count}/{sample_barcodes}")
                            all_passed = False
    
    return all_passed


def compare_summary_stats(baseline_dir, unsorted_dir, verbose=False):
    """Compare summary statistics files"""
    print("\nComparing summary statistics...")
    
    summary_files = [
        "Solo.out/Gene/Summary.csv",
        "Solo.out/GeneFull/Summary.csv"
    ]
    
    all_passed = True
    
    for summary_file in summary_files:
        baseline_file = baseline_dir / summary_file
        unsorted_file = unsorted_dir / summary_file
        
        if not baseline_file.exists():
            if verbose:
                print(f"Skipping {summary_file} (not present in baseline)")
            continue
        
        if not unsorted_file.exists():
            print(f"FAIL: {summary_file} (missing in unsorted)")
            all_passed = False
            continue
        
        print(f"Comparing {summary_file}:")
        
        try:
            # Read CSV files
            baseline_df = pd.read_csv(baseline_file, index_col=0)
            unsorted_df = pd.read_csv(unsorted_file, index_col=0)
            
            # Compare shapes
            if baseline_df.shape != unsorted_df.shape:
                print(f"  FAIL: Different shapes ({baseline_df.shape} vs {unsorted_df.shape})")
                all_passed = False
                continue
            
            # Compare values
            if baseline_df.equals(unsorted_df):
                print(f"  PASS: Summary statistics match")
            else:
                print(f"  FAIL: Summary statistics differ")
                
                # Show differences
                diff_mask = baseline_df != unsorted_df
                if diff_mask.any().any():
                    print(f"    Differing entries:")
                    for idx in baseline_df.index:
                        for col in baseline_df.columns:
                            if pd.notna(diff_mask.loc[idx, col]) and diff_mask.loc[idx, col]:
                                baseline_val = baseline_df.loc[idx, col]
                                unsorted_val = unsorted_df.loc[idx, col]
                                print(f"      {idx}, {col}: {baseline_val} vs {unsorted_val}")
                
                all_passed = False
        
        except Exception as e:
            print(f"  ERROR: Failed to compare {summary_file}: {e}")
            all_passed = False
    
    return all_passed


def main():
    args = parse_args()
    
    print("=== STARsolo Unsorted CB/UB Testing - Stage 3.4: Compare Read Distributions ===")
    print(f"Timestamp: {os.popen('date').read().strip()}")
    print(f"Baseline directory: {args.baseline_dir}")
    print(f"Unsorted directory: {args.unsorted_dir}")
    print()
    
    # Check input directories
    if not args.baseline_dir.exists():
        print(f"ERROR: Baseline directory not found: {args.baseline_dir}")
        return 1
    
    if not args.unsorted_dir.exists():
        print(f"ERROR: Unsorted directory not found: {args.unsorted_dir}")
        return 1
    
    all_passed = True
    
    # Compare matrices and barcode distributions
    if not compare_matrices(
        args.baseline_dir, 
        args.unsorted_dir, 
        sample_barcodes=args.sample_barcodes,
        verbose=args.verbose
    ):
        all_passed = False
    
    # Compare summary statistics
    if not compare_summary_stats(args.baseline_dir, args.unsorted_dir, args.verbose):
        all_passed = False
    
    print()
    if all_passed:
        print("=== Stage 3.4 Complete: Read Distributions Comparison PASSED ===")
        print("All barcode counts and distributions match between runs")
        return 0
    else:
        print("=== Stage 3.4 FAILED: Read Distributions Comparison Failed ===")
        print("Differences found in barcode counts or distributions")
        return 1


if __name__ == "__main__":
    sys.exit(main())
