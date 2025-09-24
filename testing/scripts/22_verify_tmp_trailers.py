#!/usr/bin/env python3
"""
Stage 2.2 - Verify temp file trailers
Purpose: validate the intermediate tmp file layout before pass 2 runs
"""

import os
import sys
import struct
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Verify STARsolo temp file trailer format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir",
        default=Path(__file__).parent / "../output/unsorted_twopass",
        type=Path,
        help="Directory containing the temp file"
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=10000,
        help="Maximum number of records to verify"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    return parser.parse_args()


def verify_tmp_file(tmp_path, max_records=10000, verbose=False):
    """
    Verify the tmp file format: [size][payload][trailer]
    
    Expected format per record:
    - uint32_t block_size_with_len: size of BAM record
    - block_size_with_len bytes: BAM record data
    - uint64_t trailer: (iReadAll << 32)
    """
    print(f"=== STARsolo Temp File Verification ===")
    print(f"File: {tmp_path}")
    print(f"Max records to check: {max_records}")
    print()
    
    if not tmp_path.exists():
        print(f"ERROR: Temp file not found: {tmp_path}")
        return False
    
    file_size = tmp_path.stat().st_size
    print(f"File size: {file_size:,} bytes ({file_size / (1024**2):.1f} MB)")
    
    records_checked = 0
    bytes_read = 0
    prev_iread = -1
    trailer_errors = 0
    size_errors = 0
    
    with open(tmp_path, 'rb') as f:
        while bytes_read < file_size and records_checked < max_records:
            # Read block size
            size_data = f.read(4)
            if len(size_data) < 4:
                if len(size_data) > 0:
                    print(f"ERROR: Incomplete size field at offset {bytes_read}")
                    return False
                break
            
            block_size = struct.unpack('<I', size_data)[0]  # Little-endian uint32
            bytes_read += 4
            
            if block_size == 0 or block_size > 1000000:  # Sanity check
                print(f"ERROR: Suspicious block size {block_size} at record {records_checked}")
                size_errors += 1
                if size_errors > 10:
                    print("Too many size errors, aborting")
                    return False
                continue
            
            # Read payload
            payload = f.read(block_size)
            if len(payload) < block_size:
                print(f"ERROR: Incomplete payload at record {records_checked} (expected {block_size}, got {len(payload)})")
                return False
            bytes_read += block_size
            
            # Read trailer
            trailer_data = f.read(8)
            if len(trailer_data) < 8:
                print(f"ERROR: Incomplete trailer at record {records_checked}")
                return False
            
            trailer = struct.unpack('<Q', trailer_data)[0]  # Little-endian uint64
            bytes_read += 8
            
            # Extract iread from trailer (should be trailer >> 32)
            iread = trailer >> 32
            
            # Verify iread is monotonic non-decreasing
            if prev_iread >= 0 and iread < prev_iread:
                print(f"ERROR: Non-monotonic iread at record {records_checked}: {iread} < {prev_iread}")
                trailer_errors += 1
            
            if verbose and records_checked < 10:
                print(f"Record {records_checked}: size={block_size}, iread={iread}, trailer=0x{trailer:016x}")
            
            prev_iread = iread
            records_checked += 1
            
            if records_checked % 1000 == 0:
                print(f"Verified {records_checked} records...")
    
    print()
    print(f"Verification complete:")
    print(f"  Records checked: {records_checked:,}")
    print(f"  Bytes processed: {bytes_read:,} / {file_size:,} ({100 * bytes_read / file_size:.1f}%)")
    print(f"  Size errors: {size_errors}")
    print(f"  Trailer errors: {trailer_errors}")
    print(f"  Final iread: {prev_iread}")
    
    if size_errors > 0 or trailer_errors > 0:
        print("FAILED: Errors detected in temp file format")
        return False
    
    if records_checked == 0:
        print("WARNING: No records found in temp file")
        return False
    
    print("PASSED: Temp file format appears correct")
    return True


def main():
    args = parse_args()
    
    print("=== STARsolo Unsorted CB/UB Testing - Stage 2.2: Verify Temp Trailers ===")
    print(f"Timestamp: {os.popen('date').read().strip()}")
    print()
    
    # Look for temp file
    tmp_file = args.input_dir / "Aligned.out.unsorted.solo.tmp"
    
    if not tmp_file.exists():
        print(f"Temp file not found: {tmp_file}")
        print("This is expected if KEEP_SOLO_TMP=1 was not set during the run")
        print("To preserve temp files, run 20_run_unsorted_two_pass.sh --keep-tmp")
        return 0
    
    success = verify_tmp_file(tmp_file, args.max_records, args.verbose)
    
    print()
    if success:
        print("=== Stage 2.2 Complete: Temp File Verification PASSED ===")
        return 0
    else:
        print("=== Stage 2.2 FAILED: Temp File Verification Failed ===")
        return 1


if __name__ == "__main__":
    sys.exit(main())
