#!/usr/bin/env python3
"""Split paired-end FASTQ into R1 and R2 files."""

import sys


def split_fastq(input_path, r1_path, r2_path):
    """Split a paired-end FASTQ file into R1 and R2."""
    print(f"Reading {input_path}...")
    
    r1_count = 0
    r2_count = 0
    
    with open(input_path, 'r') as infile, open(r1_path, 'w') as r1, open(r2_path, 'w') as r2:
        while True:
            # Read 4 lines per read (FASTQ format)
            header1 = infile.readline()
            if not header1:
                break
            seq1 = infile.readline()
            plus1 = infile.readline()
            qual1 = infile.readline()
            
            header2 = infile.readline()
            seq2 = infile.readline()
            plus2 = infile.readline()
            qual2 = infile.readline()
            
            # Write R1
            r1.write(header1)
            r1.write(seq1)
            r1.write(plus1)
            r1.write(qual1)
            r1_count += 1
            
            # Write R2
            r2.write(header2)
            r2.write(seq2)
            r2.write(plus2)
            r2.write(qual2)
            r2_count += 1
    
    print(f"Split complete!")
    print(f"  R1: {r1_count} reads -> {r1_path}")
    print(f"  R2: {r2_count} reads -> {r2_path}")
    print(f"  Total pairs: {r1_count}")


if __name__ == '__main__':
    # Default paths
    input_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_reads.fq'
    r1_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_reads_R1.fq'
    r2_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_reads_R2.fq'
    
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
    if len(sys.argv) > 2:
        r1_path = sys.argv[2]
    if len(sys.argv) > 3:
        r2_path = sys.argv[3]
    
    split_fastq(input_path, r1_path, r2_path)
