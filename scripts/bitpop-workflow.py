#!/usr/bin/env python3
"""
Bit-Pop Multi-Index Workflow Tool

Automates the workflow for large genome datasets (>2GB):
1. Split large FASTA files into chunks
2. Build FM-index for each chunk in parallel
3. Map reads against all chunks in parallel
4. Merge SAM results into a single output

Usage:
    python bitpop-workflow.py split <input.fna> -o <output_dir> [--max-size 2000]
    python bitpop-workflow.py build <input_dir> -o <output_dir> [--threads 8]
    python bitpop-workflow.py map <input_dir> <reads.fastq> -o <output_dir> [--threads 8]
    python bitpop-workflow.py merge <input_dir> -o <output.sam>
    python bitpop-workflow.py full <genome.fna> <reads.fastq> -o <output_dir> [--threads 8]
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import time
from pathlib import Path
from typing import List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field


# ─── Constants ────────────────────────────────────────────────────────────────

DEFAULT_MAX_SIZE_MB = 2000  # 2GB per chunk
BITPOP_CMD = "bit-pop"
CHUNK_PREFIX = "chunk"
INDEX_SUFFIX = ".bitpop"
SAM_SUFFIX = ".sam"


# ─── Data Classes ─────────────────────────────────────────────────────────────

@dataclass
class ChunkInfo:
    """Information about a genome chunk."""
    name: str
    path: str
    size_mb: float


@dataclass
class WorkflowConfig:
    """Workflow configuration."""
    threads: int = 4
    max_size_mb: int = DEFAULT_MAX_SIZE_MB
    bitpop_cmd: str = BITPOP_CMD
    work_dir: str = "."
    output_dir: str = "."


# ─── FASTA Parsing ────────────────────────────────────────────────────────────

def parse_fasta_headers(filepath: str) -> List[Tuple[str, int, int]]:
    """Parse FASTA file and return list of (header, start_offset, end_offset).
    
    Returns entries with byte offsets for efficient splitting.
    """
    entries = []
    with open(filepath, 'rb') as f:
        content = f.read()
    
    lines = content.decode('ascii', errors='ignore').split('\n')
    current_header = None
    current_start = None
    current_pos = 0
    
    for line in lines:
        if line.startswith('>'):
            if current_header is not None:
                entries.append((current_header, current_start, current_pos))
            current_header = line[1:].strip().split()[0]  # Take first word after >
            current_start = current_pos
        current_pos += len(line) + 1  # +1 for newline
    
    if current_header is not None:
        entries.append((current_header, current_start, current_pos))
    
    return entries


def split_fasta_by_headers(filepath: str, output_dir: str, max_size_mb: int = DEFAULT_MAX_SIZE_MB) -> List[ChunkInfo]:
    """Split a FASTA file by headers, ensuring each chunk <= max_size_mb.
    
    Groups small headers together to avoid too many small files.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    entries = parse_fasta_headers(filepath)
    total_size = os.path.getsize(filepath)
    max_size_bytes = max_size_mb * 1024 * 1024
    
    chunks = []
    current_chunk = []
    current_size = 0
    
    for header, start, end in entries:
        entry_size = end - start
        if entry_size > max_size_bytes:
            # Single header exceeds max size - write it alone
            if current_chunk:
                chunk_name = f"{CHUNK_PREFIX}_{len(chunks)+1}"
                chunk_path = write_chunk(filepath, current_chunk, output_dir, chunk_name)
                chunks.append(ChunkInfo(chunk_name, chunk_path, current_size / (1024*1024)))
                current_chunk = []
                current_size = 0
            
            chunk_name = f"{CHUNK_PREFIX}_{len(chunks)+1}"
            chunk_path = write_chunk(filepath, [(header, start, end)], output_dir, chunk_name)
            chunks.append(ChunkInfo(chunk_name, chunk_path, entry_size / (1024*1024)))
        elif current_size + entry_size > max_size_bytes and current_chunk:
            # Current chunk would exceed limit - write it and start new one
            chunk_name = f"{CHUNK_PREFIX}_{len(chunks)+1}"
            chunk_path = write_chunk(filepath, current_chunk, output_dir, chunk_name)
            chunks.append(ChunkInfo(chunk_name, chunk_path, current_size / (1024*1024)))
            current_chunk = [(header, start, end)]
            current_size = entry_size
        else:
            current_chunk.append((header, start, end))
            current_size += entry_size
    
    # Write remaining chunk
    if current_chunk:
        chunk_name = f"{CHUNK_PREFIX}_{len(chunks)+1}"
        chunk_path = write_chunk(filepath, current_chunk, output_dir, chunk_name)
        chunks.append(ChunkInfo(chunk_name, chunk_path, current_size / (1024*1024)))
    
    return chunks


def write_chunk(filepath: str, entries: List[Tuple[str, int, int]], output_dir: str, chunk_name: str) -> str:
    """Write a chunk of FASTA entries to a new file."""
    chunk_path = os.path.join(output_dir, f"{chunk_name}.fna")
    
    with open(filepath, 'rb') as fin:
        with open(chunk_path, 'wb') as fout:
            for header, start, end in entries:
                fin.seek(start)
                # Read until next header or end of file
                pos = start
                while pos < end:
                    line_end = fin.read().find(b'\n', pos - start)
                    if line_end == -1:
                        fout.write(fin.read())
                        break
                    line_end += 1
                    fout.write(fin.read()[:line_end])
                    pos += line_end
    
    return chunk_path


def split_fasta_by_size(filepath: str, output_dir: str, max_size_mb: int = DEFAULT_MAX_SIZE_MB) -> List[ChunkInfo]:
    """Split a FASTA file into roughly equal-sized chunks by byte size.
    
    Splits at header boundaries to keep complete sequences.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    max_size_bytes = max_size_mb * 1024 * 1024
    total_size = os.path.getsize(filepath)
    
    with open(filepath, 'rb') as f:
        content = f.read()
    
    # Find all header positions
    headers = []
    for i, line in in_lines(content.decode('ascii', errors='ignore').split('\n')):
        if line.startswith('>'):
            headers.append(i)
    
    # Split at header boundaries
    chunks = []
    chunk_start = 0
    chunk_num = 1
    
    for i, header_idx in enumerate(headers):
        if i > 0:
            # Check if previous chunk is large enough
            prev_chunk_size = header_idx - chunk_start
            if prev_chunk_size >= max_size_bytes:
                chunk_path = write_chunk_by_lines(filepath, chunk_start, header_idx, output_dir, f"{CHUNK_PREFIX}_{chunk_num}")
                chunks.append(ChunkInfo(f"{CHUNK_PREFIX}_{chunk_num}", chunk_path, prev_chunk_size / (1024*1024)))
                chunk_num += 1
                chunk_start = header_idx
    
    # Write last chunk
    if chunk_start < len(content):
        chunk_path = write_chunk_by_lines(filepath, chunk_start, len(content), output_dir, f"{CHUNK_PREFIX}_{chunk_num}")
        chunks.append(ChunkInfo(f"{CHUNK_PREFIX}_{chunk_num}", chunk_path, (len(content) - chunk_start) / (1024*1024)))
    
    return chunks


def in_lines(text: str):
    """Helper to enumerate lines."""
    for i, line in enumerate(text.split('\n')):
        yield i, line


def write_chunk_by_lines(filepath: str, start_line: int, end_line: int, output_dir: str, chunk_name: str) -> str:
    """Write lines from start_line to end_line to a new file."""
    chunk_path = os.path.join(output_dir, f"{chunk_name}.fna")
    
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    
    with open(chunk_path, 'w') as fout:
        fout.writelines(lines[start_line:end_line])
    
    return chunk_path


# ─── Build Index ──────────────────────────────────────────────────────────────

def build_single_chunk(chunk_path: str, output_dir: str, threads: int, bitpop_cmd: str) -> Tuple[str, bool, str]:
    """Build FM-index for a single chunk."""
    chunk_name = Path(chunk_path).stem
    output_path = os.path.join(output_dir, f"{chunk_name}{INDEX_SUFFIX}")
    
    cmd = [
        bitpop_cmd, "build",
        "-f", chunk_path,
        "-o", output_path,
        "--threads", str(threads),
        "--auto-k"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode == 0:
            return (chunk_name, True, "Success")
        else:
            return (chunk_name, False, result.stderr[:500])
    except subprocess.TimeoutExpired:
        return (chunk_name, False, "Timeout after 1 hour")
    except Exception as e:
        return (chunk_name, False, str(e)[:500])


def build_chunks(chunk_dir: str, output_dir: str, threads: int = 4, bitpop_cmd: str = BITPOP_CMD) -> List[ChunkInfo]:
    """Build FM-index for all chunks in parallel."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all chunk files
    chunk_files = sorted([
        os.path.join(chunk_dir, f) for f in os.listdir(chunk_dir)
        if f.endswith('.fna')
    ])
    
    if not chunk_files:
        print(f"Error: No .fna files found in {chunk_dir}")
        sys.exit(1)
    
    print(f"Building {len(chunk_files)} indexes with {threads} threads...")
    start_time = time.time()
    
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(build_single_chunk, chunk_path, output_dir, threads, bitpop_cmd): chunk_path
            for chunk_path in chunk_files
        }
        
        for future in as_completed(futures):
            chunk_path = futures[future]
            chunk_name = Path(chunk_path).stem
            print(f"  Building {chunk_name}...", end=" ")
            
            try:
                name, success, msg = future.result()
                if success:
                    print(f"OK ({os.path.getsize(os.path.join(output_dir, f'{name}{INDEX_SUFFIX}')) / (1024*1024):.1f}MB)")
                    results.append(ChunkInfo(name, os.path.join(output_dir, f'{name}{INDEX_SUFFIX}'), 
                                            os.path.getsize(os.path.join(output_dir, f'{name}{INDEX_SUFFIX}')) / (1024*1024)))
                else:
                    print(f"FAILED: {msg}")
            except Exception as e:
                print(f"ERROR: {e}")
    
    elapsed = time.time() - start_time
    print(f"\nBuild complete: {len(results)}/{len(chunk_files)} indexes in {elapsed:.1f}s")
    
    return results


# ─── Map Reads ────────────────────────────────────────────────────────────────

def map_single_index(index_path: str, reads_path: str, output_dir: str, threads: int, bitpop_cmd: str) -> Tuple[str, bool, str]:
    """Map reads against a single index."""
    index_name = Path(index_path).stem
    output_path = os.path.join(output_dir, f"{index_name}{SAM_SUFFIX}")
    
    cmd = [
        bitpop_cmd, "run",
        index_path,
        reads_path,
        "-o", output_path,
        "--threads", str(threads),
        "--min-score", "0.5"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode == 0:
            return (index_name, True, "Success")
        else:
            return (index_name, False, result.stderr[:500])
    except subprocess.TimeoutExpired:
        return (index_name, False, "Timeout after 2 hours")
    except Exception as e:
        return (index_name, False, str(e)[:500])


def map_chunks(index_dir: str, reads_path: str, output_dir: str, threads: int = 4, bitpop_cmd: str = BITPOP_CMD) -> List[ChunkInfo]:
    """Map reads against all indexes in parallel."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all index files
    index_files = sorted([
        os.path.join(index_dir, f) for f in os.listdir(index_dir)
        if f.endswith(INDEX_SUFFIX)
    ])
    
    if not index_files:
        print(f"Error: No .bitpop files found in {index_dir}")
        sys.exit(1)
    
    print(f"Mapping reads against {len(index_files)} indexes with {threads} threads...")
    start_time = time.time()
    
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(map_single_index, index_path, reads_path, output_dir, threads, bitpop_cmd): index_path
            for index_path in index_files
        }
        
        for future in as_completed(futures):
            index_path = futures[future]
            index_name = Path(index_path).stem
            print(f"  Mapping {index_name}...", end=" ")
            
            try:
                name, success, msg = future.result()
                if success:
                    sam_path = os.path.join(output_dir, f'{name}{SAM_SUFFIX}')
                    size_mb = os.path.getsize(sam_path) / (1024*1024) if os.path.exists(sam_path) else 0
                    print(f"OK ({size_mb:.1f}MB)")
                    results.append(ChunkInfo(name, sam_path, size_mb))
                else:
                    print(f"FAILED: {msg}")
            except Exception as e:
                print(f"ERROR: {e}")
    
    elapsed = time.time() - start_time
    print(f"\nMapping complete: {len(results)}/{len(index_files)} results in {elapsed:.1f}s")
    
    return results


# ─── Merge SAM ────────────────────────────────────────────────────────────────

def merge_sam_files(sam_files: List[str], output_path: str) -> int:
    """Merge multiple SAM files into one, deduplicating by read name.
    
    Returns the number of unique reads in the merged output.
    """
    if not sam_files:
        print("Error: No SAM files to merge")
        return 0
    
    # Read header from first file
    header_lines = []
    header_found = False
    
    # Collect all data lines, indexed by read name
    read_data = {}  # read_name -> (header_line, data_line)
    total_lines = 0
    duplicated = 0
    
    for sam_file in sam_files:
        if not os.path.exists(sam_file):
            print(f"  Warning: {sam_file} not found, skipping")
            continue
        
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    # Header line - collect unique headers
                    if not header_found:
                        header_lines.append(line)
                        if not line.startswith('@HD') and not header_found:
                            header_found = True
                else:
                    total_lines += 1
                    # Extract read name (first field)
                    read_name = line.split('\t')[0]
                    if read_name in read_data:
                        duplicated += 1
                    else:
                        read_data[read_name] = line
    
    # Write merged output
    with open(output_path, 'w') as fout:
        # Write header
        for line in header_lines:
            fout.write(line)
        
        # Write data lines (deduplicated)
        for read_name, data_line in read_data.items():
            fout.write(data_line)
    
    unique_reads = len(read_data)
    
    print(f"  Merged {len(sam_files)} files")
    print(f"  Total lines: {total_lines}")
    print(f"  Duplicates removed: {duplicated}")
    print(f"  Unique reads: {unique_reads}")
    print(f"  Output: {output_path} ({os.path.getsize(output_path) / (1024*1024):.1f}MB)")
    
    return unique_reads


# ─── Full Workflow ────────────────────────────────────────────────────────────

def run_full_workflow(genome_path: str, reads_path: str, output_dir: str, 
                      threads: int = 4, max_size_mb: int = DEFAULT_MAX_SIZE_MB,
                      bitpop_cmd: str = BITPOP_CMD, cleanup: bool = False):
    """Run the complete split-build-map-merge workflow.
    
    1. Split genome into chunks
    2. Build index for each chunk
    3. Map reads against each index
    4. Merge SAM results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    chunk_dir = os.path.join(output_dir, "chunks")
    index_dir = os.path.join(output_dir, "indexes")
    map_dir = os.path.join(output_dir, "mapped")
    
    total_start = time.time()
    
    # Step 1: Split
    print("=" * 60)
    print("Step 1/4: Splitting genome...")
    print("=" * 60)
    chunks = split_fasta_by_headers(genome_path, chunk_dir, max_size_mb)
    print(f"Split into {len(chunks)} chunks ({sum(c.size_mb for c in chunks):.1f}MB total)")
    
    # Step 2: Build
    print("=" * 60)
    print("Step 2/4: Building indexes...")
    print("=" * 60)
    indexes = build_chunks(chunk_dir, index_dir, threads, bitpop_cmd)
    
    # Step 3: Map
    print("=" * 60)
    print("Step 3/4: Mapping reads...")
    print("=" * 60)
    mapped = map_chunks(index_dir, reads_path, map_dir, threads, bitpop_cmd)
    
    # Step 4: Merge
    print("=" * 60)
    print("Step 4/4: Merging results...")
    print("=" * 60)
    sam_files = [m.path for m in mapped if os.path.exists(m.path)]
    output_sam = os.path.join(output_dir, "final.sam")
    unique_reads = merge_sam_files(sam_files, output_sam)
    
    # Cleanup
    if cleanup:
        print("\nCleaning up intermediate files...")
        shutil.rmtree(chunk_dir, ignore_errors=True)
        shutil.rmtree(index_dir, ignore_errors=True)
        shutil.rmtree(map_dir, ignore_errors=True)
    
    total_elapsed = time.time() - total_start
    print("\n" + "=" * 60)
    print(f"Workflow complete!")
    print(f"  Unique reads mapped: {unique_reads}")
    print(f"  Output: {output_sam}")
    print(f"  Total time: {total_elapsed:.1f}s")
    print("=" * 60)


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Bit-Pop Multi-Index Workflow Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split genome into chunks
  python bitpop-workflow.py split genome.fna -o chunks/ --max-size 2000

  # Build indexes for all chunks
  python bitpop-workflow.py build chunks/ -o indexes/ --threads 8

  # Map reads against all indexes
  python bitpop-workflow.py map indexes/ reads.fastq -o mapped/ --threads 8

  # Merge SAM results
  python bitpop-workflow.py merge mapped/ -o final.sam

  # Full workflow (all steps)
  python bitpop-workflow.py full genome.fna reads.fastq -o output/ --threads 8
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Split command
    split_parser = subparsers.add_parser('split', help='Split FASTA into chunks')
    split_parser.add_argument('genome', help='Input FASTA file')
    split_parser.add_argument('-o', '--output', required=True, help='Output directory')
    split_parser.add_argument('--max-size', type=int, default=DEFAULT_MAX_SIZE_MB, help=f'Max chunk size in MB (default: {DEFAULT_MAX_SIZE_MB})')
    
    # Build command
    build_parser = subparsers.add_parser('build', help='Build indexes for all chunks')
    build_parser.add_argument('input_dir', help='Directory with chunk .fna files')
    build_parser.add_argument('-o', '--output', required=True, help='Output directory for indexes')
    build_parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    # Map command
    map_parser = subparsers.add_parser('map', help='Map reads against all indexes')
    map_parser.add_argument('index_dir', help='Directory with .bitpop index files')
    map_parser.add_argument('reads', help='Input FASTQ reads file')
    map_parser.add_argument('-o', '--output', required=True, help='Output directory for SAM files')
    map_parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    # Merge command
    merge_parser = subparsers.add_parser('merge', help='Merge SAM results')
    merge_parser.add_argument('input_dir', help='Directory with .sam files')
    merge_parser.add_argument('-o', '--output', required=True, help='Output SAM file')
    
    # Full workflow command
    full_parser = subparsers.add_parser('full', help='Run complete workflow (split + build + map + merge)')
    full_parser.add_argument('genome', help='Input FASTA genome file')
    full_parser.add_argument('reads', help='Input FASTQ reads file')
    full_parser.add_argument('-o', '--output', required=True, help='Output directory')
    full_parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    full_parser.add_argument('--max-size', type=int, default=DEFAULT_MAX_SIZE_MB, help=f'Max chunk size in MB (default: {DEFAULT_MAX_SIZE_MB})')
    full_parser.add_argument('--no-cleanup', action='store_true', help='Keep intermediate files')
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    if args.command == 'split':
        chunks = split_fasta_by_headers(args.genome, args.output, args.max_size)
        print(f"\nSplit into {len(chunks)} chunks:")
        for chunk in chunks:
            print(f"  {chunk.name}: {chunk.size_mb:.1f}MB")
    
    elif args.command == 'build':
        build_chunks(args.input_dir, args.output, args.threads)
    
    elif args.command == 'map':
        map_chunks(args.index_dir, args.reads, args.output, args.threads)
    
    elif args.command == 'merge':
        sam_files = sorted([
            os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir)
            if f.endswith(SAM_SUFFIX)
        ])
        merge_sam_files(sam_files, args.output)
    
    elif args.command == 'full':
        cleanup = not args.no_cleanup
        run_full_workflow(
            genome_path=args.genome,
            reads_path=args.reads,
            output_dir=args.output,
            threads=args.threads,
            max_size_mb=args.max_size,
            cleanup=cleanup
        )


if __name__ == "__main__":
    main()
