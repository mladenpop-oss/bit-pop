#!/usr/bin/env python3
"""
consensus_base.py - Simple consensus using standalone `bit-pop map` for each index.
Runs `bit-pop map` for each index, then combines SAM files with consensus voting.
"""
import subprocess
import sys
import os
import tempfile
import argparse
from collections import defaultdict

def parse_fastq(fastq_path):
    """Parse FASTQ and return dict: read_name -> sequence"""
    sequences = {}
    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # quality
            name = header[1:]  # remove @
            sequences[name] = seq
    return sequences

def parse_sam(sam_path):
    """Parse SAM file and return dict: read_name -> list of (genome_name, genome_id, score, cigar, pos, is_reverse, rarity)"""
    results = defaultdict(list)
    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            read_name = fields[0]
            flag = int(fields[1])
            genome_name = fields[2]
            if genome_name == '*':
                continue
            pos = int(fields[3]) - 1  # 1-based to 0-based
            cigar = fields[5]
            
            score = 0.0
            rarity = 0.0
            genome_id = 0
            is_reverse = bool(flag & 0x10)
            
            for tag in fields[11:]:
                if tag.startswith('AS:f:'):
                    score = float(tag[5:])
                elif tag.startswith('RK:f:'):
                    rarity = float(tag[5:])
                elif tag.startswith('KK:i:'):
                    genome_id = int(tag[5:])
            
            results[read_name].append({
                'genome_name': genome_name,
                'genome_id': genome_id,
                'score': score,
                'cigar': cigar,
                'pos': pos,
                'is_reverse': is_reverse,
                'rarity': rarity,
            })
    return results

def build_name_to_id(sam_path):
    """Build genome_name -> genome_id mapping from SAM @SQ lines."""
    name_to_id = {}
    gid = 0
    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@SQ'):
                fields = line.strip().split('\t')
                name = None
                for field in fields:
                    if field.startswith('SN:'):
                        name = field[3:]
                if name:
                    name_to_id[name] = gid
                    gid += 1
    return name_to_id

def vote_weighted(results):
    """Weighted score consensus: sum scores per genome, return sorted by total score."""
    genome_scores = defaultdict(lambda: {'score': 0.0, 'count': 0, 'best': None})
    for r in results:
        g = r['genome_name']
        genome_scores[g]['score'] += r['score']
        genome_scores[g]['count'] += 1
        if genome_scores[g]['best'] is None or r['score'] > genome_scores[g]['best']['score']:
            genome_scores[g]['best'] = r
    sorted_genomes = sorted(genome_scores.items(), key=lambda x: x[1]['score'], reverse=True)
    return sorted_genomes

def vote_best_score(results):
    """Best score: return the single best result across all k-values."""
    best = max(results, key=lambda r: r['score'])
    return [(best['genome_name'], {'score': best['score'], 'count': 1, 'best': best})]

def vote_majority(results):
    """Majority vote: count votes per genome, return sorted by vote count."""
    genome_votes = defaultdict(lambda: {'count': 0, 'best': None})
    for r in results:
        g = r['genome_name']
        genome_votes[g]['count'] += 1
        if genome_votes[g]['best'] is None or r['score'] > genome_votes[g]['best']['score']:
            genome_votes[g]['best'] = r
    sorted_genomes = sorted(genome_votes.items(), key=lambda x: x[1]['count'], reverse=True)
    return sorted_genomes

def main():
    parser = argparse.ArgumentParser(description='Multi-k consensus using bit-pop map')
    parser.add_argument('-i', '--indexes', nargs='+', required=True,
                        help='Index files (e.g., cami_k10.bitpop cami_k13.bitpop)')
    parser.add_argument('-r', '--reads', required=True, help='Reads file (FASTQ/FASTA)')
    parser.add_argument('-o', '--output', required=True, help='Output SAM file')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Threads per map')
    parser.add_argument('--top-n', type=int, default=1, help='Top-N anchors')
    parser.add_argument('--strategy', choices=['weighted_score', 'best_score', 'majority'],
                        default='weighted_score', help='Consensus strategy')
    parser.add_argument('--min-score', type=float, default=0.0, help='Min score filter')
    parser.add_argument('--min-k', type=int, default=1, help='Min k-values that must map')
    parser.add_argument('--context-window', type=int, default=50, help='Context window')
    parser.add_argument('--bit-pop', default=None, help='Path to bit-pop executable (auto: target/release/bit-pop.exe)')
    args = parser.parse_args()
    
    # Auto-detect bit-pop path
    if args.bit_pop is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)  # scripts/ -> repo root
        bit_pop_exe = os.path.join(repo_root, 'target', 'release', 'bit-pop.exe')
        if os.path.exists(bit_pop_exe):
            args.bit_pop = bit_pop_exe
        else:
            print(f"Error: bit-pop.exe not found at {bit_pop_exe}")
            sys.exit(1)

    # Phase 1: Run bit-pop map for each index
    temp_files = []
    print("Phase 1: Mapping each index with bit-pop map")
    print("=" * 50)
    
    for index_path in args.indexes:
        temp_sam = tempfile.NamedTemporaryFile(suffix='.sam', delete=False, mode='w')
        temp_sam.close()
        temp_files.append(temp_sam.name)
        
        cmd = [
            args.bit_pop, 'map',
            '-i', index_path,
            '-r', args.reads,
            '-o', temp_sam.name,
            '-a', 'xor',
            '--top-n', str(args.top_n),
            '-t', str(args.threads),
        ]
        print(f"\n  Mapping: {os.path.basename(index_path)}")
        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            print(f"Error mapping {index_path}")
            sys.exit(1)

    # Phase 2: Parse reads and combine
    print("\n" + "=" * 50)
    print("Phase 2: Combining results")
    
    # Parse reads to get sequences
    print("  Parsing reads...")
    read_seqs = parse_fastq(args.reads)
    print(f"  Loaded {len(read_seqs)} reads")
    
    # Build genome name->id from first temp SAM
    name_to_id = build_name_to_id(temp_files[0])
    
    # Parse all temp SAMs
    all_results = {}  # read_name -> list of results from all k
    for i, temp_sam in enumerate(temp_files):
        print(f"  Parsing {os.path.basename(args.indexes[i])}...")
        parsed = parse_sam(temp_sam)
        for read_name, results in parsed.items():
            if read_name not in all_results:
                all_results[read_name] = []
            all_results[read_name].extend(results)
    
    # Clean up temp files
    for f in temp_files:
        os.unlink(f)
    
    # Phase 3: Write consensus SAM
    print("\nPhase 3: Writing consensus SAM")
    
    with open(args.output, 'w') as out:
        out.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, gid in name_to_id.items():
            out.write(f"@SQ\tSN:{name}\tLN:0\n")
        
        mapped = 0
        total = len(read_seqs)
        
        for read_name, seq in read_seqs.items():
            results = all_results.get(read_name, [])
            if len(results) < args.min_k:
                continue
            
            # Filter by min_score
            filtered = [r for r in results if r['score'] >= args.min_score]
            if not filtered:
                continue
            
            # Consensus voting
            if args.strategy == 'weighted_score':
                voted = vote_weighted(filtered)
            elif args.strategy == 'best_score':
                voted = vote_best_score(filtered)
            else:
                voted = vote_majority(filtered)
            
            # Write best result
            if voted:
                best_name, best_data = voted[0]
                best_r = best_data['best']
                
                flag = 0x10 if best_r['is_reverse'] else 0
                mapq = 0
                
                # Tags
                as_tag = f"\tAS:f:{best_r['score']:.4f}"
                xs_tag = f"\tXS:Z:{best_r['rarity']:.4f}"
                kk_tag = f"\tKK:i:{best_data['count']}"
                
                out.write(f"{read_name}\t{flag}\t{best_name}\t{best_r['pos']+1}\t{mapq}\t{best_r['cigar']}\t*\t0\t0\t{seq}\t*{as_tag}{xs_tag}{kk_tag}\n")
                mapped += 1
        
        print(f"  Mapped: {mapped} / {total} reads")
    
    print(f"\nDone! Output: {args.output}")

if __name__ == '__main__':
    main()
