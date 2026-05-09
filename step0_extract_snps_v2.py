"""
Step 0: Direct SNP Extraction v2 — seed-and-extend approach.

Uses minimizer-based seeding for robust alignment of evo vs parent genomes.
Handles assembly artifacts by finding collinear seed chains.

Usage:
    python step0_extract_snps_v2.py [options]
"""

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path


# ─── FASTA Reader ────────────────────────────────────────────────────────────

def read_fasta(filepath):
    """Read FASTA file, return dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


# ─── Minimizer-Based Seeding ────────────────────────────────────────────────

def compute_minimizers(sequence, k=21, w=10):
    """
    Compute minimizers of a sequence.
    
    minimizers: for each window of w consecutive k-mers, select the lexicographically
    smallest k-mer as the minimizer. Returns list of (position, minimizer_hash).
    
    This gives sparse but representative seeding.
    """
    if len(sequence) < k:
        return []
    
    minimizers = []
    
    # Use simple hash for speed
    def kmer_hash(kmer):
        h = 0
        for c in kmer:
            h = (h * 4 + ord(c)) % (2**30)
        return h
    
    # Sliding window of minimizers
    window = []  # (hash, position)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        h = kmer_hash(kmer)
        window.append((h, i))
        
        # Remove entries outside window
        while window and window[0][1] < i - w + 1:
            window.pop(0)
        
        # Add minimizer
        if len(window) == w or i >= len(sequence) - k:
            minimizers.append((window[0][1], window[0][0]))
    
    return minimizers


def find_seed_matches(query, reference, k=21, w=10):
    """
    Find seed matches using minimizers.
    Returns list of (query_pos, ref_pos) pairs.
    """
    # Compute minimizers for both sequences
    query_minims = compute_minimizers(query, k, w)
    ref_minims = compute_minimizers(reference, k, w)
    
    # Build index of reference minimizers
    minimizer_index = defaultdict(list)
    for pos, h in ref_minims:
        minimizer_index[h].append(pos)
    
    # Find matches
    matches = []
    for q_pos, q_h in query_minims:
        if q_h in minimizer_index:
            for r_pos in minimizer_index[q_h]:
                matches.append((q_pos, r_pos))
    
    return matches


# ─── Seed Chaining ───────────────────────────────────────────────────────────

def chain_seeds(matches, max_gap=10000):
    """
    Chain collinear seed matches using dynamic programming.
    
    For each match, either extend the chain ending at the previous match
    (if collinear) or start a new chain.
    
    Returns list of (q_start, q_end, r_start, r_end, chain_score) tuples.
    """
    if not matches:
        return []
    
    # Sort by query position, then by reference position
    matches.sort(key=lambda x: (x[0], x[1]))
    
    # DP: for each match, find the best previous match to extend
    # score[i] = best chain score ending at match i
    # prev[i] = index of previous match in chain
    
    n = len(matches)
    score = [1] * n
    prev = [-1] * n
    
    for i in range(n):
        qi, ri = matches[i]
        for j in range(i):
            qj, rj = matches[j]
            gap_q = qi - qj
            gap_r = ri - rj
            
            # Check collinearity: gaps should be similar
            if 0 <= gap_q - gap_r <= max_gap and gap_q >= 0:
                new_score = score[j] + 1
                if new_score > score[i]:
                    score[i] = new_score
                    prev[i] = j
    
    # Find best ending match
    if not score:
        return []
    
    best_end = max(range(n), key=lambda i: score[i])
    
    # Backtrack to find chain
    chain = []
    i = best_end
    while i != -1:
        chain.append(matches[i])
        i = prev[i]
    
    chain.reverse()
    
    if len(chain) < 10:
        return []
    
    q_start = chain[0][0]
    q_end = chain[-1][0] + 1
    r_start = chain[0][1]
    r_end = chain[-1][1] + 1
    
    return [(q_start, q_end, r_start, r_end, len(chain))]


def find_best_alignment(query, reference, k=21, w=10, min_seeds=30):
    """
    Find best alignment between query and reference using seed-and-extend.
    
    Returns (q_start, q_end, r_start, r_end, identity, coverage) or None.
    """
    ref_len = len(reference)
    query_len = len(query)
    
    if ref_len < 10000 or query_len < 10000:
        return None
    
    # Find seed matches
    matches = find_seed_matches(query, reference, k, w)
    
    if len(matches) < min_seeds:
        return None
    
    # Chain seeds
    chains = chain_seeds(matches, max_gap=5000)
    
    if not chains:
        return None
    
    # Evaluate each chain
    best_result = None
    best_score = 0
    
    for q_start, q_end, r_start, r_end, num_seeds in chains:
        aligned_len = q_end - q_start
        coverage = aligned_len / query_len
        
        # Count exact matches in aligned region
        exact_matches = 0
        for i in range(min(aligned_len, 100000)):  # Sample up to 100kb
            if query[q_start + i] == reference[r_start + i]:
                exact_matches += 1
        
        identity = exact_matches / max(1, aligned_len)
        
        # Score: prefer high identity and good coverage
        score = identity * coverage * num_seeds
        
        if score > best_score and identity >= 0.95 and coverage >= 0.01:
            best_score = score
            best_result = (q_start, q_end, r_start, r_end, identity, coverage)
    
    return best_result


# ─── Main SNP Extraction ─────────────────────────────────────────────────────

def extract_snps_direct(camilow_path, k=21, w=10):
    """Extract SNPs by comparing evo genomes vs parent genomes."""
    genomes_dir = Path(camilow_path) / 'source_genomes_low' / 'source_genomes'
    
    if not genomes_dir.exists():
        print(f"ERROR: Genomes directory not found: {genomes_dir}")
        return {}
    
    # Identify parent-evo relationships
    parent_map = {}
    parent_files = []
    evo_files = []
    
    for f in genomes_dir.iterdir():
        if f.suffix in ('.fna', '.fasta'):
            name = f.name
            if name.startswith('evo_'):
                evo_files.append(f)
            elif 'run' in name or (name[0].isdigit() and '.gt1kb.fasta' in name):
                parent_files.append(f)
    
    # Build parent strain prefix map
    for pf in parent_files:
        name = pf.name
        if 'run' in name:
            base = name.split('_run')[0].split('.')[0]
            parent_map[base] = pf
        elif '.gt1kb.fasta' in name:
            base = name.split('.gt1kb')[0]
            parent_map[base] = pf
    
    print(f"Found {len(parent_files)} parent genomes")
    print(f"Found {len(evo_files)} evo genomes")
    print(f"Parent map: {len(parent_map)} entries")
    
    # Map evo strains to parent genomes
    evo_to_parent = {}
    for ef in evo_files:
        name = ef.name
        base = name.replace('.fna', '').replace('.fasta', '')
        strain_id = base[4:]  # Remove "evo_" prefix
        
        parts = strain_id.rsplit('.', 1)
        if len(parts) == 2:
            parent_prefix = parts[0]
            strain_suffix = parts[1]
            
            if parent_prefix in parent_map:
                evo_to_parent[strain_id] = {
                    'parent_file': parent_map[parent_prefix],
                    'parent_prefix': parent_prefix,
                    'strain_suffix': strain_suffix,
                    'evo_file': ef
                }
    
    print(f"Mapped {len(evo_to_parent)} evo strains to parents\n")
    
    # Load parent genomes
    parent_genomes_cache = {}
    for strain_id, mapping in evo_to_parent.items():
        pf = mapping['parent_file']
        parent_prefix = mapping['parent_prefix']
        
        if parent_prefix not in parent_genomes_cache:
            file_size_mb = pf.stat().st_size / (1024 * 1024)
            print(f"  Loading parent: {pf.name} ({file_size_mb:.1f}MB)")
            parent_genomes_cache[parent_prefix] = read_fasta(pf)
    
    # Extract SNPs for each evo strain
    all_snps = {}
    
    for strain_id, mapping in sorted(evo_to_parent.items()):
        evo_file = mapping['evo_file']
        parent_prefix = mapping['parent_prefix']
        
        print(f"\nProcessing: {strain_id}")
        print(f"  Evo file: {evo_file.name}")
        
        evo_genomes = read_fasta(evo_file)
        if not evo_genomes:
            print(f"  SKIP: No sequences in evo file")
            continue
        
        evo_contig_name = list(evo_genomes.keys())[0]
        evo_seq = evo_genomes[evo_contig_name]
        print(f"  Evo contig: {evo_contig_name}, length: {len(evo_seq)}")
        
        parent_scaffolds = parent_genomes_cache.get(parent_prefix, {})
        if not parent_scaffolds:
            print(f"  SKIP: No parent scaffolds loaded")
            continue
        
        print(f"  Parent scaffolds: {len(parent_scaffolds)}")
        
        # Find best parent scaffold
        best_result = None
        best_score = 0
        best_scaffold_name = None
        
        for scaffold_name, scaffold_seq in parent_scaffolds.items():
            if len(scaffold_seq) < 10000:
                continue
            
            result = find_best_alignment(evo_seq, scaffold_seq, k, w)
            
            if result:
                q_start, q_end, r_start, r_end, identity, coverage = result
                score = identity * coverage
                
                if score > best_score:
                    best_score = score
                    best_result = result
                    best_scaffold_name = scaffold_name
        
        if best_result is None:
            print(f"  SKIP: No alignment found")
            continue
        
        q_start, q_end, r_start, r_end, identity, coverage = best_result
        aligned_length = q_end - q_start
        
        print(f"  Best scaffold: {best_scaffold_name}")
        print(f"  Alignment: query[{q_start}:{q_end}] vs parent[{r_start}:{r_end}]")
        print(f"  Aligned length: {aligned_length}bp")
        print(f"  Identity: {identity*100:.4f}%")
        print(f"  Coverage: {coverage*100:.1f}%")
        
        # Extract SNPs
        snps = []
        min_len = min(aligned_length, len(evo_seq) - q_start, len(parent_scaffolds[best_scaffold_name]) - r_start)
        
        for i in range(min_len):
            q_base = evo_seq[q_start + i]
            r_base = parent_scaffolds[best_scaffold_name][r_start + i]
            
            if q_base != r_base and q_base != 'N' and r_base != 'N':
                snps.append({
                    'evo_position': q_start + i,
                    'parent_position': r_start + i,
                    'scaffold': best_scaffold_name,
                    'strain': strain_id,
                    'ref_base': r_base,
                    'alt_base': q_base
                })
        
        if snps:
            print(f"  SNPs found: {len(snps)}")
            all_snps[strain_id] = snps
        else:
            print(f"  No SNPs found (perfect match)")
    
    return all_snps


# ─── Statistics ───────────────────────────────────────────────────────────────

def compute_statistics(all_snps):
    """Compute summary statistics from SNP data."""
    stats = {
        'total_strains': len(all_snps),
        'total_snps': sum(len(snps) for snps in all_snps.values()),
        'strains_with_snps': sum(1 for snps in all_snps.values() if len(snps) > 0),
        'strain_details': []
    }
    
    for strain_id in sorted(all_snps.keys()):
        snps = all_snps[strain_id]
        if snps:
            base_changes = defaultdict(int)
            for snp in snps:
                change = f"{snp['ref_base']}>{snp['alt_base']}"
                base_changes[change] += 1
            
            stats['strain_details'].append({
                'strain': strain_id,
                'snp_count': len(snps),
                'base_changes': dict(base_changes)
            })
    
    return stats


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Direct SNP extraction v2: seed-and-extend')
    parser.add_argument('--camilow', default='../CAMI_low', help='Path to CAMI_low directory')
    parser.add_argument('--k', type=int, default=21, help='K-mer size for seeding')
    parser.add_argument('--w', type=int, default=10, help='Window size for minimizers')
    parser.add_argument('--output', default='snp_profiles_direct.json', help='Output JSON file')
    parser.add_argument('--stats', default='snp_stats_direct.txt', help='Output statistics file')
    args = parser.parse_args()
    
    camilow_path = Path(args.camilow)
    if not camilow_path.exists():
        print(f"ERROR: CAMI_low directory not found: {camilow_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("Step 0: Direct SNP Extraction v2 (seed-and-extend)")
    print("=" * 60)
    print(f"CAMI_low path: {camilow_path.absolute()}")
    print(f"K-mer size: {args.k}, Window: {args.w}")
    print()
    
    start_time = time.time()
    
    all_snps = extract_snps_direct(camilow_path, k=args.k, w=args.w)
    
    elapsed = time.time() - start_time
    
    stats = compute_statistics(all_snps)
    
    output_data = {
        'snps': all_snps,
        'statistics': stats,
        'parameters': {
            'kmer_size': args.k,
            'window_size': args.w,
            'extraction_time_seconds': round(elapsed, 2)
        }
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    with open(args.stats, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("Direct SNP Extraction Statistics (v2)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total strains analyzed: {stats['total_strains']}\n")
        f.write(f"Total SNPs found: {stats['total_snps']}\n")
        f.write(f"Strains with SNPs: {stats['strains_with_snps']}\n")
        f.write(f"Extraction time: {elapsed:.2f}s\n\n")
        
        f.write("-" * 60 + "\n")
        f.write("Per-Strain Breakdown\n")
        f.write("-" * 60 + "\n\n")
        
        for detail in stats['strain_details']:
            f.write(f"  {detail['strain']}: {detail['snp_count']} SNPs\n")
            for change, count in sorted(detail['base_changes'].items()):
                f.write(f"    {change}: {count}\n")
            f.write("\n")
    
    print(f"\n{'=' * 60}")
    print(f"Results: {args.output}")
    print(f"Stats: {args.stats}")
    print(f"Time: {elapsed:.2f}s")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
