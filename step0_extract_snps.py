"""
Step 0: Direct SNP Extraction — evo vs parent genome comparison.

Compares each evo genome sequence against its parent genome using sliding window
to find strain-specific SNPs. This is the foundation for SNP-aware classification.

Usage:
    python step0_extract_snps.py [options]

Options:
    --camilow PATH       Path to CAMI_low directory (default: ../CAMI_low)
    --k MER_SIZE         K-mer size for exact matching (default: 15)
    --min-support N      Minimum support for SNP calls (default: 2)
    --output FILE        Output SNP profiles JSON (default: snp_profiles_direct.json)
    --stats FILE         Output statistics file (default: snp_stats_direct.txt)
    --max-parent-size MB Max parent genome size in MB to load (default: 50)
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
                current_header = line[1:].split()[0]  # Take first word only
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


# ─── Sliding Window Exact Match ──────────────────────────────────────────────

def find_exact_matches(query, reference, k=25):
    """
    Find all exact k-mer matches between query and reference.
    Uses non-overlapping k-mers for speed (step=k).
    Returns list of (query_pos, ref_pos) pairs.
    """
    matches = []
    ref_len = len(reference)
    query_len = len(query)
    
    # Build index of k-mers in reference (non-overlapping, step=1 for precision)
    # But only store positions for k-mers that appear <= 100 times (avoid repetitive)
    kmer_positions = defaultdict(list)
    for i in range(ref_len - k + 1):
        kmer = reference[i:i+k]
        if len(kmer_positions[kmer]) < 100:  # Filter repetitive k-mers
            kmer_positions[kmer].append(i)
    
    # Find matches in query (non-overlapping, step=k for speed)
    for i in range(0, query_len - k + 1, k):
        kmer = query[i:i+k]
        if kmer in kmer_positions and len(kmer_positions[kmer]) < 10:
            for ref_pos in kmer_positions[kmer]:
                matches.append((i, ref_pos))
    
    return matches


def find_best_local_alignment(query, reference, k=25, min_identity=0.95):
    """
    Find the best local alignment between query and reference.
    
    Strategy: scan reference with a sliding window, find region with highest
    k-mer identity to the query. For >99.9% identical genomes, this finds
    the correct alignment even with assembly artifacts.
    
    Returns (q_start, q_end, r_start, r_end, identity, coverage) or None.
    """
    ref_len = len(reference)
    query_len = len(query)
    
    if ref_len < 10000 or query_len < 10000:
        return None
    
    # Use a large window size for scanning (100kb)
    window_size = min(100000, ref_len // 4)
    step_size = window_size // 4
    
    best_result = None
    best_identity = 0.0
    
    # Build query k-mer set (non-overlapping)
    query_kmers = set()
    for i in range(0, query_len - k + 1, k):
        query_kmers.add(query[i:i+k])
    
    # Scan reference windows
    for r_start in range(0, ref_len - window_size + 1, step_size):
        r_end = r_start + window_size
        ref_window = reference[r_start:r_end]
        
        # Count matching k-mers
        matches = 0
        total = 0
        
        for i in range(0, len(ref_window) - k + 1, k):
            total += 1
            kmer = ref_window[i:i+k]
            if kmer in query_kmers:
                matches += 1
        
        identity = matches / max(1, total)
        
        if identity > best_identity:
            best_identity = identity
            best_result = (r_start, r_end)
            
            # Early stop if we find a very good match
            if identity > 0.99:
                break
    
    if best_result is None or best_identity < min_identity:
        return None
    
    r_start, r_end = best_result
    window_len = r_end - r_start
    
    # Now refine: find exact boundaries using bidirectional extension
    # Start from the window center
    center_r = (r_start + r_end) // 2
    center_q = center_r  # Assume collinear (diagonal ~0)
    
    # Extend left
    left_r = center_r
    left_q = center_q
    while left_r > 0 and left_q > 0 and left_r - 1000 < left_q:
        left_r -= 1
        left_q -= 1
    
    # Extend right
    right_r = center_r
    right_q = center_q
    while right_r < ref_len - 1 and right_q < query_len - 1:
        right_r += 1
        right_q += 1
    
    # Calculate final alignment stats
    final_q_start = max(0, left_q)
    final_q_end = min(query_len, right_q)
    final_r_start = max(0, left_r)
    final_r_end = min(ref_len, right_r)
    
    aligned_len = final_q_end - final_q_start
    if aligned_len < 1000:
        return None
    
    # Count exact matches in refined region
    exact_matches = 0
    for i in range(aligned_len):
        if query[final_q_start + i] == reference[final_r_start + i]:
            exact_matches += 1
    
    final_identity = exact_matches / aligned_len
    
    if final_identity < min_identity:
        return None
    
    coverage = aligned_len / query_len
    
    return (final_q_start, final_q_end, final_r_start, final_r_end, final_identity, coverage)


def extract_snps_from_alignment(query, reference, q_start, r_start, k=15):
    """
    Compare query vs reference within aligned region.
    Returns list of SNP dicts: {position, ref_base, alt_base, support}
    """
    snps = []
    length = min(len(query) - q_start, len(reference) - r_start)
    
    for i in range(length):
        q_base = query[q_start + i]
        r_base = reference[r_start + i]
        
        if q_base != r_base and q_base != 'N' and r_base != 'N':
            snps.append({
                'evo_position': q_start + i,
                'parent_position': r_start + i,
                'ref_base': r_base,
                'alt_base': q_base
            })
    
    return snps


# ─── Main SNP Extraction ─────────────────────────────────────────────────────

def extract_snps_direct(camilow_path, k=15, max_parent_size_mb=50):
    """
    Extract SNPs by directly comparing evo genomes vs parent genomes.
    
    Returns: dict of {strain_name: [snp_list]}
    """
    genomes_dir = Path(camilow_path) / 'source_genomes_low' / 'source_genomes'
    
    if not genomes_dir.exists():
        print(f"ERROR: Genomes directory not found: {genomes_dir}")
        return {}
    
    # Identify parent-evo relationships
    # Parent genomes: contain "run" in name (e.g., 1286_AP_run191_run197.final.scaffolds.gt1kb.fasta)
    # or numeric IDs (e.g., 1035930.gt1kb.fasta)
    # Evo genomes: start with "evo_" (e.g., evo_1286_AP.026.fna)
    
    parent_map = {}  # Maps strain prefix -> parent genome file
    
    # First, find all parent genomes
    parent_files = []
    evo_files = []
    
    for f in genomes_dir.iterdir():
        if f.suffix == '.fna' or f.suffix == '.fasta':
            name = f.name
            if name.startswith('evo_'):
                evo_files.append(f)
            elif 'run' in name or (name[0].isdigit() and '.gt1kb.fasta' in name):
                parent_files.append(f)
    
    # Build parent strain prefix map
    # Parent: 1286_AP_run191_run197.final.scaffolds.gt1kb.fasta -> "1286_AP"
    # Parent: 1035930.gt1kb.fasta -> "1035930"
    for pf in parent_files:
        name = pf.name
        # Extract base name
        if 'run' in name:
            # 1286_AP_run191_run197.final.scaffolds.gt1kb.fasta -> 1286_AP
            base = name.split('_run')[0]
            # Remove trailing dots/characters
            base = base.split('.')[0]
            parent_map[base] = pf
        elif '.gt1kb.fasta' in name:
            # 1035930.gt1kb.fasta -> 1035930
            base = name.split('.gt1kb')[0]
            parent_map[base] = pf
    
    print(f"Found {len(parent_files)} parent genomes")
    print(f"Found {len(evo_files)} evo genomes")
    print(f"Parent map: {len(parent_map)} entries")
    
    # Map evo strains to parent genomes
    evo_to_parent = {}
    for ef in evo_files:
        name = ef.name  # evo_1286_AP.026.fna
        # Extract strain identifier: evo_1286_AP.026
        base = name.replace('.fna', '').replace('.fasta', '')
        # Remove "evo_" prefix
        strain_id = base[4:]  # 1286_AP.026
        
        # Find parent: everything before the last ".NNN"
        parts = strain_id.rsplit('.', 1)
        if len(parts) == 2:
            parent_prefix = parts[0]  # 1286_AP
            strain_suffix = parts[1]  # 026
            
            # Look up parent
            if parent_prefix in parent_map:
                evo_to_parent[strain_id] = {
                    'parent_file': parent_map[parent_prefix],
                    'parent_prefix': parent_prefix,
                    'strain_suffix': strain_suffix,
                    'evo_file': ef
                }
            else:
                print(f"  WARNING: No parent found for {strain_id} (prefix: {parent_prefix})")
        else:
            print(f"  WARNING: Could not parse evo strain ID: {strain_id}")
    
    print(f"\nMapped {len(evo_to_parent)} evo strains to parents\n")
    
    # Load parent genomes (with size limit)
    parent_genomes_cache = {}
    for strain_id, mapping in evo_to_parent.items():
        pf = mapping['parent_file']
        parent_prefix = mapping['parent_prefix']
        
        if parent_prefix not in parent_genomes_cache:
            # Check file size
            file_size_mb = pf.stat().st_size / (1024 * 1024)
            if file_size_mb > max_parent_size_mb:
                print(f"  SKIP: {pf.name} is {file_size_mb:.1f}MB (max {max_parent_size_mb}MB)")
                continue
            
            print(f"  Loading parent: {pf.name} ({file_size_mb:.1f}MB)")
            parent_genomes_cache[parent_prefix] = read_fasta(pf)
    
    # Extract SNPs for each evo strain
    all_snps = {}
    
    for strain_id, mapping in sorted(evo_to_parent.items()):
        evo_file = mapping['evo_file']
        parent_prefix = mapping['parent_prefix']
        
        print(f"\nProcessing: {strain_id}")
        print(f"  Evo file: {evo_file.name}")
        
        # Load evo genome
        evo_genomes = read_fasta(evo_file)
        if not evo_genomes:
            print(f"  SKIP: No sequences in evo file")
            continue
        
        evo_contig_name = list(evo_genomes.keys())[0]
        evo_seq = evo_genomes[evo_contig_name]
        print(f"  Evo contig: {evo_contig_name}, length: {len(evo_seq)}")
        
        # Get parent scaffolds
        parent_scaffolds = parent_genomes_cache.get(parent_prefix, {})
        if not parent_scaffolds:
            print(f"  SKIP: No parent scaffolds loaded")
            continue
        
        print(f"  Parent scaffolds: {len(parent_scaffolds)}")
        
        # Find best parent scaffold for this evo strain
        best_result = None
        best_coverage = 0.0
        best_scaffold_name = None
        best_identity = 0.0
        
        for scaffold_name, scaffold_seq in parent_scaffolds.items():
            if len(scaffold_seq) < 10000:
                continue  # Skip very short scaffolds
            
            # Find best local alignment
            result = find_best_local_alignment(evo_seq, scaffold_seq, k=25, min_identity=0.95)
            
            if result:
                q_start, q_end, r_start, r_end, identity, coverage = result
                
                # Prefer highest identity, then coverage
                if identity > best_identity or (identity == best_identity and coverage > best_coverage):
                    best_identity = identity
                    best_coverage = coverage
                    best_result = (q_start, q_end, r_start, r_end)
                    best_scaffold_name = scaffold_name
        
        if best_result is None:
            print(f"  SKIP: No high-quality alignment found (need >=95% identity)")
            continue
        
        q_start, q_end, r_start, r_end = best_result
        aligned_length = q_end - q_start
        
        print(f"  Best scaffold: {best_scaffold_name}")
        print(f"  Alignment: query[{q_start}:{q_end}] vs parent[{r_start}:{r_end}]")
        print(f"  Aligned length: {aligned_length}bp")
        print(f"  Identity: {best_identity*100:.2f}%")
        print(f"  Coverage: {best_coverage*100:.1f}%")
        
        # Extract SNPs from aligned region
        # Get evo sequence portion and parent sequence portion
        evo_aligned_seq = evo_seq[q_start:q_end]
        parent_aligned_seq = parent_scaffolds[best_scaffold_name][r_start:r_end]
        
        # Compare base by base
        snps = []
        min_len = min(len(evo_aligned_seq), len(parent_aligned_seq))
        
        for i in range(min_len):
            q_base = evo_aligned_seq[i]
            r_base = parent_aligned_seq[i]
            
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
            # Add scaffold info
            for snp in snps:
                snp['scaffold'] = best_scaffold_name
                snp['strain'] = strain_id
            
            # Count SNPs
            identity = 100.0 - (len(snps) / aligned_length * 100)
            print(f"  SNPs found: {len(snps)}")
            print(f"  Identity: {identity:.4f}%")
            
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
            # Count base changes
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
    parser = argparse.ArgumentParser(description='Direct SNP extraction: evo vs parent genome comparison')
    parser.add_argument('--camilow', default='../CAMI_low', help='Path to CAMI_low directory')
    parser.add_argument('--k', type=int, default=15, help='K-mer size for exact matching')
    parser.add_argument('--min-support', type=int, default=2, help='Minimum support for SNP calls')
    parser.add_argument('--output', default='snp_profiles_direct.json', help='Output JSON file')
    parser.add_argument('--stats', default='snp_stats_direct.txt', help='Output statistics file')
    parser.add_argument('--max-parent-size', type=int, default=50, help='Max parent genome size in MB')
    args = parser.parse_args()
    
    camilow_path = Path(args.camilow)
    if not camilow_path.exists():
        print(f"ERROR: CAMI_low directory not found: {camilow_path}")
        print(f"Current working directory: {Path.cwd()}")
        sys.exit(1)
    
    print("=" * 60)
    print("Step 0: Direct SNP Extraction")
    print("=" * 60)
    print(f"CAMI_low path: {camilow_path.absolute()}")
    print(f"K-mer size: {args.k}")
    print(f"Max parent size: {args.max_parent_size}MB")
    print()
    
    start_time = time.time()
    
    # Extract SNPs
    all_snps = extract_snps_direct(
        camilow_path,
        k=args.k,
        max_parent_size_mb=args.max_parent_size
    )
    
    elapsed = time.time() - start_time
    
    # Compute statistics
    stats = compute_statistics(all_snps)
    
    # Save results
    output_data = {
        'snps': all_snps,
        'statistics': stats,
        'parameters': {
            'kmer_size': args.k,
            'max_parent_size_mb': args.max_parent_size,
            'extraction_time_seconds': round(elapsed, 2)
        }
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    # Save stats
    with open(args.stats, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("Direct SNP Extraction Statistics\n")
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
    print(f"Results saved to: {args.output}")
    print(f"Statistics saved to: {args.stats}")
    print(f"Total time: {elapsed:.2f}s")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
