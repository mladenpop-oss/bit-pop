"""
Step 0 v5: SNP Extraction - align EACH evo contig to parent scaffolds.

Evo genomes are PARTITIONS of parent genomes (e.g., contig_4_8 = partition 4 of 8).
Each evo contig is a subset of the parent, not the full genome.

Strategy: align each evo contig individually against parent scaffolds,
then merge SNP results per strain.
"""

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path


def read_fasta(filepath):
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


def find_best_scaffold_for_contig(evo_contig, parent_scaffolds, k=15):
    """
    Find the best parent scaffold for a single evo contig.
    
    Returns list of (scaffold_name, scaffold_seq, evo_start, evo_end, 
                     scaffold_start, scaffold_end, identity, coverage)
    sorted by identity descending.
    """
    evo_len = len(evo_contig)
    
    # Build evo k-mer set
    evo_kmers = set()
    for i in range(evo_len - k + 1):
        kmer = evo_contig[i:i+k]
        if 'N' not in kmer:
            evo_kmers.add(kmer)
    
    results = []
    
    for sname, sseq in parent_scaffolds.items():
        s_len = len(sseq)
        if s_len < 5000:
            continue
        
        # Sliding window scan
        window = min(50000, s_len // 5)
        step = window // 4
        
        best_id = 0
        best_pos = None
        
        for s_start in range(0, s_len - window + 1, step):
            s_end = s_start + window
            matches = 0
            total = 0
            
            for i in range(s_start, s_end - k + 1, k):
                kmer = sseq[i:i+k]
                if 'N' not in kmer:
                    total += 1
                    if kmer in evo_kmers:
                        matches += 1
            
            if total > 0:
                ident = matches / total
                if ident > best_id:
                    best_id = ident
                    best_pos = s_start
        
        if best_id >= 0.50 and best_pos is not None:
            # Refine: find exact boundaries
            s_start = best_pos
            # Expand left
            while s_start > 0:
                test = sseq[s_start-1:s_start+window]
                m = sum(1 for i in range(len(test)-k+1) if 'N' not in test[i:i+k] and test[i:i+k] in evo_kmers)
                t = sum(1 for i in range(len(test)-k+1) if 'N' not in test[i:i+k])
                if t > 0 and m/t >= 0.70:
                    s_start -= 1
                else:
                    break
            
            # Expand right
            s_end = best_pos + window
            while s_end < s_len:
                test = sseq[max(0,s_end-window):s_end+1]
                m = sum(1 for i in range(len(test)-k+1) if 'N' not in test[i:i+k] and test[i:i+k] in evo_kmers)
                t = sum(1 for i in range(len(test)-k+1) if 'N' not in test[i:i+k])
                if t > 0 and m/t >= 0.70:
                    s_end += 1
                else:
                    break
            
            # Calculate precise identity
            region = sseq[s_start:s_end]
            matches = 0
            total = 0
            for i in range(min(len(region), evo_len)):
                if region[i] != 'N' and i < evo_len:
                    total += 1
                    if evo_contig[i] == region[i]:
                        matches += 1
            
            identity = matches / max(1, total)
            coverage = min(len(region), evo_len) / evo_len
            
            results.append((sname, sseq, 0, min(len(region), evo_len), s_start, s_end, identity, coverage))
    
    results.sort(key=lambda x: x[6], reverse=True)
    return results


def extract_snps_v5(camilow_path, max_parent_size_mb=500):
    genomes_dir = Path(camilow_path) / 'source_genomes_low' / 'source_genomes'
    
    if not genomes_dir.exists():
        print(f"ERROR: {genomes_dir} not found")
        return {}
    
    # Identify files
    parent_files = []
    evo_files = []
    
    for f in genomes_dir.iterdir():
        if f.suffix in ('.fna', '.fasta'):
            name = f.name
            if name.startswith('evo_'):
                evo_files.append(f)
            elif 'run' in name or (name[0].isdigit() and '.gt1kb.fasta' in name):
                parent_files.append(f)
    
    # Build parent map
    parent_map = {}
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
    
    # Map evo to parent
    evo_to_parent = {}
    for ef in evo_files:
        name = ef.name
        base = name.replace('.fna', '').replace('.fasta', '')
        strain_id = base[4:]
        parts = strain_id.rsplit('.', 1)
        if len(parts) == 2:
            parent_prefix = parts[0]
            if parent_prefix in parent_map:
                evo_to_parent[strain_id] = {
                    'parent_file': parent_map[parent_prefix],
                    'parent_prefix': parent_prefix,
                    'evo_file': ef
                }
    
    print(f"Mapped {len(evo_to_parent)} evo strains to parents\n")
    
    # Process each strain
    all_snps = {}
    
    for strain_id, mapping in sorted(evo_to_parent.items()):
        evo_file = mapping['evo_file']
        parent_prefix = mapping['parent_prefix']
        
        print(f"\n{'='*60}")
        print(f"Strain: {strain_id}")
        print(f"{'='*60}")
        
        # Load evo (may have multiple contigs/partitions)
        evo_seqs = read_fasta(evo_file)
        print(f"  Evo contigs: {len(evo_seqs)}")
        for h, seq in evo_seqs.items():
            print(f"    {h}: {len(seq)}bp")
        
        # Load parent
        parent_file = mapping['parent_file']
        file_size_mb = parent_file.stat().st_size / (1024 * 1024)
        if file_size_mb > max_parent_size_mb:
            print(f"  SKIP: parent too large ({file_size_mb:.1f}MB)")
            continue
        
        parent_seqs = read_fasta(parent_file)
        print(f"  Parent scaffolds: {len(parent_seqs)}")
        
        # Align EACH evo contig against parent scaffolds
        all_strain_snps = []
        
        for contig_name, contig_seq in evo_seqs.items():
            print(f"\n  Aligning contig {contig_name} ({len(contig_seq)}bp)...")
            
            results = find_best_scaffold_for_contig(contig_seq, parent_seqs, k=15)
            
            if not results:
                print(f"    No scaffold match found")
                continue
            
            # Show top 3 matches
            for i, (sname, sseq, es, ee, ss, se, ident, cov) in enumerate(results[:3]):
                print(f"    #{i+1} {sname}: identity={ident*100:.2f}%, coverage={cov*100:.1f}%")
            
            # Extract SNPs from best match (if identity >= 80%)
            if results[0][6] >= 0.80:
                sname, sseq, es, ee, ss, se, identity, coverage = results[0]
                
                snps = []
                for i in range(es, ee):
                    if i < len(contig_seq) and ss + (i - es) < len(sseq):
                        evo_base = contig_seq[i]
                        par_base = sseq[ss + (i - es)]
                        if evo_base != par_base and evo_base != 'N' and par_base != 'N':
                            snps.append({
                                'evo_position': i,
                                'scaffold': sname,
                                'scaffold_pos': ss + (i - es),
                                'contig': contig_name,
                                'strain': strain_id,
                                'ref_base': par_base,
                                'alt_base': evo_base
                            })
                
                print(f"    SNPs from {contig_name}: {len(snps)} (identity={identity*100:.2f}%)")
                all_strain_snps.extend(snps)
            else:
                print(f"    Skipping (identity {results[0][6]*100:.1f}% < 80%)")
        
        if all_strain_snps:
            all_snps[strain_id] = all_strain_snps
            print(f"\n  TOTAL SNPs for {strain_id}: {len(all_strain_snps)}")
        else:
            print(f"\n  No SNPs found for {strain_id}")
    
    return all_snps


def compute_statistics(all_snps):
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


def main():
    parser = argparse.ArgumentParser(description='Step 0 v5: Per-contig SNP extraction')
    parser.add_argument('--camilow', default='../CAMI_low')
    parser.add_argument('--output', default='snp_profiles_v5.json')
    parser.add_argument('--stats', default='snp_stats_v5.txt')
    parser.add_argument('--max-parent-size', type=int, default=500)
    args = parser.parse_args()
    
    camilow_path = Path(args.camilow)
    if not camilow_path.exists():
        print(f"ERROR: {camilow_path} not found")
        sys.exit(1)
    
    print("=" * 60)
    print("Step 0 v5: SNP Extraction (Per-Contig Alignment)")
    print("=" * 60)
    print()
    
    start_time = time.time()
    all_snps = extract_snps_v5(camilow_path, args.max_parent_size)
    elapsed = time.time() - start_time
    
    stats = compute_statistics(all_snps)
    
    output_data = {
        'snps': all_snps,
        'statistics': stats,
        'parameters': {
            'method': 'per_contig_alignment_v5',
            'extraction_time_seconds': round(elapsed, 2)
        }
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    with open(args.stats, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("SNP Extraction Statistics (v5 - Per-Contig)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total strains: {stats['total_strains']}\n")
        f.write(f"Total SNPs: {stats['total_snps']}\n")
        f.write(f"Strains with SNPs: {stats['strains_with_snps']}\n")
        f.write(f"Time: {elapsed:.2f}s\n\n")
        
        for detail in stats['strain_details']:
            f.write(f"  {detail['strain']}: {detail['snp_count']} SNPs\n")
            for change, count in sorted(detail['base_changes'].items()):
                f.write(f"    {change}: {count}\n")
            f.write("\n")
    
    print(f"\n{'='*60}")
    print(f"Results: {args.output}")
    print(f"Stats: {args.stats}")
    print(f"Time: {elapsed:.2f}s")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
