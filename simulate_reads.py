#!/usr/bin/env python3
"""Simple read simulator for DNA sequencing."""

import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 2-bit encoding: A=00, C=01, G=10, T=11
BASES = ['A', 'C', 'G', 'T']
MUTATIONS = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 
             'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}

def simulate_errors(seq, error_rate, quality_base=30):
    """Introduce substitutions based on error rate, return (mutated_seq, qualities)."""
    seq = str(seq).upper()
    mutated = []
    qualities = []
    
    for base in seq:
        if random.random() < error_rate:
            # Mutation
            new_base = random.choice(MUTATIONS[base])
            mutated.append(new_base)
            # Lower quality for mutated bases
            qual = max(10, quality_base - random.randint(5, 15))
        else:
            mutated.append(base)
            # High quality for correct bases
            qual = random.randint(quality_base, quality_base + 10)
        qualities.append(qual)
    
    return ''.join(mutated), qualities

def simulate_reads(genome_fasta, num_reads, read_len, error_rate, 
                   output_fastq, paired_end=False, insert_size=300, insert_std=50):
    """Simulate sequencing reads from a genome."""
    # Load genome
    record = SeqIO.read(genome_fasta, "fasta-pearson")
    genome = str(record.seq).upper()
    genome_len = len(genome)
    
    print(f"Genome: {record.id} ({genome_len} bp)")
    print(f"Simulating {num_reads} reads of {read_len} bp")
    print(f"Error rate: {error_rate:.4f}")
    
    reads = []
    for i in range(num_reads):
        # Random position
        pos = random.randint(0, genome_len - read_len)
        seq = genome[pos:pos + read_len]
        
        # Add errors
        mutated_seq, qualities = simulate_errors(seq, error_rate)
        
        name = f"sim_{record.id}_{i+1}"
        reads.append((name, mutated_seq, qualities))
    
    # Write FASTQ
    with open(output_fastq, 'w') as f:
        for name, seq, quals in reads:
            f.write(f"@{name}\n{seq}\n+\n{''.join(chr(q + 33) for q in quals)}\n")
    
    print(f"Wrote {len(reads)} reads to {output_fastq}")
    
    # Paired-end
    if paired_end:
        output_r2 = output_fastq.replace('.fastq', '_R2.fastq')
        reads_r2 = []
        for i in range(num_reads):
            # Forward read position
            pos_f = random.randint(0, genome_len - read_len)
            # Reverse read position (insert_size away)
            pos_r = min(pos_f + insert_size + random.gauss(0, insert_std), genome_len - read_len)
            pos_r = max(pos_r, read_len)
            
            # Forward read
            seq_f = genome[pos_f:pos_f + read_len]
            mutated_f, quals_f = simulate_errors(seq_f, error_rate)
            
            # Reverse read (from opposite strand)
            seq_r = genome[pos_r:pos_r + read_len]
            seq_r_revcomp = str(Seq(seq_r).reverse_complement())
            mutated_r, quals_r = simulate_errors(seq_r_revcomp, error_rate)
            
            name = f"sim_{record.id}_{i+1}"
            reads_r2.append((f"{name}/2", mutated_r, quals_r))
        
        with open(output_r2, 'w') as f:
            for name, seq, quals in reads_r2:
                f.write(f"@{name}\n{seq}\n+\n{''.join(chr(q + 33) for q in quals)}\n")
        
        print(f"Wrote {len(reads_r2)} R2 reads to {output_r2}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate DNA sequencing reads')
    parser.add_argument('genome', help='Input FASTA genome')
    parser.add_argument('-n', '--num-reads', type=int, default=10000, help='Number of reads')
    parser.add_argument('-l', '--read-len', type=int, default=100, help='Read length')
    parser.add_argument('-e', '--error-rate', type=float, default=0.001, help='Error rate (default: 0.1%)')
    parser.add_argument('-q', '--quality', type=int, default=30, help='Base quality score')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    parser.add_argument('--paired', action='store_true', help='Paired-end mode')
    parser.add_argument('--insert-size', type=int, default=300, help='Insert size (paired-end)')
    parser.add_argument('--insert-std', type=int, default=50, help='Insert size std dev')
    
    args = parser.parse_args()
    simulate_reads(args.genome, args.num_reads, args.read_len, args.error_rate,
                   args.output, args.paired, args.insert_size, args.insert_std)
