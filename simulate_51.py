#!/usr/bin/env python3
"""Simulate PacBio HiFi reads from 51 genomes - optimized."""
import os, random, sys
from Bio import SeqIO
from Bio.Seq import Seq

GENOMES_DIR = "G:\\TEST\\BENCHMARK70\\genomes"
OUTDIR = "G:\\TEST\\BENCHMARK70\\reads"
os.makedirs(OUTDIR, exist_ok=True)

READ_LEN_MIN = 8000
READ_LEN_MAX = 20000
ERROR_RATE = 0.001
READS_PER_GENOME = 1500

BASES = ['A', 'C', 'G', 'T']
MUTATIONS = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'],
             'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}

genome_files = sorted([f for f in os.listdir(GENOMES_DIR) if f.endswith('.fa')])
print(f"Found {len(genome_files)} genomes")

all_reads = []
ground_truth = []

for gfile in genome_files:
    gpath = os.path.join(GENOMES_DIR, gfile)
    if os.path.getsize(gpath) < 100000:
        continue
    try:
        record = SeqIO.read(gpath, 'fasta')
    except:
        record = SeqIO.read(gpath, 'fasta-pearson')
    
    genome = str(record.seq).upper()
    gname = gfile.replace('.fa', '')
    print(f"{gname}: {len(genome)} bp")
    
    for i in range(READS_PER_GENOME):
        rlen = random.randint(READ_LEN_MIN, READ_LEN_MAX)
        if rlen > len(genome):
            rlen = len(genome) - 100
        pos = random.randint(0, len(genome) - rlen)
        seq = genome[pos:pos + rlen]
        
        # Fast error injection
        mutated = list(seq)
        quals = [35] * len(mutated)
        for j in range(len(mutated)):
            if mutated[j] in BASES and random.random() < ERROR_RATE:
                mutated[j] = random.choice(MUTATIONS[mutated[j]])
                quals[j] = 20
        
        read_name = f"sim_{gname}_{i+1}"
        all_reads.append((read_name, ''.join(mutated), ''.join(chr(q + 33) for q in quals)))
        ground_truth.append(f"{read_name}\t{gname}")

# Shuffle reads
random.shuffle(all_reads)

# Write FASTQ
fastq_path = os.path.join(OUTDIR, "simulated_51genomes.fastq")
with open(fastq_path, 'w') as f:
    for name, seq, qual in all_reads:
        f.write(f"@{name}\n{seq}\n+\n{qual}\n")

# Write ground truth
gt_path = os.path.join(OUTDIR, "ground_truth.tsv")
with open(gt_path, 'w') as f:
    f.write("read_name\tgenome\n")
    for line in ground_truth:
        f.write(line + "\n")

print(f"\nDone!")
print(f"FASTQ: {fastq_path} ({len(all_reads)} reads)")
print(f"Ground truth: {gt_path}")
