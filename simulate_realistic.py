#!/usr/bin/env python3
"""Simulate realistic PacBio HiFi reads - optimized."""
import os, random, sys
from Bio import SeqIO

GENOMES_DIR = "G:\\TEST\\BENCHMARK70\\genomes"
OUTDIR = "G:\\TEST\\BENCHMARK70\\reads"
os.makedirs(OUTDIR, exist_ok=True)

READ_LEN_MIN = 8000
READ_LEN_MAX = 20000
ERROR_RATE = 0.001
HOMOPOLYMER_ERROR_RATE = 0.02
CHIMERA_RATE = 0.01
COVERAGE_STD = 0.5

BASES = {'A', 'C', 'G', 'T'}
MUTATIONS = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'],
             'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}

genome_files = sorted([f for f in os.listdir(GENOMES_DIR) if f.endswith('.fa')])
print(f"Found {len(genome_files)} genomes")

# Calculate target reads per genome with coverage variation
total_target = 76500
genome_sizes = []

for gfile in genome_files:
    gpath = os.path.join(GENOMES_DIR, gfile)
    if os.path.getsize(gpath) < 100000:
        continue
    try:
        record = SeqIO.read(gpath, 'fasta')
    except:
        record = SeqIO.read(gpath, 'fasta-pearson')
    genome_sizes.append((gfile, len(str(record.seq).upper())))

# Calculate reads per genome with coverage variation
reads_per_genome = {}
for gfile, size in genome_sizes:
    variation = random.gauss(1.0, COVERAGE_STD)
    variation = max(0.3, min(2.0, variation))
    n_reads = int((total_target / len(genome_sizes)) * variation)
    reads_per_genome[gfile] = n_reads
    print(f"{gfile}: {size} bp, {n_reads} reads")

all_reads = []
ground_truth = []

for gfile, n_reads in reads_per_genome.items():
    gpath = os.path.join(GENOMES_DIR, gfile)
    try:
        record = SeqIO.read(gpath, 'fasta')
    except:
        record = SeqIO.read(gpath, 'fasta-pearson')
    
    genome = str(record.seq).upper()
    gname = gfile.replace('.fa', '')
    
    for i in range(n_reads):
        rlen = random.randint(READ_LEN_MIN, READ_LEN_MAX)
        if rlen > len(genome):
            rlen = len(genome) - 100
        pos = random.randint(0, len(genome) - rlen)
        seq = genome[pos:pos + rlen]
        
        # Fast error injection with homopolymer detection
        mutated = []
        quals = []
        run_len = 0
        prev_base = ''
        
        for base in seq:
            if base not in BASES:
                mutated.append(base)
                quals.append(30)
                run_len = 0
                prev_base = base
                continue
            
            # Track homopolymer runs
            if base == prev_base:
                run_len += 1
            else:
                run_len = 1
                prev_base = base
            
            # Higher error in homopolymers (run >= 4)
            err_rate = HOMOPOLYMER_ERROR_RATE if run_len >= 4 else ERROR_RATE
            
            if random.random() < err_rate:
                mutated.append(random.choice(MUTATIONS[base]))
                quals.append(15 if run_len >= 4 else 20)
            else:
                mutated.append(base)
                quals.append(35)
        
        read_seq = ''.join(mutated)
        read_qual = ''.join(chr(q + 33) for q in quals)
        
        # Chimera: 1% chance
        is_chimera = False
        if random.random() < CHIMERA_RATE and len(genome) > rlen * 2:
            is_chimera = True
            pos2 = random.randint(0, len(genome) - rlen // 2)
            seq2 = genome[pos2:pos2 + rlen // 2]
            
            # Simple error injection for chimera part
            m2 = []
            q2 = []
            for base in seq2:
                if base in BASES and random.random() < ERROR_RATE:
                    m2.append(random.choice(MUTATIONS[base]))
                    q2.append(20)
                else:
                    m2.append(base)
                    q2.append(35)
            
            read_seq = read_seq[:rlen//2] + ''.join(m2)
            read_qual = read_qual[:rlen//2] + ''.join(chr(q + 33) for q in q2)
        
        read_name = f"sim_{gname}_{i+1}"
        if is_chimera:
            read_name += "_chimera"
        
        all_reads.append((read_name, read_seq, read_qual))
        ground_truth.append(f"{read_name}\t{gname}")

# Shuffle reads
random.shuffle(all_reads)

# Write FASTQ
fastq_path = os.path.join(OUTDIR, "simulated_realistic.fastq")
with open(fastq_path, 'w') as f:
    for name, seq, qual in all_reads:
        f.write(f"@{name}\n{seq}\n+\n{qual}\n")

# Write ground truth
gt_path = os.path.join(OUTDIR, "ground_truth_realistic.tsv")
with open(gt_path, 'w') as f:
    f.write("read_name\tgenome\n")
    for line in ground_truth:
        f.write(line + "\n")

print(f"\nDone!")
print(f"FASTQ: {fastq_path} ({len(all_reads)} reads)")
print(f"Ground truth: {gt_path}")
