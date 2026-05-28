#!/usr/bin/env python3
"""Download bacterial genomes from NCBI for benchmark."""
import os, sys, time
from Bio import Entrez, SeqIO

Entrez.email = "mladen.popovic@gmail.com"
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "./benchmark_genomes"
os.makedirs(OUTDIR, exist_ok=True)

# Read accessions from file
with open("benchmark_accessions.txt") as f:
    ACC = [line.strip() for line in f if line.strip() and not line.startswith('#')]

print(f"Downloading {len(ACC)} genomes...")

for i, acc in enumerate(ACC):
    outpath = os.path.join(OUTDIR, f"{acc}.fa")
    if os.path.exists(outpath):
        print(f"[{i+1}/{len(ACC)}] SKIP {acc} (exists)")
        continue
    try:
        print(f"[{i+1}/{len(ACC)}] FETCH {acc}...", end=" ", flush=True)
        handle = Entrez.efetch(db='nuccore', id=acc, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        handle.close()
        SeqIO.write(record, outpath, 'fasta')
        print(f"OK ({len(record.seq)} bp)")
        time.sleep(1)  # Be nice to NCBI
    except Exception as e:
        print(f"FAIL: {e}")

print(f"\nDone. Check {OUTDIR}/")
