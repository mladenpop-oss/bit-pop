import sys
from collections import defaultdict

sam_file = sys.argv[1] if len(sys.argv) > 1 else "D:/CAMI_low/mapped_se_tn4.sam"
gt_file = sys.argv[2] if len(sys.argv) > 2 else "D:/CAMI_low/gs_read_mapping.binning/gs_read_mapping.tsv"

print("=== Unmapped Reads Analysis ===")
print()

# Parse ground truth - format: READID\tBINID\tTAXID\tPOSITION
print("Parsing ground truth...")
gt = {}
with open(gt_file, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            read_id = parts[0]  # 1st column is read ID
            genome = parts[1]   # 2nd column is genome
            gt[read_id] = genome

print(f"  Loaded {len(gt)} unique ground truth reads")

# Parse SAM - get mapped read names (strip /1, /2 suffixes)
print("Parsing SAM output...")
mapped_names = set()
with open(sam_file, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        read_name = fields[0]
        # Strip /1, /2, :1, :2 suffixes
        if read_name.endswith('/1') or read_name.endswith('/2'):
            read_name = read_name[:-2]
        elif read_name.endswith(':1') or read_name.endswith(':2'):
            read_name = read_name[:-2]
        flag = int(fields[1])
        if flag & 4:
            continue
        mapped_names.add(read_name)

print(f"  Mapped reads: {len(mapped_names)}")

# Check overlap
common = gt.keys() & mapped_names
print(f"  Common reads: {len(common)}")
unmapped_gt = set(gt.keys()) - mapped_names
print(f"  Unmapped reads: {len(unmapped_gt)}")
print(f"  Unmapped rate: {len(unmapped_gt)/len(gt)*100:.2f}%")
print()

# Name format comparison
print("=== Name Format ===")
print(f"GT sample: {list(gt.keys())[:3]}")
print(f"SAM sample: {list(mapped_names)[:3]}")
print()

# Unmapped by genome
unmapped_by_genome = defaultdict(int)
for name in unmapped_gt:
    genome = gt[name]
    unmapped_by_genome[genome] += 1

print("=== Unmapped by Genome (Top 30) ===")
sorted_genomes = sorted(unmapped_by_genome.items(), key=lambda x: x[1], reverse=True)
for genome, count in sorted_genomes[:30]:
    gt_total = sum(1 for g in gt.values() if g == genome)
    unmapped_rate = count / gt_total * 100 if gt_total > 0 else 0
    print(f"  {genome:30s}  Unmapped: {count:7d}  ({unmapped_rate:5.1f}%)  [GT total: {gt_total}]")
