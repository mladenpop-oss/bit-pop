import sys
from collections import defaultdict

sam_file = sys.argv[1] if len(sys.argv) > 1 else "D:\\CAMI_low\\mapped_base.sam"
gt_file = sys.argv[2] if len(sys.argv) > 2 else "D:\\CAMI_low\\gs_read_mapping.binning\\gs_read_mapping.tsv"

print("=== CAMI Accuracy Test ===")
print()

# Parse ground truth
print("Parsing ground truth...")
gt = {}
with open(gt_file, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            gt[parts[0]] = parts[1]
print(f"  Loaded {len(gt)} ground truth mappings")

# Parse SAM - keep only first mapping per read
print("Parsing SAM output...")
mapped = {}
lines = 0
with open(sam_file, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        lines += 1
        fields = line.split()
        if len(fields) < 3:
            continue
        
        read_name = fields[0]
        flag = int(fields[1])
        
        # Skip unmapped
        if flag & 4:
            continue
        
        # Strip /1 /2 suffix
        if read_name.endswith('/1') or read_name.endswith('/2'):
            read_name = read_name[:-2]
        
        genome = fields[2]
        
        # Keep first mapping per read
        if read_name not in mapped:
            mapped[read_name] = genome

print(f"  Processed {lines} SAM lines")
print(f"  Loaded {len(mapped)} unique read mappings")

# Compare - only reads that are in both datasets
print()
print("=== Accuracy Results ===")

correct = 0
wrong = 0
unmapped = 0
gt_reads = 0
wrong_predictions = defaultdict(lambda: defaultdict(int))

# Only compare reads that are in both ground truth and our SAM output
common_reads = set(gt.keys()) & set(mapped.keys())
print(f"Common reads (in both GT and SAM): {len(common_reads)}")
print()

for read_id in common_reads:
    gt_reads += 1
    expected = gt[read_id]
    predicted = mapped[read_id]
    
    if predicted == expected:
        correct += 1
    else:
        wrong += 1
        wrong_predictions[expected][predicted] += 1

print(f"Common reads:          {gt_reads}")
print()
print(f"Correct predictions:   {correct} ({correct/gt_reads*100:.2f}%)")
print(f"Wrong predictions:     {wrong} ({wrong/gt_reads*100:.2f}%)")

# Per-genome breakdown - only common reads
print()
print("=== Per-Genome Breakdown ===")
genome_stats = defaultdict(lambda: {'total': 0, 'correct': 0})

for read_id in common_reads:
    expected = gt[read_id]
    genome_stats[expected]['total'] += 1
    if mapped[read_id] == expected:
        genome_stats[expected]['correct'] += 1

print(f"{'Genome':<40} {'Total':>8} {'Correct':>10} {'Accuracy':>10}")
print("-" * 75)
for genome in sorted(genome_stats.keys()):
    s = genome_stats[genome]
    acc = s['correct'] / s['total'] * 100
    print(f"{genome:<40} {s['total']:>8} {s['correct']:>10} {acc:>9.2f}%")

# Wrong predictions breakdown
if wrong_predictions:
    print()
    print("=== Top Wrong Predictions ===")
    for expected in list(wrong_predictions.keys())[:10]:
        preds = wrong_predictions[expected]
        top = sorted(preds.items(), key=lambda x: x[1], reverse=True)[:3]
        print(f"{expected} -> {', '.join(f'{p}({c})' for p, c in top)}")
