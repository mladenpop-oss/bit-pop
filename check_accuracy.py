#!/usr/bin/env python3
"""Measure mapping accuracy against ground truth."""
from collections import defaultdict

# Load ground truth
gt = {}
with open("G:\\TEST\\BENCHMARK70\\reads\\ground_truth.tsv") as f:
    next(f)  # skip header
    for line in f:
        read, genome = line.strip().split('\t')
        gt[read] = genome

# Parse SAM
correct = 0
wrong = 0
unmapped = 0
total = 0

with open("G:\\TEST\\BENCHMARK70\\mapped_k13.sam") as f:
    for line in f:
        if line.startswith('@'):
            continue
        parts = line.strip().split('\t')
        read_name = parts[0]
        flag = int(parts[1])
        genome = parts[2].split()[0]  # Get just accession
        
        if flag & 4:  # unmapped
            unmapped += 1
            continue
        
        total += 1
        expected = gt.get(read_name)
        if expected:
            if genome == expected:
                correct += 1
            else:
                wrong += 1

print(f"Total mapped: {total}")
print(f"Correct: {correct}")
print(f"Wrong: {wrong}")
print(f"Unmapped: {unmapped}")
print(f"Accuracy: {correct/total*100:.2f}%")
