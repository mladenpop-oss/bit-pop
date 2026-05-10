#!/usr/bin/env python3
"""Analyze discordant paired-end reads from CAMI benchmark."""

import sys
from collections import defaultdict


def parse_sam(sam_path):
    """Parse SAM file and extract R1/R2 mappings."""
    r1_maps = {}
    r2_maps = {}

    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            qname_full = parts[0]
            flag = int(parts[1])
            rname = parts[2]

            # Determine R1 or R2 from name suffix
            is_mapped = not bool(flag & 0x4)
            if not is_mapped:
                continue

            if qname_full.endswith('/1'):
                qname = qname_full[:-2]
                r1_maps[qname] = rname
            elif qname_full.endswith('/2'):
                qname = qname_full[:-2]
                r2_maps[qname] = rname

    return r1_maps, r2_maps


def parse_ground_truth(gt_path):
    """Parse ground truth TSV."""
    gt = {}
    with open(gt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            # Format: read_name genome
            read_name = parts[0].rstrip('/1').rstrip('/2')
            genome = parts[1]
            if read_name not in gt:
                gt[read_name] = genome
    return gt


def analyze(r1_maps, r2_maps, gt):
    """Analyze discordant pairs."""
    total_pairs = 0
    both_mapped = 0
    concordant = 0
    discordant = 0

    discordant_evo = 0
    discordant_other = 0

    # Track specific misclassification patterns
    misclass_patterns = defaultdict(int)

    for read_name in r1_maps:
        if read_name not in r2_maps:
            continue

        total_pairs += 1
        r1_genome = r1_maps[read_name]
        r2_genome = r2_maps[read_name]

        if r1_genome == r2_genome:
            concordant += 1
        else:
            discordant += 1

            # Check ground truth
            if read_name in gt:
                true_genome = gt[read_name]
                r1_correct = (r1_genome == true_genome)
                r2_correct = (r2_genome == true_genome)

                # Classify discordant type
                if 'evo_' in true_genome:
                    discordant_evo += 1
                    if not r1_correct and not r2_correct:
                        misclass_patterns[f"{true_genome}: {r1_genome} vs {r2_genome}"] += 1
                    elif not r1_correct:
                        misclass_patterns[f"{true_genome}: R1={r1_genome} (R2 correct)"] += 1
                    elif not r2_correct:
                        misclass_patterns[f"{true_genome}: R2={r2_genome} (R1 correct)"] += 1
                else:
                    discordant_other += 1

    print("=" * 70)
    print("DISCORDANT PAIR ANALYSIS")
    print("=" * 70)
    print(f"\nTotal pairs in SAM: {total_pairs}")
    print(f"Concordant (R1==R2): {concordant} ({concordant/total_pairs*100:.1f}%)")
    print(f"Discordant (R1!=R2): {discordant} ({discordant/total_pairs*100:.1f}%)")
    print(f"\n  evo_* discordant: {discordant_evo}")
    print(f"  other discordant: {discordant_other}")

    if discordant_evo > 0:
        print(f"\n{'=' * 70}")
        print(f"EVO_* DISCORDANT PATTERNS (top 20)")
        print(f"{'=' * 70}")
        sorted_patterns = sorted(misclass_patterns.items(), key=lambda x: -x[1])
        for pattern, count in sorted_patterns[:20]:
            print(f"  {count:5d}x  {pattern}")

    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"{'=' * 70}")
    print(f"Discordant pairs are mostly evo_*: {discordant_evo/discordant*100:.1f}%")
    print(f"Recommendation: Focus EM/alignment refinement on these {discordant_evo} reads")


if __name__ == '__main__':
    sam_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_mapped.sam'
    gt_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_ground_truth.tsv'

    print("Parsing SAM...")
    r1_maps, r2_maps = parse_sam(sam_path)
    print(f"  R1 mappings: {len(r1_maps)}")
    print(f"  R2 mappings: {len(r2_maps)}")

    print("Parsing ground truth...")
    gt = parse_ground_truth(gt_path)
    print(f"  GT entries: {len(gt)}")

    print("Analyzing...")
    analyze(r1_maps, r2_maps, gt)
