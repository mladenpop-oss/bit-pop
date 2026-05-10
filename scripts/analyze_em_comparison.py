#!/usr/bin/env python3
"""Compare bit-pop baseline vs EM results, focusing on evo_* reads."""

import sys
from collections import defaultdict


def parse_sam_genome(sam_path):
    """Parse SAM and extract read_name -> best genome mapping."""
    mappings = {}
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
            
            # Skip unmapped
            if flag & 0x4:
                continue
            
            # Get read name (strip /1 or /2)
            if qname_full.endswith('/1'):
                qname = qname_full[:-2]
            elif qname_full.endswith('/2'):
                qname = qname_full[:-2]
            else:
                qname = qname_full
            
            # Only take primary alignment (first occurrence per read)
            if qname not in mappings:
                mappings[qname] = rname
    
    return mappings


def parse_ground_truth(gt_path):
    """Parse ground truth TSV."""
    gt = {}
    with open(gt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            read_name = parts[0].rstrip('/1').rstrip('/2')
            genome = parts[1]
            if read_name not in gt:
                gt[read_name] = genome
    return gt


def analyze_evo_only(baseline_maps, em_maps, gt):
    """Analyze only evo_* reads."""
    # Categories
    baseline_correct = 0
    baseline_wrong = 0
    
    em_fixed = 0
    em_broke = 0
    em_still_wrong = 0
    em_unchanged_wrong = 0
    
    # Specific misclassifications
    original_errors = defaultdict(int)
    em_errors = defaultdict(int)
    
    for read_name, true_genome in gt.items():
        if 'evo_' not in true_genome:
            continue
        
        if read_name not in baseline_maps or read_name not in em_maps:
            continue
        
        baseline_genome = baseline_maps[read_name]
        em_genome = em_maps[read_name]
        
        baseline_ok = (baseline_genome == true_genome)
        em_ok = (em_genome == true_genome)
        
        if baseline_ok:
            baseline_correct += 1
            if not em_ok:
                em_broke += 1
                em_errors[f"{true_genome}: {em_genome}"] += 1
        else:
            baseline_wrong += 1
            original_errors[f"{true_genome} -> {baseline_genome}"] += 1
            
            if em_ok:
                em_fixed += 1
            elif em_genome == baseline_genome:
                em_unchanged_wrong += 1
            else:
                em_errors[f"{true_genome}: {baseline_genome} -> {em_genome}"] += 1
                if em_ok:
                    em_fixed += 1
    
    print("=" * 70)
    print("EVO_* READS: BASELINE vs EM ANALYSIS")
    print("=" * 70)
    print(f"\nTotal evo_* reads with GT: {len(gt)}")
    print(f"\nBaseline results:")
    print(f"  Correct: {baseline_correct}")
    print(f"  Wrong: {baseline_wrong}")
    if baseline_correct + baseline_wrong > 0:
        print(f"  Accuracy: {baseline_correct/(baseline_correct+baseline_wrong)*100:.1f}%")
    
    print(f"\nEM impact on WRONG baseline predictions ({baseline_wrong} reads):")
    print(f"  Fixed: {em_fixed} ({em_fixed/baseline_wrong*100:.1f}%)")
    print(f"  Broke (was correct): {em_broke}")
    print(f"  Still wrong (same): {em_unchanged_wrong}")
    print(f"  Still wrong (changed): {baseline_wrong - em_fixed - em_unchanged_wrong}")
    
    if em_fixed > 0:
        print(f"\n{'=' * 70}")
        print(f"FIXED MISCLASSIFICATIONS (top 15)")
        print(f"{'=' * 70}")
        sorted_fixed = sorted([(k, v) for k, v in em_errors.items() if '->' in k], 
                             key=lambda x: -x[1])
        for pattern, count in sorted_fixed[:15]:
            print(f"  {count:5d}x  {pattern}")
    
    print(f"\n{'=' * 70}")
    print(f"BROKE CORRECT PREDICTIONS (top 15)")
    print(f"{'=' * 70}")
    broke_patterns = [(k, v) for k, v in em_errors.items() if '->' not in k]
    for pattern, count in sorted(broke_patterns, key=lambda x: -x[1])[:15]:
        print(f"  {count:5d}x  {pattern}")
    
    # Overall evo_* accuracy
    total_evo = baseline_correct + baseline_wrong
    em_correct = baseline_correct + em_fixed - em_broke
    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"{'=' * 70}")
    print(f"Baseline evo_* accuracy: {baseline_correct/total_evo*100:.1f}%")
    print(f"EM evo_* accuracy: {em_correct/total_evo*100:.1f}%")
    print(f"Delta: {(em_correct/total_evo - baseline_correct/total_evo)*100:+.1f}pp")


if __name__ == '__main__':
    import os
    
    baseline_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_mapped.sam'
    em_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_em_mapped.sam'
    gt_path = r'C:\Users\Daddy\Documents\GitHub\CAMI_low\cami_20k_ground_truth.tsv'
    
    if len(sys.argv) > 1:
        baseline_path = sys.argv[1]
    if len(sys.argv) > 2:
        em_path = sys.argv[2]
    if len(sys.argv) > 3:
        gt_path = sys.argv[3]
    
    # Check if EM file exists
    if not os.path.exists(em_path):
        print(f"EM SAM file not found: {em_path}")
        print("Run: bit-pop em -i cami_20k_mapped.sam -o cami_20k_em_mapped.sam")
        print("\nOr specify paths:")
        print("  python analyze_em_comparison.py <baseline.sam> <em.sam> <ground_truth.tsv>")
        exit(1)
    
    print("Parsing baseline SAM...")
    baseline_maps = parse_sam_genome(baseline_path)
    print(f"  Mappings: {len(baseline_maps)}")
    
    print("Parsing EM SAM...")
    em_maps = parse_sam_genome(em_path)
    print(f"  Mappings: {len(em_maps)}")
    
    print("Parsing ground truth...")
    gt = parse_ground_truth(gt_path)
    print(f"  GT entries: {len(gt)}")
    
    print("\nAnalyzing evo_* reads...")
    analyze_evo_only(baseline_maps, em_maps, gt)
