"""
EM Algorithm for soft-assignment read classification.

Improves strain-level resolution for highly similar genomes (evo_* strains)
by using probabilistic assignments instead of hard classification.

Algorithm:
  Phase 1: Soft assignment - each read gets probabilities for ALL genomes
  Phase 2: EM iterations - E-step (compute probabilities) + M-step (update abundances)
  Phase 3: Back-propagate - hard assignment from final probabilities

Usage:
    python scripts/em_classifier.py \
        --sam CAMI_low/cami_20k_mapped.sam \
        --gt CAMI_low/cami_20k_ground_truth.tsv \
        --novelty CAMI_low/novelty_complete.tsv \
        --unique-common CAMI_low/unique_common.tsv \
        --output-report cami_em_report.txt \
        --max-iterations 50 \
        --convergence-threshold 0.001

This is a POST-PROCESSING tool that reads existing bit-pop SAM output
and improves classifications using EM algorithm.
"""
import argparse
import sys
import os
import math
import time
from collections import OrderedDict, defaultdict
from typing import Optional


def parse_args():
    parser = argparse.ArgumentParser(
        description="EM Algorithm for soft-assignment read classification"
    )
    parser.add_argument("--sam", required=True, help="bit-pop SAM output file")
    parser.add_argument("--gt", required=True, help="Ground truth TSV")
    parser.add_argument("--novelty", default=None, help="novelty_complete.tsv for genome classification")
    parser.add_argument("--unique-common", default=None, help="unique_common.tsv for strain classification")
    parser.add_argument("--output-report", default="cami_em_report.txt", help="Output report file")
    parser.add_argument("--max-iterations", type=int, default=50, help="Max EM iterations")
    parser.add_argument("--convergence-threshold", type=float, default=0.001, help="EM convergence threshold")
    parser.add_argument("--min-abundance", type=float, default=1e-6, help="Minimum abundance to keep a genome")
    parser.add_argument("--temperature", type=float, default=0.1, help="Temperature for score normalization (lower = sharper)")
    parser.add_argument("--min-score", type=float, default=0.5, help="Minimum alignment score to consider a mapping")
    parser.add_argument("--top-k", type=int, default=10, help="Top-K genomes per read for EM")
    parser.add_argument("--debug", action="store_true", help="Print debug info during EM iterations")
    return parser.parse_args()


def load_ground_truth(gt_path):
    """Load ground truth: read_name -> (bin_id, tax_id)"""
    gt = OrderedDict()
    with open(gt_path, "r") as f:
        for line in f:
            if line.startswith("@@SEQUENCEID") or line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                seq_id = parts[0]
                bin_id = parts[1]
                tax_id = parts[2] if len(parts) > 2 else "?"
                gt[seq_id] = (bin_id, tax_id)
    return gt


def load_sam_with_scores(sam_path, min_score=0.5):
    """
    Load SAM mappings with scores.
    Returns:
        mappings: read_name -> list of (genome_name, score, is_mapped)
        conflicts: number of paired-end conflicts
        both_unmapped: number of reads with both reads unmapped
    """
    # read_name -> list of (genome_name, score, is_mapped)
    raw_mappings = OrderedDict()

    with open(sam_path, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue

            qname = parts[0].rstrip("/1").rstrip("/2")
            flag = int(parts[1])
            rname = parts[2]

            # Extract MAPQ as a proxy for score
            mapq = int(parts[4]) if parts[4] != "0" else 0
            score = mapq / 60.0  # Normalize MAPQ to 0-1 range

            # Extract NM tag from optional fields
            nm = 0
            for field in parts[11:]:
                if field.startswith("NM:i:"):
                    try:
                        nm = int(field.split(":")[2])
                    except (ValueError, IndexError):
                        pass
                    break

            # Weight by flag: primary (0/16) gets full score, supplementary (2048) gets reduced
            is_supplementary = (flag & 0x800) != 0
            if is_supplementary:
                score *= 0.5  # Supplementary alignments get half score

            # NM-based scoring: lower NM = better match
            # Use read length from SEQ field to compute mismatch rate
            seq = parts[9] if len(parts) > 9 else ""
            read_len = len(seq)
            if read_len > 0 and nm > 0:
                mismatch_rate = nm / read_len
                nm_score = 1.0 - mismatch_rate  # 1.0 = perfect, 0.0 = all mismatches
                # Blend MAPQ score and NM score (70% MAPQ, 30% NM)
                score = score * 0.7 + nm_score * 0.3

            if qname not in raw_mappings:
                raw_mappings[qname] = []

            raw_mappings[qname].append({
                "genome": rname if rname != "*" else None,
                "score": score,
                "nm": nm,
                "flag": flag,
                "pos": parts[3],
                "cigar": parts[5],
            })

    # Resolve paired-end: collect all genome mappings per read
    resolved = OrderedDict()
    conflicts = 0
    both_unmapped = 0

    for read_name, entries in raw_mappings.items():
        genomes = []
        for e in entries:
            if e["genome"] is not None:
                genomes.append((e["genome"], e["score"]))

        if not genomes:
            resolved[read_name] = []
            both_unmapped += 1
        else:
            resolved[read_name] = genomes

        # Check for paired-end conflict
        unique_genomes = set(g[0] for g in genomes)
        if len(unique_genomes) > 1:
            conflicts += 1

    return resolved, conflicts, both_unmapped


def classify_genome(genome_name, novelty_path=None, unique_common_path=None):
    """Classify genome by type based on naming patterns."""
    if genome_name.startswith("evo_"):
        return "evo_* (similar strains)"
    if genome_name.startswith("Sample"):
        return "Sample* (single-contig)"
    if genome_name.startswith("1") and "_" not in genome_name:
        return "numeric (NCBI ID)"

    if novelty_path and os.path.exists(novelty_path):
        with open(novelty_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 1 and parts[-1] == genome_name:
                    novelty_class = parts[0] if len(parts) > 1 else "unknown"
                    if "virus" in genome_name.lower() or "virus" in novelty_class.lower():
                        return "virus"
                    return novelty_class

    if unique_common_path and os.path.exists(unique_common_path):
        with open(unique_common_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[-1] == genome_name:
                    return parts[0] if len(parts) > 0 else "unknown"

    return "other"


def normalize_to_probabilities(genome_scores, temperature=0.1, top_k=10):
    """
    Convert genome scores to probability distribution using softmax with temperature.

    Args:
        genome_scores: list of (genome_name, score) tuples
        temperature: softmax temperature (lower = sharper distribution)
        top_k: only keep top-K genomes

    Returns:
        dict: genome_name -> probability
    """
    if not genome_scores:
        return {}

    # Sort by score descending
    sorted_scores = sorted(genome_scores, key=lambda x: -x[1])

    # Take top-K
    top_scores = sorted_scores[:top_k]

    # Filter by minimum score
    top_scores = [(g, s) for g, s in top_scores if s >= 0.0]

    if not top_scores:
        return {}

    # Softmax with temperature
    logit_scores = [s for _, s in top_scores]

    # Numerical stability: subtract max
    max_logit = max(logit_scores) if logit_scores else 0
    adjusted = [(g, (s - max_logit) / temperature) for g, s in top_scores]

    # Exp and normalize
    exps = [(g, math.exp(min(50, max(-50, a)))) for g, a in adjusted]
    total = sum(e for _, e in exps)

    if total == 0:
        # Uniform distribution as fallback
        n = len(exps)
        return {g: 1.0 / n for g, _ in exps}

    return {g: e / total for g, e in exps}


class EMClassifier:
    """
    EM Algorithm for soft-assignment read classification.

    Works on top of existing bit-pop hard classifications.
    Improves strain-level resolution by using probabilistic assignments.
    """

    def __init__(self, min_abundance=1e-6, max_iterations=50, convergence_threshold=0.001):
        self.min_abundance = min_abundance
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold

        # All genomes seen across all reads
        self.all_genomes = set()

        # read_id -> {genome_id: probability}
        self.soft_assignments = OrderedDict()

        # Genome abundances (theta)
        self.abundances = {}

        # Statistics
        self.iteration_history = []

    def initialize(self, mappings):
        """
        Initialize soft assignments from hard classifications.

        Args:
            mappings: read_name -> list of (genome_name, score)
        """
        # Collect all genomes
        for read_name, genome_scores in mappings.items():
            for genome_name, score in genome_scores:
                self.all_genomes.add(genome_name)

        # Initialize soft assignments using softmax
        for read_name, genome_scores in mappings.items():
            probs = normalize_to_probabilities(genome_scores, temperature=0.1, top_k=10)
            if probs:
                self.soft_assignments[read_name] = probs

        # Initialize uniform abundances
        n_genomes = len(self.all_genomes)
        if n_genomes > 0:
            self.abundances = {g: 1.0 / n_genomes for g in self.all_genomes}

    def e_step(self):
        """
        E-step: Compute P(genome_i | read) = P(read | genome_i) * theta_i / Z

        Uses current abundance estimates as priors.
        """
        for read_name, genome_scores in self.soft_assignments.items():
            if not genome_scores:
                continue

            # Apply abundance priors
            weighted_probs = {}
            for genome, prob in genome_scores.items():
                abundance = self.abundances.get(genome, self.min_abundance)
                weighted_probs[genome] = prob * abundance

            # Normalize
            total = sum(weighted_probs.values())
            if total > 0:
                self.soft_assignments[read_name] = {
                    g: p / total for g, p in weighted_probs.items()
                }
            else:
                # Fallback to uniform
                n = len(genome_scores)
                self.soft_assignments[read_name] = {g: 1.0 / n for g in genome_scores}

    def m_step(self):
        """
        M-step: Update theta_i = sum(P(genome_i | read_r)) / total_reads

        Recomputes genome abundances from soft assignments.
        """
        # Sum probabilities for each genome
        sums = defaultdict(float)
        total_reads = 0

        for probs in self.soft_assignments.values():
            for genome, prob in probs.items():
                sums[genome] += prob
            total_reads += 1

        if total_reads == 0:
            return

        # Normalize to get new abundances
        new_abundances = {}
        for genome in self.all_genomes:
            abundance = sums.get(genome, 0) / total_reads
            if abundance >= self.min_abundance:
                new_abundances[genome] = abundance

        # Renormalize
        total = sum(new_abundances.values())
        if total > 0:
            self.abundances = {g: a / total for g, a in new_abundances.items()}

    def run(self, mappings, debug=False):
        """
        Run EM algorithm.

        Args:
            mappings: read_name -> list of (genome_name, score)
            debug: print iteration info

        Returns:
            dict: iteration statistics
        """
        self.initialize(mappings)
        n_reads = len(self.soft_assignments)

        if debug:
            print(f"EM initialized: {n_reads} reads, {len(self.all_genomes)} genomes")
            print(f"Initial uniform abundance: {1.0/len(self.all_genomes):.6f}")

        prev_abundances = None

        for iteration in range(1, self.max_iterations + 1):
            # E-step
            self.e_step()

            # M-step
            self.m_step()

            # Check convergence (KL divergence between iterations)
            if prev_abundances is not None:
                kl_divergence = self._compute_kl_divergence(prev_abundances, self.abundances)
            else:
                kl_divergence = float('inf')

            # Record statistics
            max_abundance = max(self.abundances.values()) if self.abundances else 0
            min_abundance = min(self.abundances.values()) if self.abundances else 0
            active_genomes = sum(1 for a in self.abundances.values() if a > self.min_abundance)

            stats = {
                "iteration": iteration,
                "kl_divergence": kl_divergence,
                "max_abundance": max_abundance,
                "min_abundance": min_abundance,
                "active_genomes": active_genomes,
            }
            self.iteration_history.append(stats)

            if debug:
                print(f"  Iter {iteration:3d}: KL={kl_divergence:.6f}, "
                      f"active={active_genomes}, max_abund={max_abundance:.4f}")

            # Check convergence
            if kl_divergence < self.convergence_threshold:
                if debug:
                    print(f"  Converged at iteration {iteration} (KL={kl_divergence:.6f})")
                break

            prev_abundances = dict(self.abundances)

        return {
            "iterations": len(self.iteration_history),
            "final_kl": self.iteration_history[-1]["kl_divergence"] if self.iteration_history else 0,
            "active_genomes": self.iteration_history[-1]["active_genomes"] if self.iteration_history else 0,
        }

    def _compute_kl_divergence(self, p, q):
        """Compute KL divergence: KL(p || q) = sum(p * log(p/q))"""
        kl = 0.0
        for genome in self.all_genomes:
            p_val = p.get(genome, self.min_abundance)
            q_val = q.get(genome, self.min_abundance)
            if p_val > 0:
                kl += p_val * math.log(p_val / q_val)
        return kl

    def back_propagate(self):
        """
        Convert soft assignments back to hard classifications.

        Returns:
            dict: read_name -> best_genome (or None if no assignment)
        """
        hard_assignments = OrderedDict()

        for read_name, probs in self.soft_assignments.items():
            if not probs:
                hard_assignments[read_name] = None
            else:
                best_genome = max(probs, key=probs.get)
                hard_assignments[read_name] = best_genome

        return hard_assignments

    def get_abundance_report(self):
        """Get sorted abundance report."""
        sorted_genomes = sorted(
            self.abundances.items(),
            key=lambda x: -x[1]
        )
        return sorted_genomes


def compare_results_with_em(hard_assignments, gt):
    """Compare EM hard assignments against ground truth."""
    total = 0
    mapped = 0
    correct = 0
    wrong = 0
    unmapped_in_gt = 0

    genome_stats = defaultdict(lambda: {"total": 0, "correct": 0})
    type_stats = defaultdict(lambda: {"total": 0, "correct": 0})
    confusion = defaultdict(int)

    for read_name, predicted in hard_assignments.items():
        if read_name not in gt:
            continue

        total += 1
        true_genome, tax_id = gt[read_name]
        true_type = classify_genome(true_genome)

        type_stats[true_type]["total"] += 1
        genome_stats[true_genome]["total"] += 1

        if predicted is None:
            unmapped_in_gt += 1
            continue

        mapped += 1

        if predicted == true_genome:
            correct += 1
            type_stats[true_type]["correct"] += 1
            genome_stats[true_genome]["correct"] += 1
        else:
            wrong += 1
            confusion[(true_genome, predicted)] += 1

    return {
        "total": total,
        "mapped": mapped,
        "correct": correct,
        "wrong": wrong,
        "unmapped_in_gt": unmapped_in_gt,
        "genome_stats": dict(genome_stats),
        "type_stats": dict(type_stats),
        "confusion": dict(confusion),
    }


def generate_report(em_stats, original_stats, em_conflicts, original_conflicts,
                    abundance_report, report_path):
    """Generate comprehensive comparison report."""
    lines = []

    lines.append("=" * 70)
    lines.append("CAMI EM CLASSIFIER REPORT - bit-pop")
    lines.append("=" * 70)
    lines.append("")

    # Summary comparison
    lines.append("## SUMMARY COMPARISON (EM vs Original)")
    lines.append("-" * 70)
    lines.append(f"  {'Metric':<35s} {'Original':>12s} {'EM':>12s} {'Delta':>12s}")
    lines.append(f"  {'-'*35} {'-'*12} {'-'*12} {'-'*12}")

    def pct(n, d):
        return f"{n/d*100:.2f}%" if d > 0 else "0.00%"

    # Overall accuracy
    orig_acc = original_stats['correct'] / original_stats['total'] * 100 if original_stats['total'] > 0 else 0
    em_acc = em_stats['correct'] / em_stats['total'] * 100 if em_stats['total'] > 0 else 0

    metrics = [
        ("Total reads (with GT)", original_stats['total'], em_stats['total']),
        ("Mapped", original_stats['mapped'], em_stats['mapped']),
        ("Mapping rate", pct(original_stats['mapped'], original_stats['total']),
         pct(em_stats['mapped'], em_stats['total'])),
        ("Correct", original_stats['correct'], em_stats['correct']),
        ("Accuracy (of mapped)", pct(original_stats['correct'], original_stats['mapped']),
         pct(em_stats['correct'], em_stats['mapped'])),
        ("Accuracy (of total)", f"{orig_acc:.2f}%", f"{em_acc:.2f}%"),
        ("Paired-end conflicts", original_conflicts, em_conflicts),
    ]

    for name, orig_val, em_val in metrics:
        if isinstance(orig_val, str):
            delta = f"+{float(em_val.rstrip('%')) - float(orig_val.rstrip('%')):.2f}pp" if '%' in orig_val else "?"
            lines.append(f"  {name:<35s} {orig_val:>12s} {em_val:>12s} {delta:>12s}")
        else:
            delta = em_val - orig_val
            delta_str = f"+{delta}" if delta >= 0 else str(delta)
            lines.append(f"  {name:<35s} {orig_val:>12d} {em_val:>12d} {delta_str:>12s}")

    lines.append("")

    # Per-genome-type comparison
    lines.append("## PER-GENOME-TYPE COMPARISON")
    lines.append("-" * 70)
    lines.append(f"  {'Genome Type':<30s} {'Orig Acc':>10s} {'EM Acc':>10s} {'Delta':>10s}")
    lines.append(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10}")

    all_types = set(list(original_stats['type_stats'].keys()) + list(em_stats['type_stats'].keys()))
    for type_name in sorted(all_types, key=lambda t: -original_stats['type_stats'].get(t, {}).get('total', 0)):
        orig_ts = original_stats['type_stats'].get(type_name, {"total": 0, "correct": 0})
        em_ts = em_stats['type_stats'].get(type_name, {"total": 0, "correct": 0})
        orig_acc = orig_ts["correct"] / orig_ts["total"] * 100 if orig_ts["total"] > 0 else 0
        em_acc = em_ts["correct"] / em_ts["total"] * 100 if em_ts["total"] > 0 else 0
        delta = em_acc - orig_acc
        delta_str = f"+{delta:.2f}pp" if delta >= 0 else f"{delta:.2f}pp"
        lines.append(f"  {type_name:<30s} {orig_acc:>9.2f}% {em_acc:>9.2f}% {delta_str:>10s}")

    lines.append("")

    # Evo_* detailed analysis
    evo_orig = original_stats['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
    evo_em = em_stats['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
    evo_orig_acc = evo_orig["correct"] / evo_orig["total"] * 100 if evo_orig["total"] > 0 else 0
    evo_em_acc = evo_em["correct"] / evo_em["total"] * 100 if evo_em["total"] > 0 else 0

    lines.append("## EVO_* STRAIN ANALYSIS")
    lines.append("-" * 70)
    lines.append(f"  Original accuracy: {evo_orig_acc:.2f}% ({evo_orig['correct']}/{evo_orig['total']})")
    lines.append(f"  EM accuracy:       {evo_em_acc:.2f}% ({evo_em['correct']}/{evo_em['total']})")
    lines.append(f"  Improvement:       {evo_em_acc - evo_orig_acc:+.2f}pp")
    lines.append("")

    # Abundance report
    lines.append("## GENOME ABUNDANCES (EM result)")
    lines.append("-" * 70)
    lines.append(f"  {'Rank':>4s} {'Genome':<30s} {'Abundance':>12s} {'Type':<25s}")
    lines.append(f"  {'-'*4} {'-'*30} {'-'*12} {'-'*25}")

    for rank, (genome, abundance) in enumerate(abundance_report[:20], 1):
        genome_type = classify_genome(genome)
        lines.append(f"  {rank:>4d} {genome:<30s} {abundance:>11.4f} {genome_type:<25s}")

    if len(abundance_report) > 20:
        lines.append(f"  ... and {len(abundance_report) - 20} more genomes")

    lines.append("")

    # EM iteration history
    if em_stats.get("iterations", 0) > 0:
        lines.append("## EM CONVERGENCE")
        lines.append("-" * 70)
        lines.append(f"  Iterations: {em_stats['iterations']}")
        lines.append(f"  Final KL divergence: {em_stats['final_kl']:.6f}")
        lines.append(f"  Active genomes at end: {em_stats['active_genomes']}")
        lines.append("")

    # Files
    lines.append("## FILES")
    lines.append("-" * 70)
    lines.append(f"  Report: {report_path}")
    lines.append("")

    report_text = "\n".join(lines)

    with open(report_path, "w") as f:
        f.write(report_text)

    print(f"Report written to: {report_path}")
    return report_text


def main():
    args = parse_args()

    print("=" * 70)
    print("CAMI EM Classifier - Soft Assignment for Strain Resolution")
    print("=" * 70)
    print()

    # Load ground truth
    print(f"Loading ground truth: {args.gt}")
    gt = load_ground_truth(args.gt)
    print(f"  Loaded {len(gt)} read-genome mappings")

    # Load SAM with scores
    print(f"Loading SAM mappings: {args.sam}")
    mappings, sam_conflicts, both_unmapped = load_sam_with_scores(
        args.sam, min_score=args.min_score
    )
    print(f"  Loaded {len(mappings)} resolved read mappings")
    print(f"  Paired-end conflicts: {sam_conflicts}")

    # Run EM
    print(f"\nRunning EM algorithm...")
    print(f"  Max iterations: {args.max_iterations}")
    print(f"  Convergence threshold: {args.convergence_threshold}")
    print(f"  Temperature: {args.temperature}")
    print(f"  Top-K genomes per read: {args.top_k}")
    print()

    em = EMClassifier(
        min_abundance=args.min_abundance,
        max_iterations=args.max_iterations,
        convergence_threshold=args.convergence_threshold,
    )

    start_time = time.time()
    em_stats = em.run(mappings, debug=args.debug)
    elapsed = time.time() - start_time

    print(f"\nEM completed in {elapsed:.1f}s")
    print(f"  Iterations: {em_stats['iterations']}")
    print(f"  Active genomes: {em_stats['active_genomes']}")
    print()

    # Back-propagate to hard assignments
    print("Back-propagating soft assignments to hard classifications...")
    em_hard_assignments = em.back_propagate()

    # Compare results
    print("\nComparing EM results against ground truth...")
    em_comparison = compare_results_with_em(em_hard_assignments, gt)

    # Load original results for comparison
    print("Loading original hard classifications for comparison...")

    def load_original_hard_assignments(mappings, min_score=0.5):
        """Load original hard assignments (best score per read)."""
        hard = OrderedDict()
        for read_name, genome_scores in mappings.items():
            if not genome_scores:
                hard[read_name] = None
            else:
                best = max(genome_scores, key=lambda x: x[1])
                hard[read_name] = best[0]
        return hard

    original_hard = load_original_hard_assignments(mappings, args.min_score)
    original_comparison = compare_results_with_em(original_hard, gt)

    # Get abundance report
    abundance_report = em.get_abundance_report()

    # Generate report
    print("\nGenerating report...")
    report = generate_report(
        em_comparison, original_comparison,
        sam_conflicts, sam_conflicts,
        abundance_report, args.output_report
    )

    print("\n" + report)

    # Print key findings
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    evo_orig = original_comparison['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
    evo_em = em_comparison['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
    evo_orig_acc = evo_orig["correct"] / evo_orig["total"] * 100 if evo_orig["total"] > 0 else 0
    evo_em_acc = evo_em["correct"] / evo_em["total"] * 100 if evo_em["total"] > 0 else 0

    orig_acc = original_comparison['correct'] / original_comparison['total'] * 100 if original_comparison['total'] > 0 else 0
    em_acc = em_comparison['correct'] / em_comparison['total'] * 100 if em_comparison['total'] > 0 else 0

    print(f"  Overall accuracy: {orig_acc:.2f}% -> {em_acc:.2f}% ({em_acc - orig_acc:+.2f}pp)")
    print(f"  Evo_* accuracy:   {evo_orig_acc:.2f}% -> {evo_em_acc:.2f}% ({evo_em_acc - evo_orig_acc:+.2f}pp)")

    if evo_em_acc - evo_orig_acc > 5:
        print(f"\n  SUCCESS: EM improved evo_* accuracy by {evo_em_acc - evo_orig_acc:.2f}pp!")
    elif evo_em_acc - evo_orig_acc > 0:
        print(f"\n  MODERATE: EM improved evo_* accuracy by {evo_em_acc - evo_orig_acc:.2f}pp")
    else:
        print(f"\n  NO IMPROVEMENT: EM did not improve evo_* accuracy")
        print(f"  Possible reasons:")
        print(f"    - Temperature too high/low")
        print(f"    - Not enough EM iterations")
        print(f"    - SAM MAPQ scores don't reflect true confidence")
        print(f"    - Need to use k=12 index for more strain-specific k-mers")


if __name__ == "__main__":
    main()
