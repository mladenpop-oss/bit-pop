"""
Multi-feature weighted classifier prototype for Bit-Pop SAM output.

Uses all available SAM tags (AS, RK, NM, GM, HF, XS, MQ, MD) instead of
just MAPQ-derived scores. Compares against the existing EM classifier.

Usage:
    python scripts/weighted_classifier.py \
        --sam CAMI_low/cami_20k_mapped.sam \
        --gt CAMI_low/cami_20k_ground_truth.tsv \
        --output-report cami_weighted_report.txt
"""
import argparse
import sys
import os
import math
import re
import time
from collections import OrderedDict, defaultdict
from typing import Optional


def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-feature weighted classifier for Bit-Pop SAM output"
    )
    parser.add_argument("--sam", required=True, help="bit-pop SAM output file")
    parser.add_argument("--gt", required=True, help="Ground truth TSV")
    parser.add_argument("--novelty", default=None, help="novelty_complete.tsv")
    parser.add_argument("--unique-common", default=None, help="unique_common.tsv")
    parser.add_argument("--output-report", default="cami_weighted_report.txt", help="Output report file")
    parser.add_argument("--mode", choices=["weighted", "rf", "em"], default="weighted",
                        help="Classification mode: weighted (default), rf (random forest), em (existing EM)")
    parser.add_argument("--weights", default=None,
                        help="Custom weights as comma-separated: AS,RK,NM,GM,HF,XS (default: auto)")
    parser.add_argument("--temperature", type=float, default=0.1, help="Softmax temperature")
    parser.add_argument("--top-k", type=int, default=10, help="Top-K genomes per read")
    parser.add_argument("--min-score", type=float, default=0.5, help="Minimum alignment score")
    parser.add_argument("--debug", action="store_true", help="Print debug info")
    return parser.parse_args()


def parse_sam_tags(fields):
    """Parse optional SAM tags from a line's fields.

    Handles concatenated tags like:
    NM:i:0\tMD:Z:ACT1G2TAS:f:0.84RK:f:0.001HF:f:0.0
    """
    tags = {}

    # Combine all fields (optional fields already sliced to parts[11:])
    all_optional = " ".join(fields)

    # Parse NM (edit distance) - always integer
    nm_match = re.search(r'NM:i:(\d+)', all_optional)
    if nm_match:
        tags["nm"] = int(nm_match.group(1))

    # Parse AS (alignment score) - float
    as_match = re.search(r'AS:f:([\d.]+)', all_optional)
    if as_match:
        tags["as"] = float(as_match.group(1))

    # Parse RK (k-mer rarity) - float
    rk_match = re.search(r'RK:f:([\d.]+)', all_optional)
    if rk_match:
        tags["rk"] = float(rk_match.group(1))

    # Parse MD (mismatch string) - ends before AS: or RK: or end of string
    md_match = re.search(r'MD:Z:([A-Z0-9]+?)(?:AS:f:|RK:f:|GM:f:|HF:f:|XS:f:|MQ:f:|$)', all_optional)
    if md_match:
        tags["md"] = md_match.group(1)

    # Parse GM (Gaussian insert size) - float
    gm_match = re.search(r'GM:f:([\d.]+)', all_optional)
    if gm_match:
        tags["gm"] = float(gm_match.group(1))

    # Parse HF (homopolymer fingerprint) - float
    hf_match = re.search(r'HF:f:([\d.]+)', all_optional)
    if hf_match:
        tags["hf"] = float(hf_match.group(1))

    # Parse XS (suboptimal score) - float
    xs_match = re.search(r'XS:f:([\d.]+)', all_optional)
    if xs_match:
        tags["xs"] = float(xs_match.group(1))

    # Parse MQ (quality penalty) - float
    mq_match = re.search(r'MQ:f:([\d.]+)', all_optional)
    if mq_match:
        tags["mq"] = float(mq_match.group(1))

    return tags


def compute_md_score(md_str, read_len):
    """Compute alignment score from MD tag."""
    if not md_str:
        return 0.0
    matches = 0
    mismatches = 0
    for char in md_str:
        if char.isdigit():
            matches += int(char)
        elif char.isalpha():
            mismatches += 1
    total = matches + mismatches
    if total == 0:
        return 0.0
    return matches / total


def load_ground_truth(gt_path, max_entries=None):
    """Load ground truth: read_name -> (bin_id, tax_id)"""
    gt = OrderedDict()
    count = 0
    with open(gt_path, "r") as f:
        for line in f:
            if line.startswith("@@SEQUENCEID") or line.startswith("@"):
                continue
            if max_entries and count >= max_entries:
                break
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                seq_id = parts[0]
                bin_id = parts[1]
                tax_id = parts[2] if len(parts) > 2 else "?"
                gt[seq_id] = (bin_id, tax_id)
                count += 1
    return gt


def load_sam_full_features(sam_path, min_score=0.5):
    """
    Load SAM mappings with ALL available features.

    Returns:
        mappings: read_name -> list of {genome, features, flags}
        conflicts: paired-end conflict count
        both_unmapped: count of completely unmapped reads
    """
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
            seq = parts[9] if len(parts) > 9 else ""
            read_len = len(seq)

            # Parse all optional tags
            tags = parse_sam_tags(parts[11:])

            # Compute derived scores
            as_score = tags.get("as", 0.0)
            rk_score = tags.get("rk", 0.0)
            nm_score = 1.0 - (tags.get("nm", 0) / read_len) if read_len > 0 else 0.0
            gm_score = tags.get("gm", 0.0)
            hf_score = tags.get("hf", 0.0)
            xs_score = tags.get("xs", 0.0)
            mq_score = tags.get("mq", 0.0)
            md_score = compute_md_score(tags.get("md", ""), read_len)

            # MAPQ as fallback score
            mapq = int(parts[4]) if parts[4] != "0" else 0
            mapq_score = mapq / 60.0

            # Use AS if available, otherwise MAPQ
            primary_score = as_score if as_score > 0 else mapq_score

            is_supplementary = (flag & 0x800) != 0
            is_reverse = (flag & 0x10) != 0
            is_paired = (flag & 0x1) != 0

            entry = {
                "genome": rname if rname != "*" else None,
                "primary_score": primary_score,
                "as": as_score,
                "rk": rk_score,
                "nm_score": nm_score,
                "gm": gm_score,
                "hf": hf_score,
                "xs": xs_score,
                "md": md_score,
                "mq": mq_score,
                "mapq_score": mapq_score,
                "nm": tags.get("nm", 0),
                "flag": flag,
                "is_supplementary": is_supplementary,
                "is_reverse": is_reverse,
                "is_paired": is_paired,
                "read_len": read_len,
            }

            if qname not in raw_mappings:
                raw_mappings[qname] = []
            raw_mappings[qname].append(entry)

    # Resolve paired-end: collect all genome mappings per read
    resolved = OrderedDict()
    conflicts = 0
    both_unmapped = 0

    for read_name, entries in raw_mappings.items():
        genomes = []
        for e in entries:
            if e["genome"] is not None:
                genomes.append(e)

        if not genomes:
            resolved[read_name] = []
            both_unmapped += 1
        else:
            resolved[read_name] = genomes

        unique_genomes = set(g["genome"] for g in genomes)
        if len(unique_genomes) > 1:
            conflicts += 1

    return resolved, conflicts, both_unmapped


def weighted_vote(genome_entries, active_features=None):
    """
    Compute weighted score for a single genome candidate.

    active_features: list of feature names to use, e.g. ["as"], ["as", "rk"], ["as", "md"]
    """
    if active_features is None:
        active_features = ["as"]

    # Aggregate scores per genome
    genome_scores = {}
    for entry in genome_entries:
        g = entry["genome"]
        if g not in genome_scores:
            genome_scores[g] = []
        genome_scores[g].append(entry)

    # Compute final score per genome
    final_scores = {}
    for g, entries in genome_scores.items():
        # Primary: AS score (average across R1/R2)
        as_scores = [e["as"] if e["as"] > 0 else e["mapq_score"] for e in entries]
        avg_as = sum(as_scores) / len(as_scores)

        score = avg_as  # Base score is always AS

        if "rk" in active_features:
            # RK (k-mer rarity) - normalize to 0-1 range
            rk_scores = [e["rk"] for e in entries]
            avg_rk = sum(rk_scores) / len(rk_scores)
            # RK values are typically 0-0.02, normalize
            norm_rk = min(1.0, avg_rk / 0.02)
            score = score * 0.9 + norm_rk * 0.1  # 10% RK influence

        if "nm" in active_features:
            # NM score (edit distance)
            nm_scores = [e["nm_score"] for e in entries]
            avg_nm = sum(nm_scores) / len(nm_scores)
            score = score * 0.9 + avg_nm * 0.1  # 10% NM influence

        if "md" in active_features:
            # MD score (mismatch detail)
            md_scores = [e["md"] for e in entries]
            avg_md = sum(md_scores) / len(md_scores)
            score = score * 0.9 + avg_md * 0.1  # 10% MD influence

        if "gm" in active_features and entries[0]["gm"] > 0:
            # GM (Gaussian insert size)
            gm_scores = [e["gm"] for e in entries]
            avg_gm = sum(gm_scores) / len(gm_scores)
            score = score * 0.9 + avg_gm * 0.1  # 10% GM influence

        if "hf" in active_features and entries[0]["hf"] > 0:
            # HF (homopolymer fingerprint)
            hf_scores = [e["hf"] for e in entries]
            avg_hf = sum(hf_scores) / len(hf_scores)
            score = score * 0.9 + avg_hf * 0.1  # 10% HF influence

        final_scores[g] = score

    return final_scores


def softmax_normalize(genome_scores, temperature=0.1, top_k=10):
    """Normalize scores to probabilities using softmax."""
    if not genome_scores:
        return {}

    sorted_scores = sorted(genome_scores.items(), key=lambda x: -x[1])
    top_scores = sorted_scores[:top_k]

    if not top_scores:
        return {}

    max_score = max(s for _, s in top_scores)
    adjusted = [(g, (s - max_score) / temperature) for g, s in top_scores]

    exps = [(g, math.exp(min(50, max(-50, a)))) for g, a in adjusted]
    total = sum(e for _, e in exps)

    if total == 0:
        n = len(exps)
        return {g: 1.0 / n for g, _ in exps}

    return {g: e / total for g, e in exps}


def classify_genome(genome_name):
    """Classify genome by type."""
    if genome_name.startswith("evo_"):
        return "evo_* (similar strains)"
    if genome_name.startswith("Sample"):
        return "Sample* (single-contig)"
    if genome_name.startswith("1") and "_" not in genome_name:
        return "numeric (NCBI ID)"
    return "other"


def hard_assign_weighted(mappings, active_features=None):
    """
    Hard assign reads using weighted voting.

    active_features: list of feature names to use
    """
    if active_features is None:
        active_features = ["as"]

    assignments = OrderedDict()

    for read_name, entries in mappings.items():
        if not entries:
            assignments[read_name] = None
            continue

        # Group by genome
        genome_groups = defaultdict(list)
        for e in entries:
            genome_groups[e["genome"]].append(e)

        # Weighted voting per genome
        best_genome = None
        best_score = -1

        for genome, group in genome_groups.items():
            score_dict = weighted_vote(group, active_features)
            score = max(score_dict.values()) if score_dict else 0

            if score > best_score:
                best_score = score
                best_genome = genome

        assignments[read_name] = best_genome

    return assignments


def compare_results(assignments, gt, name=""):
    """Compare assignments against ground truth."""
    total = 0
    mapped = 0
    correct = 0
    wrong = 0
    unmapped_in_gt = 0

    genome_stats = defaultdict(lambda: {"total": 0, "correct": 0})
    type_stats = defaultdict(lambda: {"total": 0, "correct": 0})
    confusion = defaultdict(int)

    for read_name, predicted in assignments.items():
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
        "name": name,
        "total": total,
        "mapped": mapped,
        "correct": correct,
        "wrong": wrong,
        "unmapped_in_gt": unmapped_in_gt,
        "genome_stats": dict(genome_stats),
        "type_stats": dict(type_stats),
        "confusion": dict(confusion),
    }


def generate_report(results, report_path):
    """Generate comparison report."""
    lines = []

    lines.append("=" * 70)
    lines.append("MULTI-FEATURE WEIGHTED CLASSIFIER REPORT - bit-pop")
    lines.append("=" * 70)
    lines.append("")

    def pct(n, d):
        return f"{n/d*100:.2f}%" if d > 0 else "0.00%"

    # Summary table
    lines.append("## SUMMARY COMPARISON")
    lines.append("-" * 70)
    header = f"  {'Metric':<35s}" + "".join(f"{r['name']:>15s}" for r in results)
    lines.append(header)
    lines.append(f"  {'-'*35}" + "".join(f"{'-'*15}" for _ in results))

    for metric_name, key in [("Total reads", "total"), ("Mapped", "mapped"), ("Correct", "correct")]:
        line = f"  {metric_name:<35s}" + "".join(f"{r[key]:>15d}" for r in results)
        lines.append(line)

    line = f"  {'Accuracy (of mapped)':<35s}" + "".join(f"{r['correct']/r['mapped']*100 if r['mapped'] > 0 else 0:>14.2f}%" for r in results)
    lines.append(line)

    line = f"  {'Accuracy (of total)':<35s}" + "".join(f"{r['correct']/r['total']*100 if r['total'] > 0 else 0:>14.2f}%" for r in results)
    lines.append(line)
    lines.append("")

    # Per-genome-type comparison
    lines.append("## PER-GENOME-TYPE COMPARISON")
    lines.append("-" * 70)

    all_types = set()
    for r in results:
        all_types.update(r['type_stats'].keys())

    header = f"  {'Genome Type':<30s}" + "".join(f"{r['name']:>15s}" for r in results)
    lines.append(header)
    lines.append(f"  {'-'*30}" + "".join(f"{'-'*15}" for _ in results))

    for type_name in sorted(all_types, key=lambda t: -results[0]['type_stats'].get(t, {}).get('total', 0)):
        line = f"  {type_name:<30s}"
        for r in results:
            ts = r['type_stats'].get(type_name, {"total": 0, "correct": 0})
            acc = ts["correct"] / ts["total"] * 100 if ts["total"] > 0 else 0
            line += f"{acc:>14.2f}%"
        lines.append(line)
    lines.append("")

    # Evo_* analysis
    lines.append("## EVO_* STRAIN ANALYSIS")
    lines.append("-" * 70)
    for r in results:
        evo = r['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
        acc = evo["correct"] / evo["total"] * 100 if evo["total"] > 0 else 0
        lines.append(f"  {r['name']:<20s}: {acc:.2f}% ({evo['correct']}/{evo['total']})")
    lines.append("")

    # Top wrong predictions
    lines.append("## TOP WRONG PREDICTIONS (weighted)")
    lines.append("-" * 70)
    if results:
        weighted_result = results[0]
        confusion = weighted_result.get("confusion", {})
        top_confusions = sorted(confusion.items(), key=lambda x: -x[1])[:10]
        for (true_genome, pred_genome), count in top_confusions:
            lines.append(f"  {true_genome} -> {pred_genome}: {count} reads")
    lines.append("")

    lines.append(f"## FILES")
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
    print("Multi-Feature Weighted Classifier - bit-pop")
    print("=" * 70)
    print()

    # Load ground truth
    print(f"Loading ground truth: {args.gt}")
    gt = load_ground_truth(args.gt, max_entries=500000)
    print(f"  Loaded {len(gt)} read-genome mappings")

    # Load SAM with full features
    print(f"Loading SAM mappings: {args.sam}")
    mappings, sam_conflicts, both_unmapped = load_sam_full_features(
        args.sam, min_score=args.min_score
    )
    print(f"  Loaded {len(mappings)} resolved read mappings")
    print(f"  Paired-end conflicts: {sam_conflicts}")
    print()

    # Show feature availability stats
    print("## Feature availability across all mappings:")
    feature_sums = defaultdict(float)
    feature_counts = defaultdict(int)
    for entries in mappings.values():
        for e in entries:
            if e["genome"] is not None:
                if e["as"] > 0:
                    feature_counts["as"] += 1
                    feature_sums["as"] += e["as"]
                if e["rk"] > 0:
                    feature_counts["rk"] += 1
                    feature_sums["rk"] += e["rk"]
                if e["gm"] > 0:
                    feature_counts["gm"] += 1
                    feature_sums["gm"] += e["gm"]
                if e["hf"] > 0:
                    feature_counts["hf"] += 1
                    feature_sums["hf"] += e["hf"]
                if e["md"] > 0:
                    feature_counts["md"] += 1
                    feature_sums["md"] += e["md"]
    total_entries = sum(1 for entries in mappings.values() for e in entries if e["genome"])
    for feature in ["as", "rk", "gm", "hf", "md"]:
        count = feature_counts[feature]
        avg = feature_sums[feature] / count if count > 0 else 0
        pct_val = count / total_entries * 100 if total_entries > 0 else 0
        print(f"  {feature:>3s}: {count:>6d} ({pct_val:.1f}%) avg={avg:.4f}")
    print()

    # Systematic feature testing
    print("## SYSTEMATIC FEATURE TESTING")
    print("Testing each feature individually and in combination...")
    print()

    # Define feature combinations to test (start with key ones)
    feature_tests = [
        ["as"],           # Baseline: AS only
        ["as", "rk"],     # AS + k-mer rarity
        ["as", "nm"],     # AS + edit distance
        ["as", "md"],     # AS + mismatch detail
    ]

    results = []
    print(f"  {'Features':<30s} {'Accuracy':>10s} {'Evo_*':>10s} {'Delta':>10s}")
    print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10}")

    base_result = None
    for features in feature_tests:
        start_time = time.time()
        assignments = hard_assign_weighted(mappings, active_features=features)
        elapsed = time.time() - start_time

        result = compare_results(assignments, gt, name="+".join(features))
        acc = result['correct'] / result['mapped'] * 100 if result['mapped'] > 0 else 0
        evo = result['type_stats'].get("evo_* (similar strains)", {"total": 0, "correct": 0})
        evo_acc = evo["correct"] / evo["total"] * 100 if evo["total"] > 0 else 0

        delta = ""
        if base_result:
            base_acc = base_result['correct'] / base_result['mapped'] * 100 if base_result['mapped'] > 0 else 0
            delta = f"{acc - base_acc:+.2f}pp"

        if features == ["as"]:
            base_result = result

        features_str = " + ".join(f.upper() for f in features)
        print(f"  {features_str:<30s} {acc:>9.2f}% {evo_acc:>9.2f}% {delta:>10s}")

        results.append((features, result, elapsed))

    # Generate report with all results
    print("\nGenerating report...")
    report = generate_report([r[1] for r in results], args.output_report)

    print("\n" + report)

    # Key findings
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    base_acc = base_result['correct'] / base_result['mapped'] * 100 if base_result['mapped'] > 0 else 0
    best_result = max(results, key=lambda x: x[1]['correct'] / x[1]['mapped'] * 100 if x[1]['mapped'] > 0 else 0)
    best_acc = best_result[1]['correct'] / best_result[1]['mapped'] * 100 if best_result[1]['mapped'] > 0 else 0
    best_features = best_result[0]

    print(f"  Baseline (AS only):     {base_acc:.2f}%")
    print(f"  Best combination:       {best_acc:.2f}% ({'+'.join(f.upper() for f in best_features)})")
    print(f"  Improvement:            {best_acc - base_acc:+.2f}pp")


if __name__ == "__main__":
    main()
