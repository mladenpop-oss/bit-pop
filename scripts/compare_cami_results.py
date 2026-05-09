"""
Compare bit-pop SAM mapping results against CAMI ground truth.

Generates a detailed report with:
- Overall mapping rate and accuracy
- Per-genome accuracy
- Per-genome-type accuracy (numeric, Sample*, evo*, virus)
- Confusion matrix (top misclassifications)
- Read name consistency checks

Usage:
    python compare_cami_results.py --sam <output.sam> --gt <ground_truth.tsv> \
        --novelty <novelty_complete.tsv> --unique-common <unique_common.tsv> \
        [--output-report report.txt]
"""
import argparse
import os
import sys
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(description="Compare bit-pop results with CAMI ground truth")
    parser.add_argument("--sam", required=True, help="bit-pop SAM output file")
    parser.add_argument("--gt", required=True, help="Ground truth TSV")
    parser.add_argument("--novelty", default=None, help="novelty_complete.tsv for genome classification")
    parser.add_argument("--unique-common", default=None, help="unique_common.tsv for strain classification")
    parser.add_argument("--output-report", default="cami_test_report.txt", help="Output report file")
    return parser.parse_args()


def load_ground_truth(gt_path):
    """Load ground truth: read_name -> (BINID, TAXID)"""
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


def load_sam_mappings(sam_path):
    """
    Load SAM mappings, handling paired-end reads.
    Returns: read_name -> best_genome (consensus of R1/R2)
    Also returns: read_name -> list of (genome, is_mapped) for debug
    """
    # read_name -> list of (genome_name, is_mapped)
    mappings = OrderedDict()

    # Track which reads we've seen (for paired-end)
    # For paired-end: R1 and R2 have SAME read name (normalized)
    # We take consensus: if both map to same genome -> that genome
    # If they differ -> conflict

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

            if qname not in mappings:
                mappings[qname] = []

            mappings[qname].append({
                "genome": rname if rname != "*" else None,
                "flag": flag,
                "pos": parts[3],
                "mapq": parts[4],
                "cigar": parts[5],
            })

    # Resolve paired-end: take consensus genome per read name
    resolved = OrderedDict()
    conflicts = 0
    both_unmapped = 0

    for read_name, entries in mappings.items():
        genomes = [e["genome"] for e in entries if e["genome"] is not None]
        mapped_count = len(genomes)

        if mapped_count == 0:
            resolved[read_name] = None
            both_unmapped += 1
        elif len(set(genomes)) == 1:
            # Consensus: all mapped entries agree
            resolved[read_name] = genomes[0]
        else:
            # Conflict: R1 and R2 map to different genomes
            resolved[read_name] = genomes[0]  # take first as tiebreak
            conflicts += 1

    return resolved, conflicts, both_unmapped


def classify_genome(genome_name, novelty_path, unique_common_path):
    """Classify genome by type based on naming patterns and metadata files."""
    # Primary: use naming patterns
    if genome_name.startswith("evo_"):
        return "evo_* (similar strains)"
    if genome_name.startswith("Sample"):
        return "Sample* (single-contig)"
    if genome_name.startswith("1") and "_" not in genome_name:
        return "numeric (NCBI ID)"

    # Secondary: use metadata files
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


def compare_results(resolved, gt, novelty_path, unique_common_path):
    """Compare resolved mappings against ground truth."""
    total = 0
    mapped = 0
    correct = 0
    wrong = 0
    unmapped_in_gt = 0

    # Per-genome stats
    genome_stats = {}  # genome -> {total, correct}

    # Per-genome-type stats
    type_stats = {}  # type -> {total, correct}

    # Confusion matrix (top misclassifications)
    confusion = {}  # (true_genome, predicted_genome) -> count

    for read_name, predicted in resolved.items():
        if read_name not in gt:
            continue

        total += 1
        true_genome, tax_id = gt[read_name]

        # Classify true genome
        true_type = classify_genome(true_genome, novelty_path, unique_common_path)
        if true_type not in type_stats:
            type_stats[true_type] = {"total": 0, "correct": 0}
        type_stats[true_type]["total"] += 1

        if true_genome not in genome_stats:
            genome_stats[true_genome] = {"total": 0, "correct": 0}
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
            key = (true_genome, predicted)
            confusion[key] = confusion.get(key, 0) + 1

    return {
        "total": total,
        "mapped": mapped,
        "correct": correct,
        "wrong": wrong,
        "unmapped_in_gt": unmapped_in_gt,
        "genome_stats": genome_stats,
        "type_stats": type_stats,
        "confusion": confusion,
    }


def generate_report(stats, conflicts, both_unmapped, report_path, gt_path, sam_path):
    """Generate comprehensive text report."""
    lines = []

    lines.append("=" * 70)
    lines.append("CAMI TEST BENCHMARK REPORT - bit-pop")
    lines.append("=" * 70)
    lines.append("")

    # Summary
    lines.append("## SUMMARY")
    lines.append("-" * 70)
    lines.append(f"  Total read pairs (with ground truth): {stats['total']}")
    lines.append(f"  Mapped by bit-pop:                     {stats['mapped']}")
    lines.append(f"  Mapping rate:                          {stats['mapped']/stats['total']*100:.2f}%")
    lines.append(f"  Correct mappings:                      {stats['correct']}")
    lines.append(f"  Accuracy (of mapped):                  {stats['correct']/stats['mapped']*100:.2f}%")
    lines.append(f"  Wrong mappings:                        {stats['wrong']}")
    lines.append(f"  Unmapped (but in ground truth):        {stats['unmapped_in_gt']}")
    lines.append(f"  Paired-end conflicts (R1!=R2 genome):  {conflicts}")
    lines.append(f"  Both reads unmapped:                   {both_unmapped}")
    lines.append("")

    # Overall accuracy including unmapped as wrong
    overall_acc = stats['correct'] / stats['total'] * 100 if stats['total'] > 0 else 0
    lines.append(f"  Overall accuracy (correct/total):      {overall_acc:.2f}%")
    lines.append("")

    # Per-genome-type breakdown
    lines.append("## PER-GENOME-TYPE ACCURACY")
    lines.append("-" * 70)
    lines.append(f"  {'Genome Type':<35s} {'Total':>8s} {'Mapped':>8s} {'Correct':>8s} {'Accuracy':>10s}")
    lines.append(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*8} {'-'*10}")

    for type_name in sorted(stats["type_stats"].keys(), key=lambda t: -stats["type_stats"][t]["total"]):
        ts = stats["type_stats"][type_name]
        acc = ts["correct"] / ts["total"] * 100 if ts["total"] > 0 else 0
        lines.append(f"  {type_name:<35s} {ts['total']:>8d} {ts['total']:>8d} {ts['correct']:>8d} {acc:>9.2f}%")

    lines.append("")

    # Per-genome accuracy (sorted by total, top 30)
    lines.append("## PER-GENOME ACCURACY (Top 30 by read count)")
    lines.append("-" * 70)
    lines.append(f"  {'Genome':<35s} {'Total':>8s} {'Correct':>8s} {'Accuracy':>10s}")
    lines.append(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*10}")

    sorted_genomes = sorted(stats["genome_stats"].items(), key=lambda x: -x[1]["total"])
    for genome, gs in sorted_genomes[:30]:
        acc = gs["correct"] / gs["total"] * 100 if gs["total"] > 0 else 0
        lines.append(f"  {genome:<35s} {gs['total']:>8d} {gs['correct']:>8d} {acc:>9.2f}%")

    if len(sorted_genomes) > 30:
        lines.append(f"  ... and {len(sorted_genomes) - 30} more genomes")

    lines.append("")

    # Bottom 10 genomes
    lines.append("## BOTTOM 10 GENOMES BY ACCURACY (min 10 reads)")
    lines.append("-" * 70)
    lines.append(f"  {'Genome':<35s} {'Total':>8s} {'Correct':>8s} {'Accuracy':>10s}")
    lines.append(f"  {'-'*35} {'-'*8} {'-'*8} {'-'*10}")

    min_reads = 10
    low_acc = []
    for genome, gs in stats["genome_stats"].items():
        if gs["total"] >= min_reads:
            acc = gs["correct"] / gs["total"] * 100
            low_acc.append((genome, gs["total"], gs["correct"], acc))

    low_acc.sort(key=lambda x: x[3])
    for genome, total, correct, acc in low_acc[:10]:
        lines.append(f"  {genome:<35s} {total:>8d} {correct:>8d} {acc:>9.2f}%")

    lines.append("")

    # Confusion matrix (top misclassifications)
    lines.append("## TOP 20 MISCLASSIFICATIONS")
    lines.append("-" * 70)
    lines.append(f"  {'True Genome':<25s} {'Predicted As':<25s} {'Count':>8s}")
    lines.append(f"  {'-'*25} {'-'*25} {'-'*8}")

    sorted_confusion = sorted(stats["confusion"].items(), key=lambda x: -x[1])
    for (true_genome, pred_genome), count in sorted_confusion[:20]:
        lines.append(f"  {true_genome:<25s} {pred_genome:<25s} {count:>8d}")

    lines.append("")

    # Files used
    lines.append("## FILES")
    lines.append("-" * 70)
    lines.append(f"  SAM output:       {sam_path}")
    lines.append(f"  Ground truth:     {gt_path}")
    lines.append(f"  SAM size:         {os.path.getsize(sam_path) / 1024 / 1024:.1f} MB")
    lines.append("")

    report_text = "\n".join(lines)

    with open(report_path, "w") as f:
        f.write(report_text)

    print(f"Report written to: {report_path}")
    return report_text


def main():
    args = parse_args()

    # Load data
    print(f"Loading ground truth: {args.gt}")
    gt = load_ground_truth(args.gt)
    print(f"  Loaded {len(gt)} read-genome mappings")

    print(f"Loading SAM mappings: {args.sam}")
    resolved, conflicts, both_unmapped = load_sam_mappings(args.sam)
    print(f"  Loaded {len(resolved)} resolved read mappings")
    print(f"  Paired-end conflicts: {conflicts}")

    # Compare
    print("\nComparing results...")
    stats = compare_results(resolved, gt, args.novelty, args.unique_common)

    # Generate report
    print("\nGenerating report...")
    report = generate_report(stats, conflicts, both_unmapped, args.output_report, args.gt, args.sam)

    # Also print to stdout
    print("\n" + report)


if __name__ == "__main__":
    main()
