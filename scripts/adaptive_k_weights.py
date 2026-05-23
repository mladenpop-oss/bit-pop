"""
Adaptive k-weighting for multi-k consensus with top-n.
Tests different k-weight strategies on discordant reads.
"""
import argparse
import re
from collections import defaultdict


def adaptive_k_weight(read_len, k, strategy="linear"):
    """Calculate adaptive weight for a k-value based on read length."""
    if strategy == "linear":
        return 1.0 + (k - 10) * (read_len - 100) / 500
    elif strategy == "k_priority":
        return k / 8.0
    elif strategy == "inverse":
        return 8.0 / k
    return 1.0


def load_truth(truth_path):
    """Load ground truth mapping."""
    truth = {}
    if not truth_path:
        return truth
    with open(truth_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                read_id = parts[0]
                genome = parts[1]
                base = re.sub(r'/[12]$', '', read_id)
                truth[base] = genome
                truth[read_id] = genome
    return truth


def compare_voting(sam_path, truth_path, limit=None):
    """Compare uniform vs adaptive voting on consensus SAM with top-n."""
    print("Loading ground truth...")
    truth = load_truth(truth_path)
    print(f"  Loaded {len(truth)} truth entries")

    print("Loading SAM...")
    read_k_votes = defaultdict(list)
    read_lengths = {}

    with open(sam_path, 'r') as f:
        count = 0
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            read_name = parts[0]
            flag = int(parts[1])
            if flag & 0x100:
                continue
            read_len = len(parts[9])
            read_lengths[read_name] = read_len
            tag_str = '\t'.join(parts[11:])
            rk_tags = re.findall(r'(RK\d+):Z:(\S+)', tag_str)
            for rk_k, rk_v in rk_tags:
                k_val = int(rk_k[2:])
                read_k_votes[read_name].append((k_val, rk_v))
            count += 1
            if limit and count >= limit:
                break

    # Analyze discordant reads
    uniform_correct = 0
    uniform_total = 0
    adaptive_correct = 0
    adaptive_total = 0
    discordant_uniform = 0
    discordant_adaptive = 0

    strategies = ["uniform", "linear", "k_priority", "inverse"]
    strat_stats = {s: {'correct': 0, 'total': 0} for s in strategies}

    for read_name, votes in read_k_votes.items():
        if len(votes) < 2:
            continue

        base_name = re.sub(r'/[12]$', '', read_name)
        true_genome = truth.get(base_name) or truth.get(read_name)
        if not true_genome:
            continue

        read_len = read_lengths.get(read_name, 150)
        k_genomes = {}
        for k, genome in votes:
            if k not in k_genomes:
                k_genomes[k] = genome

        is_discordant = len(set(k_genomes.values())) > 1

        for strategy in strategies:
            genome_weights = {}
            for k, genome in k_genomes.items():
                if strategy == "uniform":
                    w = 1.0
                else:
                    w = adaptive_k_weight(read_len, k, strategy)
                genome_weights[genome] = genome_weights.get(genome, 0) + w

            pred = max(genome_weights, key=genome_weights.get)
            strat_stats[strategy]['total'] += 1
            if pred == true_genome:
                strat_stats[strategy]['correct'] += 1
            if is_discordant:
                if pred == true_genome:
                    if strategy == "uniform":
                        discordant_uniform += 1
                    elif strategy == "k_priority":
                        discordant_adaptive += 1

    # Print results
    print(f"\n{'='*60}")
    print(f"{'Strategy':<20} {'Accuracy':>10} {'Correct':>10} {'Total':>8}")
    print(f"{'='*60}")
    for strategy in strategies:
        s = strat_stats[strategy]
        if s['total'] > 0:
            acc = s['correct'] / s['total'] * 100
            print(f"{strategy:<20} {acc:>9.2f}% {s['correct']:>10} {s['total']:>8}")
    print(f"{'='*60}")

    print(f"\nDiscordant reads accuracy:")
    print(f"  Uniform:      {discordant_uniform}")
    print(f"  K-priority:   {discordant_adaptive}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    compare_voting(args.sam, args.truth, args.limit)
