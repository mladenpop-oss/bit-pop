r"""
Extract 20K read pairs from CAMI interleaved FASTQ + generate ground truth TSV.

Usage:
    python extract_cami_20k.py --reads <interleaved_fastq> --gt <ground_truth_tsv> --output <output_fastq> --output-gt <output_ground_truth.tsv> --num-pairs 20000
"""
import argparse
import os
import sys
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(description="Extract 20K CAMI reads + ground truth")
    parser.add_argument("--reads", required=True, help="Path to interleaved FASTQ file")
    parser.add_argument("--gt", required=True, help="Path to ground truth TSV")
    parser.add_argument("--output", required=True, help="Output FASTQ file")
    parser.add_argument("--output-gt", required=True, help="Output ground truth TSV")
    parser.add_argument("--num-pairs", type=int, default=20000, help="Number of read pairs to extract")
    return parser.parse_args()


def normalize_name(header_line):
    """Extract read name from FASTQ header, stripping /1 /2 suffix."""
    name = header_line.lstrip("@").strip()
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def extract_reads(fastq_path, num_pairs):
    """Extract first num_pairs read pairs from interleaved FASTQ."""
    pairs = []
    count = 0

    with open(fastq_path, "r") as f:
        while count < num_pairs:
            header1 = f.readline().strip()
            if not header1:
                break
            seq1 = f.readline().strip()
            plus1 = f.readline().strip()
            qual1 = f.readline().strip()

            header2 = f.readline().strip()
            seq2 = f.readline().strip()
            plus2 = f.readline().strip()
            qual2 = f.readline().strip()

            name1 = normalize_name(header1)
            name2 = normalize_name(header2)

            if name1 == name2:
                pairs.append({
                    "name": name1,
                    "r1_seq": seq1,
                    "r1_qual": qual1,
                    "r2_seq": seq2,
                    "r2_qual": qual2,
                })
                count += 1
            else:
                print(f"WARNING: Mismatched pair at position {count}: {name1} vs {name2}", file=sys.stderr)

    return pairs


def build_ground_truth_lookup(gt_path, read_names):
    """Build lookup: read_name -> (BINID, TAXID) from ground truth TSV."""
    lookup = OrderedDict()
    found = 0
    not_found = 0

    with open(gt_path, "r") as f:
        for line in f:
            if line.startswith("@@SEQUENCEID"):
                continue
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            seq_id = parts[0]
            bin_id = parts[1]
            tax_id = parts[2] if len(parts) > 2 else "?"

            if seq_id in read_names:
                lookup[seq_id] = (bin_id, tax_id)
                found += 1

    for name in read_names:
        if name not in lookup:
            not_found += 1

    return lookup, found, not_found


def write_fastq(pairs, output_path):
    """Write extracted pairs to FASTQ file."""
    with open(output_path, "w") as f:
        for pair in pairs:
            f.write(f"@{pair['name']}/1\n")
            f.write(f"{pair['r1_seq']}\n")
            f.write("+\n")
            f.write(f"{pair['r1_qual']}\n")
            f.write(f"@{pair['name']}/2\n")
            f.write(f"{pair['r2_seq']}\n")
            f.write("+\n")
            f.write(f"{pair['r2_qual']}\n")


def write_ground_truth(lookup, output_path):
    """Write ground truth TSV in CAMI format."""
    with open(output_path, "w") as f:
        f.write("@Version:0.9.1\n")
        f.write("@SampleID:bitpop_20k\n")
        f.write("\n")
        f.write("@@SEQUENCEID\tBINID\tTAXID\t_READID\n")
        for seq_id, (bin_id, tax_id) in lookup.items():
            f.write(f"{seq_id}\t{bin_id}\t{tax_id}\t{seq_id}\n")


def main():
    args = parse_args()

    print(f"Extracting up to {args.num_pairs} read pairs from: {args.reads}")
    pairs = extract_reads(args.reads, args.num_pairs)
    print(f"Extracted {len(pairs)} read pairs")

    if len(pairs) == 0:
        print("ERROR: No reads extracted. Check FASTQ file path.", file=sys.stderr)
        sys.exit(1)

    # Build set of read names for ground truth lookup
    read_names = {p["name"] for p in pairs}

    print(f"Building ground truth lookup from: {args.gt}")
    lookup, found, not_found = build_ground_truth_lookup(args.gt, read_names)
    print(f"Ground truth found: {found}, not found: {not_found}")

    # Write outputs
    write_fastq(pairs, args.output)
    print(f"Wrote FASTQ: {args.output} ({os.path.getsize(args.output) / 1024 / 1024:.1f} MB)")

    write_ground_truth(lookup, args.output_gt)
    print(f"Wrote ground truth: {args.output_gt} ({os.path.getsize(args.output_gt) / 1024:.1f} KB)")

    # Print per-genome distribution
    genome_counts = {}
    for seq_id, (bin_id, _) in lookup.items():
        genome_counts[bin_id] = genome_counts.get(bin_id, 0) + 1

    print(f"\nPer-genome distribution ({len(genome_counts)} genomes):")
    for genome, count in sorted(genome_counts.items(), key=lambda x: -x[1])[:20]:
        print(f"  {genome:30s} {count:6d} pairs")
    if len(genome_counts) > 20:
        print(f"  ... and {len(genome_counts) - 20} more genomes")

    print(f"\nTotal pairs with ground truth: {len(lookup)}")
    print(f"Total pairs WITHOUT ground truth: {len(pairs) - len(lookup)}")


if __name__ == "__main__":
    main()
