# Bit-Pop: Multi-Genome DNA Read Classification

[![CI](https://github.com/mladenpop-oss/bit-pop/actions/workflows/ci.yml/badge.svg)](https://github.com/mladenpop-oss/bit-pop/actions/workflows/ci.yml)
[![Docker](https://github.com/mladenpop-oss/bit-pop/actions/workflows/docker.yml/badge.svg)](https://github.com/mladenpop-oss/bit-pop/pkgs/container/bit-pop)
[![Tests](https://img.shields.io/badge/tests-312%2B%20unit%2C%205%20integration-blue)](https://github.com/mladenpop-oss/bit-pop)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043593.svg)](https://doi.org/10.5281/zenodo.20043593)

> **Ultra-fast multi-genome DNA read classification.** Maps reads across dozens of genomes using a compact FM-index built in Rust with bit-level parallelism.

**Quick benchmark** (3 genomes, 19.7 Mb): **99.3% mapping rate**, **99.9% classification accuracy**, **0.9s per 10k reads**.

**PacBio HiFi benchmark** (69 genomes, 285 Mb, 86k long reads): **95.2% accuracy**, **99% mapping rate**, **8 min** (realistic error profile: homopolymers, chimera, coverage variation).

**CAMI benchmark** (61 genomes, 157 Mb, ~1M reads):
- **Best accuracy**: k13 tn=4 + EM t=1.0, **92.29% accuracy**, 70.0% mapping
- **Best consensus**: k12+k13+k14+k15 tn2 + EM, **90.07% accuracy**, 91.1% mapping
- **Best coverage**: k18+ tn=4, **86.79% accuracy**, 99.6% mapping
- **Species-level**: ~100% accuracy, **strain-level**: 60-90% within clade

While existing aligners (Bowtie2, BWA, minimap2) map reads to single reference genomes, Bit-Pop identifies **which genome** in a collection best matches each read — making it ideal for metagenomic classification tasks.

## Quick Start

```bash
# 1. Build (requires Rust: https://rustup.rs)
git clone https://github.com/mladenpop-oss/bit-pop.git
cd bit-pop
cargo build --release

# 2. One-command workflow: build index + map reads
./target/release/bit-pop run \
  data/genomes/Ecoli_K12_MG1655.fna \
  data/reads/simulated_ecoli_10k_new.fastq

# 3. Paired-end mode
./target/release/bit-pop run \
  data/genomes/Ecoli_K12_MG1655.fna \
  -1 data/reads/R1.fastq -2 data/reads/R2.fastq

# 4. Output BAM format (binary, compressed)
./target/release/bit-pop run \
  data/genomes/Ecoli_K12_MG1655.fna \
  data/reads/simulated_ecoli_10k_new.fastq \
  -o output.bam --bam

# 5. Download from NCBI and map
./target/release/bit-pop run \
  --ncbi "Escherichia coli" \
  data/reads/simulated_ecoli_10k_new.fastq
```

See [Usage](#usage) for full documentation.

## Features

### Core

- **Multi-genome indexing**: All reference genomes indexed in a single FM-index structure
- **Speed via bit-level operations**: 2-bit XOR alignment achieving ~2.3 ns per 31-base XOR chunk operation
- **Top-N rarest k-mer anchors**: Fallback to 2nd/3rd rarest k-mers for improved mapping rate (`--top-n`)
- **Combined ranking**: Formula balancing alignment score (85%) and k-mer rarity (15%)
- **Reverse complement support**: Full RC-aware mapping with proper SAM FLAG 0x10 handling
- **Paired-end support**: Full SAM specification compliance with proper FLAG handling
- **Native BAM output**: Binary alignment map format with BGZF compression (`--bam` flag)
- **Discordant pair reconciliation**: Automatic concordant genome resolution for R1/R2 cross-genome conflicts
- **Gaussian insert size model**: Probabilistic paired-end classification using normal distribution of observed insert sizes
- **Parallel mapping**: Work-stealing scheduler using rayon for multi-core speedup
- **Parallel index build**: Multi-threaded BWT and suffix array construction
- **Memory-mapped FASTA**: Reduced memory footprint with `--mmap` flag
- **Auto index caching**: Reuses `.bitpop` files when genomes haven't changed
- **NCBI integration**: Download genomes directly from NCBI with `--ncbi` flag
- **Progress reporting**: CLI progress bars for build and mapping operations
- **Smart defaults**: Automatic output paths, index detection, and progress reporting

### Post-Processing

- **EM post-processing**: Expectation-Maximization algorithm for multi-candidate refinement (`bit-pop em` command), +0.38% accuracy on k13 CAMI (91.91% → 92.29%)
- **Taxonomic classification**: NCBI taxonomy tree with LCA algorithm for genus/phylum/class-level abundance profiles (`bit-pop tax` command)

### 🔬 Experimental Features

The following features are available for testing purposes and may be useful for specific use cases (long reads, high-error data):

| Feature | Flag | Description |
|---------|------|-------------|
| Multi-k consensus | `bit-pop consensus` / `scripts/consensus_base.py` | Combine multiple k-mer indexes (best: k12-k15 tn2 + EM = 90.07% accuracy, 91% mapping) |
| Streaming mode | `--stream` | Process large FASTQ files with fixed RAM (~3GB per chunk) |
| Two-pass mapping | `--two-pass` | Re-map unmapped reads with lower threshold (+183 correct reads, CAMI 100k sample) |
| Fuzzy k-mer matching | `--method` | Fuzzy k-mer matching for error-prone reads |
| Soft-clipping | `-a softclip` | Adapter/low-quality region detection for reads with contamination |
| Gap-aware chaining | `-a chain` | Long-read alignment for ONT/PacBio data |
| Adaptive k-mer size | `--auto-k` | Automatic k-mer size selection based on genome size |

## Comparison with Existing Tools

| Feature | Bit-Pop | Bowtie2 | BWA-MEM | minimap2 |
|---------|---------|---------|---------|----------|
| Multi-genome classification | ✅ Native | ❌ Single genome | ❌ Single genome | ⚠️ With --index |
| Speed (10k reads, 3 genomes) | **0.9s** | ~5-10s | ~8-15s | ~3-5s |
| Index size (19.7 Mb) | **~152 MB** | ~200 MB | ~250 MB | ~180 MB |
| Quality-aware alignment | ✅ Phred-scaled | ✅ | ✅ | ✅ |
| Paired-end support | ✅ | ✅ | ✅ | ✅ |
| NCBI integration | ✅ Built-in | ❌ | ❌ | ❌ |
| Rust + bit-parallel | ✅ | C++ | C | C++ |

**When to use Bit-Pop**: Fast multi-genome classification where you need to identify which genome a read belongs to, rather than precise positional alignment.

## Bit-Pop vs Kraken2 — Different Tools for Different Use Cases

### Key Differences

| | Bit-Pop | Kraken2 |
|---|---------|---------|
| Database size | MB (only your genomes) | 100GB+ (entire NCBI) |
| Customization | Add/remove genomes in seconds | Fixed database, no customization |
| Build time | 2 minutes | Hours to days |
| Index growth | Grows only with your data | Pre-built massive database with unused data |

### When to Use Bit-Pop

- **Clinical microbiology** — A hospital tracks 20 strains. Build the index once, classify every patient sample in 0.13s.
- **Outbreak detection** — A new bacterium appears. Download one genome (MB), add to index, classify immediately.
- **Targeted analysis** — You know which genomes matter and want a lightweight, customizable solution.

Kraken2 is better for: broad metagenomics where you don't know what you're looking for.
Bit-Pop is better for: **targeted searching** where you know what matters.

## Pipeline

1. **FM-index** (SA-IS via libsais) for efficient k-mer lookup
2. **Anchor-based k-mer filtering** (top-N rarest k-mer selection with fallback)
3. **2-bit XOR alignment** (~2.3 ns per 31-base chunk for exact/near-exact matches)
4. **Multi-genome ranking** with combined scoring formula
5. **Reverse complement** scoring — tries both forward and RC, returns best match

## Installation

### Homebrew (macOS/Linux)

```bash
brew tap mladenpop-oss/homebrew-bit-pop
brew install bit-pop
```

### Cargo

```bash
cargo install bit-pop
```

### Docker

```bash
# Pull from GitHub Container Registry
docker pull ghcr.io/mladenpop-oss/bit-pop:latest

# Run with mounted data
docker run --rm -v $(pwd)/data:/home/bitpop/data ghcr.io/mladenpop-oss/bit-pop:latest \
  run /home/bitpop/data/genomes/Ecoli_K12_MG1655.fna \
      /home/bitpop/data/reads/simulated_ecoli_10k_new.fastq

# Or use the interactive shell
docker run --rm -it ghcr.io/mladenpop-oss/bit-pop:latest /bin/bash
```

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/mladenpop-oss/bit-pop/releases).

### From Source

#### Prerequisites

- Rust toolchain (2021 edition)

```bash
git clone https://github.com/mladenpop-oss/bit-pop.git
cd bit-pop
cargo build --release
```

### Optional Dependencies

- Python 3.x with Biopython - only required for read simulation (`scripts/simulate_reads.py`)

## Usage

### Global Flags

| Flag | Description |
|------|-------------|
| `-v, --verbose` | Enable verbose output |

### One-Command Workflow (Recommended)

```bash
# Single-end mode
./target/release/bit-pop run genome.fna reads.fastq

# Single-end with explicit reads flag
./target/release/bit-pop run genome.fna -r reads.fastq

# Paired-end mode
./target/release/bit-pop run genome.fna -1 R1.fastq -2 R2.fastq

# Multiple genomes from folder
./target/release/bit-pop run genomes/ reads.fastq

# Download from NCBI and map
./target/release/bit-pop run --ncbi "Escherichia coli" reads.fastq

# With custom options
./target/release/bit-pop run genome.fna -r reads.fastq \
  -o output.sam \
  -k 8 \
  -q 20 \
  -t 4
```

### Advanced Commands

#### Build Index

```bash
./target/release/bit-pop build \
  -f genome1.fasta -f genome2.fasta -f genome3.fasta \
  -o index.bitpop \
  -k 10 \
  -t 4

# CAMI mode — extract genome names from filenames
./target/release/bit-pop build \
  -f genome.fasta \
  -o index.bitpop \
  --cami \
  --mmap
```

#### Map Reads

```bash
# Single-end
./target/release/bit-pop map \
  -i index.bitpop \
  -r reads.fastq \
  -o output.sam \
  -a xor \
  -t 4

# Paired-end
./target/release/bit-pop map \
  -i index.bitpop \
  --reads-1 R1.fastq \
  --reads-2 R2.fastq \
  -o output.sam \
  -a xor \
  -t 4

# PacBio chunked mapping
./target/release/bit-pop map \
  -i index.bitpop \
  -r reads.fastq \
  -o output.sam \
  --chunk-size 1000 \
  --chunk-pct 0.8 \
  --chunk-min 50 \
  --chunk-max 200 \
  --chunk-vote-threshold 0.7 \
  --chunk-top-n 3

# Streaming mode (for huge FASTQ files, limits RAM usage)
./target/release/bit-pop map \
  -i index.bitpop \
  -r large_reads.fastq \
  -o output.sam \
  --stream \
  --max-ram 32G \
  -t 16

# Two-pass mapping (re-maps unmapped reads with lower threshold + EM refinement)
./target/release/bit-pop map \
  -i index.bitpop \
  -r reads.fastq \
  -o output.sam \
  --two-pass \
  --second-pass-score 0.4 \
  --top-n 4 \
  -t 16
```

#### Show Index Statistics

```bash
./target/release/bit-pop stats -i index.bitpop
```

#### Add Genomes to Existing Index

```bash
./target/release/bit-pop load \
  -i existing.bitpop \
  -f new_genome.fasta \
  -o updated.bitpop

# With memory-mapped I/O and parallel build
./target/release/bit-pop load \
  -i existing.bitpop \
  -f new_genome.fasta \
  -o updated.bitpop \
  --mmap \
  -t 4
```

**Parameters:**
- `-i, --index`: Existing index path (required)
- `-f, --fasta`: New genome FASTA file (required)
- `-o, --output`: Updated index path (required)
- `-t, --threads`: Number of threads (default: 1)
- `--mmap`: Use memory-mapped FASTA loading (feature-gated)

#### Search NCBI

```bash
./target/release/bit-pop search \
  --organism "Escherichia coli" \
  -n 10

# Filter by molecule type
./target/release/bit-pop search \
  --organism "Escherichia coli" \
  -m "genomic DNA" \
  -n 20

# With API key for higher rate limit
./target/release/bit-pop search \
  --organism "Escherichia coli" \
  -n 20 \
  --api-key YOUR_API_KEY \
  --email user@example.com
```

**Parameters:**
- `--organism`: Organism name to search (required)
- `-n, --limit`: Maximum number of results (default: 10)
- `-m, --molecule-type`: Filter by molecule type (e.g., "genomic DNA")
- `--api-key`: NCBI API key for higher rate limit
- `--email`: Email for NCBI request tracking

#### Fetch Genome from NCBI

```bash
./target/release/bit-pop fetch \
  --accession NC_000913.3 \
  -o index.bitpop

# Output FASTA instead of building index
./target/release/bit-pop fetch \
  --accession NC_000913.3 \
  -o output.fasta \
  --fasta-only

# With auto-k, custom cache, force re-download
./target/release/bit-pop fetch \
  --accession NC_000913.3 \
  -o index.bitpop \
  --auto-k \
  --cache-dir ./cache \
  --force
```

#### Update Cached Genomes

```bash
./target/release/bit-pop update

# Check specific index for updates
./target/release/bit-pop update \
  -i index.bitpop

# Force update all genomes
./target/release/bit-pop update \
  --force \
  --cache-dir ./cache

# With API key and email
./target/release/bit-pop update \
  --force \
  --api-key YOUR_API_KEY \
  --email user@example.com
```

**Parameters:**
- `-i, --index`: Check specific index (default: all cached)
- `--force`: Force re-download all genomes
- `--cache-dir`: Custom cache directory (default: ~/.cache/bitpop)
- `--api-key`: NCBI API key for higher rate limit
- `--email`: Email for NCBI request tracking

#### EM Post-Processing

Apply Expectation-Maximization algorithm to improve multi-candidate SAM mappings:

```bash
# Run EM on a SAM file produced by `bit-pop map`
./target/release/bit-pop em \
  -i mapped.sam \
  -o em_mapped.sam \
  --convergence 0.001 \
  --max-iterations 20 \
  --temperature 0.1 \
  --top-k 10

# With confidence threshold (recommended for near-identical strains)
./target/release/bit-pop em \
  -i mapped.sam \
  -o em_mapped.sam \
  --convergence 0.001 \
  --max-iterations 50 \
  --temperature 0.1 \
  --top-k 10 \
  --confidence-threshold 0.95
```

**What it does**: When a read maps to multiple genomes with similar scores, EM uses population-level abundance signals to reassign reads to the most likely genome. Typically converges in 6-7 iterations. On CAMI dataset, EM with `t=1.0, ct=0.95` on k13 mapping adds **+0.38% accuracy** (91.91% → 92.29%, 10,622 reassignments).

**Parameters**:
- `--convergence`: KL divergence threshold for stopping (default: 0.001)
- `--max-iterations`: Maximum EM iterations (default: 50)
- `--temperature`: Softmax temperature for probability smoothing (default: 0.1)
- `--top-k`: Number of top candidates per read for EM (default: 10)
- `--confidence-threshold`: Minimum probability to apply EM reassignment (default: 0.0)

#### Multi-K Consensus Mapping

Combine multiple indexes with different k-mer sizes for improved strain resolution:

```bash
# Multi-k consensus (k=8 + k=10 + k=12)
./target/release/bit-pop consensus \
  -i index_k8.bitpop:8 -i index_k10.bitpop:10 -i index_k12.bitpop:12 \
  -r reads.fastq \
  -o output.sam \
  -t 4

# With top-N candidates (recommended for strain resolution)
./target/release/bit-pop consensus \
  -i index_k8.bitpop:8 -i index_k10.bitpop:10 -i index_k12.bitpop:12 \
  -r reads.fastq \
  -o output.sam \
  --top-n 4 \
  -t 16

# Streaming mode (for huge FASTQ files, limits RAM usage)
./target/release/bit-pop consensus \
  -i index_k8.bitpop:8 -i index_k10.bitpop:10 -i index_k12.bitpop:12 \
  -r large_reads.fastq \
  -o output.sam \
  --stream \
  --max-ram 32G \
  -t 16

# With SNP detection
./target/release/bit-pop consensus \
  -i index_k8.bitpop:8,index_k10.bitpop:10 \
  -r reads.fastq \
  -o output.sam \
  --snp-detect \
  --snp-min-support 3 \
  --snp-penalty 0.1

# With chunked long-read support
./target/release/bit-pop consensus \
  -i index_k10.bitpop:10 \
  -r reads.fastq \
  -o output.sam \
  --chunk-size 86 \
  --chunk-pct 0.0 \
  --chunk-min 50 \
  --chunk-max 200
```

**K-Priority Weighting**: By default, consensus uses **k-priority weighting** where larger k-values get higher weight (`weight = k / min_k`). For k=8/10/12: weights are 1.0x, 1.25x, 1.5x respectively.

**Parameters**:
- `-i, --indexes`: List of `index:k` pairs (comma-separated, required)
- `-r, --reads`: Reads file (FASTQ, required)
- `-o, --output`: Output SAM file (required)
- `--strategy`: Voting strategy: `weighted_score` (default) or `majority`
- `--min-score`: Minimum alignment score threshold (default: 0.5)
- `--chunk-size`: Chunk size for long reads (default: 86)
- `--chunk-pct`: Chunk size as percentage (default: 0.0)
- `--chunk-min`: Minimum chunk size in bp (default: 50)
- `--chunk-max`: Maximum chunk size in bp (default: 200)
- `--snp-detect`: Enable SNP detection (default: false)
- `--snp-min-support`: SNP minimum support count (default: 3)
- `--snp-penalty`: SNP penalty value (default: 0.1)
- `--min-k-mappings`: Minimum k-values that must find a mapping (default: 1)
- `-t, --threads`: Number of threads (default: 1)
- `--top-n`: Number of top candidates to output per read (0 = only winner, default: 1)
- `--stream`: Process reads in chunks to limit memory usage (enables streaming mode)
- `--max-ram`: Max RAM to use for streaming (e.g., "32G", "16GB"). Auto-calculates optimal chunk size. Default: uses available system memory.
- `--two-pass`: Enable two-pass mapping — unmapped reads are re-mapped with a lower threshold, then refined with EM algorithm
- `--second-pass-score`: Minimum score threshold for second pass (default: 0.4, range: 0.0-1.0). Lower values map more reads but may introduce false positives. Recommended: 0.4 for best accuracy/recall trade-off
- `--diagnose-unmapped`: Sample 1000 unmapped reads and report why they failed (useful for debugging index/parameter issues)

#### Multi-K Consensus — Python Script

For faster multi-k consensus, use the standalone Python script which calls `bit-pop map` for each index separately, then combines results. This is **faster than built-in consensus** because it uses the optimized standalone map engine:

```bash
# Build indexes with different k values
bit-pop build -f genomes/ -o cami_k10.bitpop -k 10
bit-pop build -f genomes/ -o cami_k13.bitpop -k 13

# Run consensus via Python script
python scripts/consensus_base.py \
  --indexes cami_k10.bitpop,cami_k13.bitpop \
  --reads reads.fastq \
  --output consensus.sam \
  --strategy weighted_score \
  --threads 16
```

**Performance comparison** (~1M reads, CAMI, 61 genomes):

| Method | Mapped | Mapping Rate | Accuracy | EM Accuracy |
|--------|--------|--------------|----------|-------------|
| k13 tn=4 | 697,182 | 70.0% | 91.91% | **92.29%** |
| k14 tn=4 | 742,747 | 74.6% | 91.55% | - |
| k15 tn=4 | 840,561 | 84.5% | 89.48% | - |
| **k12-k15 tn2** | **909,493** | **91.0%** | **89.71%** | **90.07%** |
| k13+k22 tn2 | 994,188 | 99.5% | 89.50% | 89.86% |
| k13-k15+k22 tn2 | 994,214 | 99.98% | 89.58% | 89.39% |

**Winner**: k12-k15 consensus gives **+212k mapped reads** vs k13 alone with **-1.8% accuracy** (before EM). With EM: k13 solo = 92.29%, k12-k15 consensus = 90.07%. Adding k22+ increases coverage to ~100% but adds noise that lowers accuracy.

**Why use the Python script**: Built-in consensus uses `map_read()` in a parallel iterator, which is slower than the standalone `map` command. The Python script calls `bit-pop map` for each index, achieving better performance.

#### Chunk-Consensus Mapping

Use a single index with multiple chunk-% configurations and require voting agreement for higher accuracy:

```bash
# 3 chunk-% configs: 1%, 10%, 50%
./target/release/bit-pop chunk-consensus \
  -i index.bitpop \
  -c 0.01,0.10,0.50 \
  -r reads.fastq \
  -o output.sam \
  -t 8

# Custom min agreement (2 out of 3 configs must agree)
./target/release/bit-pop chunk-consensus \
  -i index.bitpop \
  -c 0.01,0.10,0.50 \
  -r reads.fastq \
  -o output.sam \
  --min-agreement 2 \
  --strategy majority \
  -t 8

# Weighted scoring (lower chunk-% → higher weight)
./target/release/bit-pop chunk-consensus \
  -i index.bitpop \
  -c 0.01,0.10,0.50 \
  -r reads.fastq \
  -o output.sam \
  --strategy weighted_score \
  -t 8
```

**Concept**: Smaller chunks (1%) map more reads but with lower accuracy. Larger chunks (50%) map fewer reads but with higher accuracy. Chunk-consensus requires a read to map to the **same genome** in at least N configs (default: majority = N/2+1), trading mapping rate for accuracy higher than any single configuration.

**Parameters**:
- `-i, --index`: Index file (.bitpop, required)
- `-c, --chunk-pcts`: Chunk percentages as fractions, comma-separated (e.g. "0.01,0.10,0.50", required)
- `-r, --reads`: Reads file (FASTQ, required)
- `-o, --output`: Output SAM file (required)
- `--strategy`: Voting strategy: `majority` (default) or `weighted_score`
- `--min-score`: Minimum alignment score threshold (default: 0.5)
- `--min-agreement`: Minimum configs that must agree (default: majority = N/2+1)
- `--chunk-min`: Minimum chunk size in bp (default: 50)
- `--chunk-max`: Maximum chunk size in bp (default: 200)
- `-t, --threads`: Number of threads (default: 1)
- `--top-n`: Number of top candidates per read (default: 1)

**SAM tags** (chunk-consensus specific):

| Tag | Type | Description |
|-----|------|-------------|
| `CP{pct}:Z:{name}/{score}` | string | Per-config result (e.g. `CP1:Z:Ecoli/0.9200`) |
| `CV:Z:{name}` | string | Final consensus genome |
| `CC:i:{n}` | integer | Total configs that found a mapping |
| `VC:i:{n}` | integer | Vote count (configs agreeing on winner) |
| `AS:f:{score}` | float | Consensus score |
| `XS:f:{score}` | float | Suboptimal score (supplementary only) |

**Example SAM line:**
```
read1   0       Ecoli     101     57      50M     *       0       0       ACGT...  *  NM:i:2  CP1:Z:Ecoli/0.9200  CP10:Z:Ecoli/0.8800  CP50:Z:Ecoli/0.9500  CV:Z:Ecoli  CC:i:3  VC:i:3  AS:f:0.9167
```

#### Taxonomic Classification Report

Generate taxonomic abundance profiles from SAM output using NCBI taxonomy tree with Lowest Common Ancestor (LCA) algorithm:

```bash
# Download NCBI taxonomy (once)
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# Extract nodes.dmp and names.dmp

# Generate taxonomic report from SAM output
./target/release/bit-pop tax \
  -i mapped.sam \
  --nodes-dmp nodes.dmp \
  --names-dmp names.dmp \
  -o taxonomy_report.txt

# JSON output
./target/release/bit-pop tax \
  -i mapped.sam \
  --nodes-dmp nodes.dmp \
  --names-dmp names.dmp \
  --format json \
  -o taxonomy_report.json

# Top 20 entries per rank
./target/release/bit-pop tax \
  -i mapped.sam \
  --nodes-dmp nodes.dmp \
  --names-dmp names.dmp \
  --top-n 20
```

**What it does:**
Maps genome-level read counts to NCBI taxonomy tree using LCA, producing taxonomic abundance profiles by rank (species, genus, family, phylum, etc.). Automatically resolves ambiguous mappings by finding the lowest common ancestor in the taxonomy tree.

**Example report output:**
```
=== Taxonomic Report (150 total reads) ===

FAMILY
Name                                   Reads        %
----------------------------------------------------------------
Enterobacteriaceae                         150    100.0%

GENUS
Name                                   Reads        %
----------------------------------------------------------------
Escherichia                               100     66.7%
Salmonella                                 50     33.3%

SPECIES
Name                                   Reads        %
----------------------------------------------------------------
Escherichia coli                          100     66.7%
Salmonella enterica                        50     33.3%
```

**Parameters:**
- `-i, --input`: Input SAM file (bit-pop mapping output, required)
- `--nodes-dmp`: Path to NCBI nodes.dmp file (required)
- `--names-dmp`: Path to NCBI names.dmp file (required)
- `-o, --output`: Output file (default: stdout)
- `--top-n`: Number of top entries per rank (default: 10)
- `--format`: Output format: `text` (default), `json`

## SAM/BAM Output Format

Bit-Pop produces standard SAM 1.6 output (text) or BAM 1.6 output (binary, BGZF-compressed) with additional optional tags for enhanced analysis. Use `--bam` flag for BAM output.

**SAM vs BAM**: BAM is the binary, compressed version of SAM. It contains the same data but is smaller and faster to read/write. Both formats are compatible with samtools, bcftools, IGV, and other bioinformatics tools.

### Standard SAM Fields

| Field | Description |
|-------|-------------|
| QNAME | Read name |
| FLAG | SAM flag (0, 4=UNMAPPED, 16=REVERSE, 2048=SUPPLEMENTARY) |
| RNAME | Reference genome name |
| POS | 1-based mapping position |
| MAPQ | Mapping quality (score × 60) |
| CIGAR | Alignment operations (M, X for XOR; M, I, D for Smith-Waterman; S for soft-clipping; I/D/S for chain mode gaps) |
| RNEXT/PNEXT/TLEN | Pair information (paired-end mode) |
| SEQ | Read sequence |
| QUAL | Quality string (quality-aware mode) |

### Optional Tags

| Tag | Type | Description |
|-----|------|-------------|
| `NM:i:` | integer | Edit distance (mismatches + insertions from CIGAR) |
| `AS:f:` | float | Alignment score (0.0-1.0, raw, before rarity weighting) |
| `MD:Z:` | string | Mismatching bases string (e.g., "10A5T3" = 10 matches, mismatch A, 5 matches, mismatch T, 3 matches) |
| `RK:f:` | float | K-mer rarity score (1 / occurrence_count of read's first k-mer) |
| `HF:f:` | float | Homopolymer fingerprint similarity (0.0-1.0, only when `--hf` enabled) |
| `GM:f:` | float | Gaussian insert size confidence (0.0-1.0, paired-end only, higher = more plausible insert size) |
| `XS:f:` | float | Suboptimal score (supplementary mappings only) |
| `MQ:f:` | float | Quality penalty (quality-aware mode only) |

### Example SAM Line

```
read1   0       chr1    101     57      50M     *       0       0       ACGTACGT...  *  NM:i:2  MD:Z:20A10T15  AS:f:0.9500  RK:f:0.001000  GM:f:0.8523  HF:f:0.8723
read1   2048    chr2    200     48      50M     *       0       0       ACGTACGT...  *  NM:i:3  MD:Z:15G20C10  AS:f:0.8000  RK:f:0.003000  GM:f:0.1205  HF:f:0.4521  XS:f:0.8000
```

### MD Tag — samtools/bcftools Compatibility

The `MD:Z:` tag enables full compatibility with the standard bioinformatics ecosystem.

**Native BAM output**: Use `--bam` flag to output BAM directly (no samtools conversion needed):

```bash
# Native BAM output
bit-pop run genome.fna reads.fastq -o output.bam --bam

# Or convert SAM to BAM (if using SAM output)
samtools view -b mapped.sam | samtools sort -o mapped.bam
samtools index mapped.bam

# SNP calling with bcftools
samtools mpileup -f reference.fasta mapped.bam | bcftools call -mv -o variants.vcf

# Visualize in IGV (MD tag enables base-level verification)
igv mapped.bam
```

**Why this matters:** Without the MD tag, tools like `samtools mpileup` and `bcftools call` cannot distinguish true mutations from sequencing errors. The MD tag specifies exactly which bases differ from the reference at each position, enabling:
- **SNP calling** — identify variants with bcftools/GATK
- **Base-level verification** — IGV shows mismatches in red with MD confirmation
- **Pileup analysis** — mpileup uses MD to count matches/mismatches per position
- **Cross-tool compatibility** — works with VarScan, FreeBayes, and all SAM-dependent tools

### Quality-Aware Mode Tags

When using quality scores (from FASTQ), an additional tag is included:

| Tag | Type | Description |
|-----|------|-------------|
| `MQ:f:` | float | Phred-scaled quality penalty (negative value for high-quality mismatches) |

### EM Post-Processing Output

The `bit-pop em` command reassigns reads based on population-level abundance signals. Output SAM preserves all optional tags (`AS`, `MD`, `RK`, `XS`) and updates:
- `RNAME` — changed if EM reassigns read to different genome
- `MAPQ` — set to 40 for reassigned reads

### Score Interpretation

| Tag | Range | Meaning |
|-----|-------|---------|
| `AS:f:` | 0.0–1.0 | Pure alignment quality (independent of k-mer rarity) |
| `RK:f:` | 0.0–1.0 | K-mer specificity (higher = rarer = more informative) |
| `HF:f:` | 0.0–1.0 | Homopolymer fingerprint similarity (only when `--hf` enabled) |
| `XS:f:` | 0.0–1.0 | Alternative mapping score (supplementary only) |

**Using AS + RK together:**
| AS | RK | Interpretation |
|----|-----|----------------|
| High | High | Strong, specific match — high confidence |
| High | Low | Good match but common k-mer — repeat region |
| Low | High | Rare k-mer but poor alignment — likely noise |
| Low | Low | Poor match, common k-mer — reject |

**Using HF for strain resolution:**
| AS | HF | Interpretation |
|----|-----|----------------|
| High | High | Perfect strain match — high confidence assignment |
| High | Low | Good alignment but wrong strain — k-mer ambiguity |
| Low | High | Homopolymer match but poor alignment — partial match |
| Low | Low | No match — reject |

### `run` Command Options

| Flag | Description | Default |
|------|-------------|---------|
| `genome` | Genome file, folder, or NCBI organism | (required) |
| `-r, --reads` | Reads file for single-end mode | (required) |
| `-1, --reads-1` | R1 FASTQ for paired-end | (required with -2) |
| `-2, --reads-2` | R2 FASTQ for paired-end | (required with -1) |
| `--ncbi` | Fetch genome from NCBI | false |
| `-o, --output` | Output SAM file | `<reads_name>.sam` |
| `-k, --k` | K-mer size | 10 |
| `--auto-k` | Auto-calculate optimal k-mer size | false |
| `--read-type` | Read type: short (clamp [10,15]) / long (clamp [13,19]) | short |
| `-s, --spaced-seed` | Enable spaced seed matching | false |
| `-a, --align-mode` | Alignment mode: xor, sw, hybrid | hybrid |
| `-m, --min-score` | Minimum alignment score (0.0-1.0) | 0.7 |
| `-q, --min-quality` | Minimum Phred quality (0 = no filter) | 0 |
| `-t, --threads` | Number of threads | 1 |
| `--top-n` | Top N rarest k-mer anchors (higher = better mapping rate, slower) | 1 |
| `--reconcile-top-n` | Top N candidates per read for discordant pair reconciliation (paired-end only) | 5 |
| `--mmap` | Use memory-mapped FASTA loading | false |
| `--force` | Force rebuild index | false |
| `--method` | Fuzzy k-mer method: none, fuzzy-kmer, fuzzy-seed, neighborhood | none |
| `--fuzzy-mismatches` | Max mismatches for fuzzy matching | 1 |
| `--spaced-seed-pattern` | Custom spaced seed pattern string | None |
| `--golden-anchors` | Quality-weighted k-mer anchors for long reads | false |
| `--chunk-size` | Chunk size for PacBio long-read mapping | None |
| `--chunk-pct` | Chunk size as percentage of read length | None |
| `--chunk-min` | Minimum chunk size clamp | None |
| `--chunk-max` | Maximum chunk size clamp | None |
| `--chunk-vote-threshold` | Minimum fraction of chunks that must agree | None |
| `--chunk-top-n` | Number of top genomes to return per read in chunk mode | None |
| `--snp-detect` | Enable SNP-aware scoring | false |
| `--snp-min-support` | Minimum support count for SNP detection | 3 |
| `--hf` | Enable homopolymer fingerprint scoring | false |
| `--hf-min` | Minimum run length for homopolymer fingerprint | 3 |
| `-i, --index` | Use existing .bitpop index (skip genome loading) | (none) |
| `--em` | Apply EM post-processing after mapping | false |
| `--search-radius` | Search radius in bp around anchor (±N, default: 5, max: 200) | 5 |
| `--chunk-strategy` | Chunk anchor strategy: rarest, golden, spaced | rarest |
| `--api-key` | NCBI API key for higher rate limit | (none) |
| `--email` | Email for NCBI request tracking | (none) |
| `--bam` | Output BAM format instead of SAM | false |

### `build` Command Options

| Flag | Description | Default |
|------|-------------|---------|
| `-f, --fasta` | Input FASTA file(s) (required, can be repeated) | (required) |
| `-o, --output` | Output index path (required) | (required) |
| `-k, --k` | K-mer size | 8 |
| `--auto-k` | Auto-scale k-mer size based on genome size | false |
| `--read-type` | Read type: short (Illumina) / long (Nanopore/PacBio) | short |
| `-t, --threads` | Number of threads for parallel build | 1 |
| `--method` | Fuzzy k-mer method: none, fuzzy-kmer, fuzzy-seed, neighborhood | none |
| `--fuzzy-mismatches` | Max mismatches for fuzzy matching | 1 |
| `-s, --spaced-seed` | Enable spaced seed pattern matching | false |
| `--spaced-seed-pattern` | Custom spaced seed pattern (e.g., "11101001110111") | 11111011111111 |
| `--cami` | Extract genome name from filename (CAMI mode) | false |
| `--pacbio` | Extract genome name from filename (PacBio mode) | false |
| `--search-radius` | Search radius in bp around anchor (±N, default: 5, max: 200) | 5 |
| `--hf` | Enable homopolymer fingerprint scoring | false |
| `--hf-min` | Minimum run length for homopolymer fingerprint | 3 |
| `--mmap` | Use memory-mapped FASTA loading | false |

### `build` Command — CAMI Dataset Support

For CAMI benchmark datasets, use the `--cami` flag during index build to extract genome names from filenames instead of FASTA headers:

```bash
# CAMI dataset — genome names extracted from filenames
bit-pop build --cami -f 1036554.gt1kb.fasta -o index.bitpop
# → "1036554.gt1kb.fasta" → genome name: "1036554"

# evo_* strains — .NNN suffix preserved
# → "evo_1049056.011.fna" → genome name: "evo_1049056.011"
```

This fixes accuracy for CAMI datasets where FASTA headers don't match ground truth labels. Without `--cami`, accuracy can be as low as 1.07%; with it, accuracy reaches ~93% on the CAMI Low Complexity benchmark.

### Align Modes

- `xor`: Fast 2-bit XOR alignment only
- `sw`: Smith-Waterman refinement for all reads
- `hybrid`: XOR first, SW only when confidence < 0.9
- `softclip`: XOR with soft-clipping — slides windows across read to find optimal alignment region, emitting `S` operations in CIGAR for adapter/low-quality regions. Uses two-pass strategy (coarse scan + fine-grained refinement) for O(N * step) performance. Ideal for reads with adapter contamination or spanning repeat boundaries.
- `chain`: Gap-aware XOR chaining — true long-read alignment for ONT/PacBio. Uses minimizer-based seed chaining with XOR gap extension. Handles 5-15% error rates and long indels natively, replacing chunk-based workaround with proper minimizer chaining (similar to minimap2's approach).

**Soft-clipping example:**
```
Read:  |ADAPTER (10bp)|ACTUAL_READ (50bp)|NOISE (5bp)|
CIGAR: 10S50M5S
```

**Chain mode vs chunk mode:**
| Aspect | Chunk mode (old) | Chain mode (new) |
|--------|------------------|------------------|
| Approach | Fixed-size chunks + voting | Minimizer seeds + collinear chaining |
| Gaps/indels | Not handled | Handled via chaining tolerance |
| Error rate | Works best <5% | Handles 5-15% (ONT/PacBio) |
| Speed | O(chunks × genomes) | O(minimizers × log(hits)) |
| CIGAR | Generic M | Detailed M/X/I/D/S |

**Chain mode configuration:**
```bash
# Default chain mode (k=15, w=10, min_seeds=3)
bit-pop map -i index.bitpop -r reads.fastq -a chain

# Custom chain config for high-error ONT data
bit-pop map -i index.bitpop -r ont_reads.fastq -a chain \
  --chain-k 15 \
  --chain-w 10 \
  --chain-min-seeds 3 \
  --chain-max-gap 500 \
  --chain-gap-open -5 \
  --chain-gap-extend -0.5
```

### `map` Command Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --index` | Input index path (required) | (required) |
| `-r, --reads` | Reads file for single-end mode | (none) |
| `-1, --reads-1` | R1 FASTQ for paired-end | (none) |
| `-2, --reads-2` | R2 FASTQ for paired-end | (none) |
| `-o, --output` | Output SAM file (required) | (required) |
| `-a, --align-mode` | Alignment mode: xor, sw, hybrid | xor |
| `-m, --min-score` | Minimum alignment score (0.0-1.0) | 0.7 |
| `-q, --min-quality` | Minimum Phred quality (0 = no filter) | 0 |
| `-t, --threads` | Number of threads | 1 |
| `--top-n` | Top N rarest k-mer anchors | 1 |
| `--reconcile-top-n` | Top N candidates per read for discordant pair reconciliation (paired-end only) | 5 |
| `--method` | Fuzzy k-mer method: none, fuzzy-kmer, fuzzy-seed, neighborhood | none |
| `--fuzzy-mismatches` | Max mismatches for fuzzy matching | 1 |
| `-s, --spaced-seed` | Enable spaced seed matching | false |
| `--spaced-seed-pattern` | Custom spaced seed pattern | 11111011111111 |
| `--golden-anchors` | Use golden anchor selection (quality-weighted) | false |
| `--em` | Apply EM post-processing after mapping | false |
| `--search-radius` | Search radius in bp around anchor (±N, default: 5, max: 200) | 5 |
| `--chunk-size` | Chunk size for PacBio long-read mapping | 0 (auto) |
| `--chunk-pct` | Chunk size as percentage of read length (0.0-1.0) | 0.0 |
| `--chunk-min` | Minimum chunk size clamp for dynamic chunking | 50 |
| `--chunk-max` | Maximum chunk size clamp for dynamic chunking | 200 |
| `--chunk-vote-threshold` | Minimum fraction of chunks that must agree (0.0-1.0) | 0.0 |
| `--chunk-top-n` | Number of top genomes per read in chunk mode | 1 |
| `--chunk-strategy` | Chunk anchor strategy: rarest, golden, spaced | rarest |
| `--snp-detect` | Enable SNP-aware scoring | false |
| `--snp-min-support` | Minimum support count for SNP detection | 3 |
| `--hf` | Enable homopolymer fingerprint scoring | false |
| `--hf-min` | Minimum run length for homopolymer fingerprint | 3 |
| `--bam` | Output BAM format instead of SAM | false |

## Benchmark Results

### 3-Genome Benchmark (Simulated Reads)

**Setup**: E. coli K-12 MG1655 (4.6 Mb), S. aureus (2.9 Mb), S. cerevisiae (12.2 Mb). 20,000 simulated reads (100 bp, 0.1% error rate, Q30-Q40). k=10. Standard desktop CPU.

#### Results (k=10, top_n=1)

| Genome | Size | Mapped | Mapping Rate | Accuracy |
|--------|------|--------|--------------|----------|
| E. coli | 4.6 Mb | 9,755/10,000 | 97.6% | 99.9% |
| S. aureus | 2.9 Mb | 4,910/5,000 | 98.2% | 99.9% |
| S. cerevisiae | 12.2 Mb | 4,905/5,000 | 98.1% | 100.0% |
| **Total** | **19.7 Mb** | **19,570/20,000** | **97.9%** | **99.9%** |

#### Results (k=10, top_n=3)

| Genome | Size | Mapped | Mapping Rate | Accuracy |
|--------|------|--------|--------------|----------|
| E. coli | 4.6 Mb | 9,924/10,000 | 99.2% | 99.9% |
| S. aureus | 2.9 Mb | 4,968/5,000 | 99.4% | 99.9% |
| S. cerevisiae | 12.2 Mb | 4,970/5,000 | 99.4% | 100.0% |
| **Total** | **19.7 Mb** | **19,862/20,000** | **99.3%** | **99.9%** |

**Performance trade-off**: top_n=3 is ~3x slower than top_n=1 (2.8s vs 0.9s for E. coli). Recommended: `--top-n 2` for balance between speed and accuracy.

**Throughput**: ~1,500 reads/second (top_n=1)

### CAMI Low Complexity Benchmark (61 Genomes, 1M Reads)

**Setup**: 61 genomes (1900 sequences, ~157M bases, ~1.3GB index). ~1M reads sampled from original 30GB CAMI dataset. k=10, 16 threads, `--cami` flag for genome naming.

#### K-mer Size Sweep (xor alignment, top-n 4)

| k | Mapped | Mapping Rate | Accuracy | Time |
|---|--------|--------------|----------|------|
| 10 | 857,983 (86.2%) | 86.2% | 86.49% | ~60s |
| 11 | 816,490 (82.0%) | 82.0% | 87.01% | ~55s |
| 12 | 721,026 (72.4%) | 72.4% | 89.83% | ~55s |
| **13** | **697,182 (70.0%)** | **70.0%** | **91.91%** | **~50s** |
| 14 | 742,747 (74.6%) | 74.6% | 91.55% | ~55s |
| 15 | 840,561 (84.5%) | 84.5% | 89.48% | ~60s |
| 16 | 940,303 (94.5%) | 94.5% | 87.28% | ~65s |
| 17 | 980,328 (98.6%) | 98.6% | 86.77% | ~70s |
| 18 | 990,721 (99.6%) | 99.6% | 86.79% | ~75s |
| 19 | 993,265 (99.9%) | 99.9% | 86.82% | ~75s |
| 20 | 993,894 (99.95%) | 99.95% | 86.83% | ~75s |
| 21 | 994,079 (99.97%) | 99.97% | 86.84% | ~75s |
| 22 | 994,122 (99.97%) | 99.97% | 86.84% | ~75s |

**Pattern**: k13 gives peak accuracy (91.91%), then accuracy drops and plateaus at ~86.8% for k≥18. Coverage increases monotonically with k, reaching ~100% at k≥18.

#### Top-N Series (xor alignment, k=10)

| Config | Mapped | Mapping Rate | Accuracy | Time |
|--------|--------|--------------|----------|------|
| xor top-n 1 | 638,704 (63.9%) | 63.9% | 86.15% | 129s |
| xor top-n 2 | 748,711 (74.9%) | 74.9% | 86.31% | 189s |
| xor top-n 3 | 813,874 (81.4%) | 81.4% | 86.41% | 250s |
| xor top-n 4 | 861,762 (86.2%) | 86.2% | 86.48% | 132s |
| xor top-n 5 | 712,720 (71.3%) | 71.3% | 86.50% | ~300s |
| xor top-n 6 | 511,411 (51.2%) | 51.2% | 86.57% | ~300s |

#### EM Post-Processing (on top-n 4)

| Config | Mapped | Accuracy | Time |
|--------|--------|----------|------|
| xor tn=4 | 861,762 | 86.48% | 132s |
| xor tn=4 + EM (t=0.1) | 857,992 | 86.58% | ~132s |
| **xor tn=4 + EM (t=0.1, ct=0.95)** | **857,992** | **86.77%** | **~132s** |
| xor tn=4 + EM (t=0.5, ct=0.95) | 857,992 | 85.29% | ~132s |

#### Two-Pass Mapping (re-maps unmapped reads with lower threshold + EM)

Two-pass mapping recovers reads that fail the initial 0.7 threshold by re-mapping them with a lower threshold, then refining with EM. Best at threshold 0.4:

| Config | Mapped | Accuracy | Correct | Wrong |
|--------|--------|----------|---------|-------|
| xor tn=4 | 49,028 | 86.52% | 42,421 | 6,607 |
| xor tn=4 + 2pass (0.5) | 49,071 | 86.48% | 42,436 | 6,635 |
| **xor tn=4 + 2pass (0.4) + EM** | **49,929** | **85.33%** | **42,604** | **7,325** |
| xor tn=4 + 2pass (0.3) + EM | 50,000 | 85.08% | 42,540 | 7,460 |

**Trade-off**: Two-pass gains +183 correct reads at the cost of +718 wrong reads. Use when maximizing recall is more important than precision.

#### EM Post-Processing (k13, top-n 4)

| Config | Mapped | Accuracy | Changed | Time |
|--------|--------|----------|---------|------|
| k13 tn=4 | 697,182 | 91.91% | - | ~50s |
| **k13 tn=4 + EM (t=1.0, ct=0.95)** | **697,188** | **92.29%** | **10,622** | **~60s** |
| k13 tn=4 + EM (t=0.1, ct=0.95) | 697,188 | 91.93% | 9,980 | ~60s |

**Note**: Temperature t=1.0 is better than t=0.1 for k13 — softer probability distribution allows more correct reassignments.

#### Multi-K Consensus (consensus_base.py, top-n 4, EM t=1.0)

| Config | Mapped | Mapping Rate | Accuracy | EM Accuracy | Changed |
|--------|--------|--------------|----------|-------------|---------|
| k13 tn=4 | 697,182 | 70.0% | 91.91% | **92.29%** | 10,622 |
| k12+k13 tn2 | 768,639 | 77.0% | 89.75% | **90.10%** | 132,153 |
| **k12+k13+k14+k15 tn2** | **909,493** | **91.0%** | **89.71%** | **90.07%** | **132,153** |
| k13+k22 tn2 | 994,188 | 99.5% | 89.50% | 89.86% | 184,078 |
| k13+k14+k15+k22 tn2 | 994,214 | 99.98% | 89.58% | 89.39% | 190,020 |

**Pattern**: Adding higher k (k22+) increases coverage to ~100% but introduces noise that lowers accuracy. k12-k15 range is optimal for consensus.

#### 🏆 Best Configuration

| Goal | Config | Mapped | Mapping Rate | Accuracy | Time |
|------|--------|--------|--------------|----------|------|
| **Best accuracy** | k13 tn=4 + EM t=1.0 | 697,188 | 70.0% | **92.29%** | ~60s |
| **Best consensus** | k12-k15 tn2 + EM t=1.0 | 909,493 | 91.0% | **90.07%** | ~250s |
| **Best coverage** | k18 tn=4 | 990,721 | 99.6% | 86.79% | ~80s |

#### Per-Genome Accuracy Breakdown (k12-k15 tn2 + EM t=1.0)

- Numeric genomes (1030752, 1036554, etc.): **99-100%**
- Sample genomes: **100%**
- 1052944, 1053058, 1052947: **36-73%** (closely related group)
- 1286_AP parent: **20.87%** (confused with evo_1286_AP strains)
- evo_* strains: **24-88%** (near-identical, fundamental limitation)
  - evo_1035930.011: 78.77%
  - evo_1035930.029: 88.14%
  - evo_1035930.032: 41.26%
  - evo_1049056.011: 79.54%
  - evo_1049056.013: 30.14%
  - evo_1049056.015: 42.03%
  - evo_1049056.031: 30.02%
  - evo_1049056.039: 24.24%
  - evo_1286_AP.008: 83.67%
  - evo_1286_AP.026: 84.57%
  - evo_1286_AP.033: 42.52%
  - evo_1286_AP.037: 63.74%

**Key insight**: Misclassifications are **within-clade only**. Reads from evo_1035930.* never map to evo_1049056.* or evo_1286_AP.*. Species-level classification is ~100% accurate. Strain-level confusion only occurs between genomes sharing >99.9% identity.

#### Unmapped Reads Analysis

**137,700 unmapped reads** (13.8% of 999K FASTQ reads). Parent genomes (1139_AG, 1220_AD, 1030752) have **2x higher unmapped rate** (18-19%) vs evo_* strains (9-10%). Cause: larger genomes → more repetitive regions → more common k-mers → reads fail k-mer rarity threshold.

**Diagnosis**: Use `--diagnose-unmapped` to sample 1000 unmapped reads and report why they failed. Common causes:
- `K-mers in index, alignment failed` — read has matching k-mers but alignment score below threshold (use `--two-pass` to recover)
- `No k-mers in index` — read sequence not represented in reference genomes
- `All k-mers too repetitive` — k-mers appear too many times to be useful anchors

#### Key Findings

1. **`--top-n` is the only flag that significantly affects accuracy** — growth is linear but slow (+0.16% for tn1→tn2, +0.07% for tn4→tn5)
2. **Mapping rate is fixed per top-n** — other flags do not change how many reads map
3. **All "advanced" flags** (HF, SNP, golden, spaced-seed, search-radius, chunk-strategy) **have no effect** on accuracy
4. **Spaced seeds are catastrophic** — only 413-807 mapped reads (vs 748K baseline)
5. **SW mode is too slow** — timeout after 600s
6. **Diminishing returns** from tn=5 onward — mapped count drops drastically (712K→511K)
7. **EM adds +0.38% accuracy on k13** (91.91% → 92.29%) via reassignment within strain groups (10,622 reassignments). Temperature t=1.0 is better than t=0.1
8. **k-mer size sweep (k10-k22)**: k13 gives peak accuracy (91.91%), k18+ gives peak coverage (~100%) but accuracy drops to ~86.8% plateau
9. **Multi-k consensus: k12-k15 is optimal** — 90.07% accuracy with EM, 91% mapping. Adding k22+ increases coverage to ~100% but adds noise (accuracy drops to 89.4%)
10. **Two-pass mapping** (`--two-pass`) recovers unmapped reads at threshold 0.4: +183 correct reads, +718 wrong reads, 99.9% mapping rate
11. **Two-pass threshold**: 0.4 is optimal — 0.5 maps too few additional reads, 0.3 maps all but adds more false positives
12. **Species-level classification: ~100%** — misclassifications only occur within clades sharing >99.9% identity, never between different species
13. **Strain-level classification: 60-90%** — evo_* strains confused with parent and sibling strains. Weighted avg: 61.8%. Larger strains (evo_029, evo_008) classify better (84-88%)
14. **Unmapped reads are not "on the edge"** — scores of 0.35-0.49, far below 0.7 threshold. These are strain variants not in the reference

#### Why is overall accuracy lower than single-genome benchmarks?

The evo_* genomes are >99.9% identical strains from the same sample assembly. They share most k-mers with each other, causing reads to map to the wrong strain. This is a **fundamental limitation** of k-mer-based classification for near-identical genomes, not a bug. **Species-level classification remains ~100% accurate** — misclassifications only occur within clades, never between different species. SNP-aware weighting or ML would be required for strain-level resolution.

**See**: [docs/paper.pdf](docs/paper.pdf) for detailed analysis.

### PacBio HiFi Benchmark (69 Genomes, 86k Reads)

**Setup**: 69 bacterial genomes (51+ species, 285 Mb). 86,248 simulated PacBio HiFi reads (8-20 kb, realistic error profile: 0.1% base errors, 2% homopolymer errors, 1% chimeras, variable coverage ±50%). k=70, 16 threads.

#### Accuracy vs k-mer Size

| k | Accuracy | Mapping Rate | Map Time |
|---|----------|-------------|----------|
| 10 | 42.6% | 99.97% | - |
| 13 | 73.7% | 99.97% | - |
| 15 | 75.1% | 99.97% | - |
| 20 | 79.7% | 99.97% | - |
| 25 | 82.7% | 99.97% | 7.5 min |
| 30 | 82.9% | 99.97% | 7.9 min |
| 40 | 83.0% | 99.97% | 8.0 min |
| 70 | 83.1% | 99.97% | 7.7 min |
| **70 (no chunk)** | **95.7%** | **99.93%** | **6.9 min** |

#### Realistic vs Simple Error Profile

| Dataset | k | Accuracy | Mapping Rate |
|---------|---|----------|-------------|
| Simple (0.1% errors) | 70 | 95.7% | 99.93% |
| **Realistic (homopolymers, chimera, coverage)** | **70** | **95.2%** | **99.0%** |

**Key findings:**
- Accuracy plateaus at k≥40 for chunked mapping
- **No-chunk mode** (full read alignment) gives **+12.6% accuracy** over chunked mode
- Realistic error profile (homopolymers, chimera, coverage variation) causes only **-0.5% accuracy drop**
- 8 minutes to map 86k long reads on 16 threads
- ~100% mapping rate across all k values

**See**: [simulate_realistic.py](simulate_realistic.py) for read simulation script.

## Project Structure

```
├── src/                    # Rust source code (17 modules)
│   ├── main.rs             # CLI entry point (12 subcommands)
│   ├── lib.rs              # Core library (BitPop struct, DNA encoding)
│   ├── fm.rs               # FM-index (SA-IS, BWT, backward search)
│   ├── align.rs            # Alignment (XOR, SW, Myers)
│   ├── sam.rs              # SAM output format
│   ├── em.rs               # EM post-processing algorithm
│   ├── taxonomy.rs         # NCBI taxonomy tree + LCA algorithm
│   ├── consensus.rs        # Multi-k consensus mapping
│   ├── chunk_consensus.rs  # Multi chunk-% consensus voting
│   ├── snp.rs              # SNP detection and scoring
│   ├── fasta.rs            # FASTA parsing + memory-mapped reader
│   ├── fastq.rs            # FASTQ parsing + quality filtering
│   ├── rank.rs             # Multi-genome ranking
│   ├── ncbi.rs             # NCBI E-utilities API client
│   ├── cache.rs            # Local cache management
│   ├── index_manager.rs    # Dynamic index management
│   ├── delta.rs            # Delta encoding + VLI compression
│   ├── persisted.rs        # Advanced persistence (memmap2, format v5)
│   └── serialize.rs        # Binary serialization
├── bin/                    # Additional CLI tools
│   └── extract_seqs.rs     # Extract genome sequences from .bitpop index
├── benches/                # Criterion benchmarks (17 benchmark groups)
├── tests/                  # Integration tests (5 tests)
├── scripts/
│   ├── simulate_reads.py       # Read simulation (Biopython)
│   ├── analyze_benchmark_new.ps1 # Benchmark analysis
│   ├── bitpop-workflow.py      # Multi-index workflow tool
│   ├── consensus_base.py       # Multi-k consensus (calls standalone map, faster than built-in)
│   ├── cami_accuracy.py        # Evaluate SAM accuracy against CAMI ground truth
│   └── em_classifier.py        # Python EM prototype (reference implementation)
├── data/
│   ├── genomes/            # Reference genomes (.fna, .fasta)
│   └── reads/              # Sequencing reads (.fastq)
├── docs/
│   ├── paper.tex           # Academic paper
│   ├── paper.pdf           # Compiled paper
│   ├── references.bib      # Bibliography
│   └── CITATION.cff        # Citation metadata
├── bioconda-recipe/        # Conda package recipe (future)
├── Cargo.toml              # Rust project configuration
└── README.md               # This file
```

### Data Files

**Genomes:**
- `data/genomes/Ecoli_K12_MG1655.fna` - E. coli K-12 MG1655 (4.6 Mb)
- `data/genomes/CP029198.1.fasta` - Staphylococcus aureus (2.9 Mb)
- `data/genomes/Sac_cerevisiae_complete.fasta` - S. cerevisiae S288C (12.2 Mb)

**Simulated Reads:**
- `data/reads/simulated_ecoli_10k_new.fastq` - 10,000 E. coli reads
- `data/reads/simulated_aureus_5k_new.fastq` - 5,000 S. aureus reads
- `data/reads/simulated_cerevisiae_5k_new.fastq` - 5,000 S. cerevisiae reads

## Testing

```bash
# Run all tests (unit + integration)
cargo test

# Run only integration tests
cargo test --test integration_tests

# Run benchmarks
cargo bench
```

**Test coverage:**
- 312+ unit tests (alignment, indexing, serialization, SAM output, spaced seeds, delta encoding, persistence, EM algorithm, taxonomy, chain mode)
- 5 integration tests (build, map, multi-genome, SAM format, cache reuse)
- 17 Criterion benchmark groups (XOR, SW, Myers, FM-index, k-mer filter, full pipeline)

## Limitations

- Research tool; not validated for clinical use
- Index file sizes ~152 MB for 19.7 Mb genome
- Chunked reads (>31bp) use generic CIGAR without per-base mismatch detail
- **Strain-level resolution**: Genomes that are >99.9% identical (same sample, different strains) share most k-mers. Reads may map to the wrong strain. This is a **fundamental information-theoretic limitation** of k-mer rarity-based classification, not a bug. EM post-processing can consolidate reads within strain groups (+0.38% on k13) but cannot fully resolve sibling strains that share >99.9% of their k-mers. **Species-level classification is ~100% accurate** — misclassifications only occur within clades, never between different species. SNP-aware weighting or ML would be required for full strain-level resolution.

## Large Genome Support

**Limitation:** FM-index construction uses libsais which has a ~2GB limit per index (~2.1B characters).

**Solution for large genomes (>2GB):** Use the workflow tool to automatically split, build, map, and merge:

```bash
# Full workflow (all steps automatic)
python scripts/bitpop-workflow.py full genome.fna reads.fastq -o output/ --threads 8

# Or manual step-by-step:
python scripts/bitpop-workflow.py split genome.fna -o chunks/
python scripts/bitpop-workflow.py build chunks/ -o indexes/ --threads 8
python scripts/bitpop-workflow.py map indexes/ reads.fastq -o mapped/ --threads 8
python scripts/bitpop-workflow.py merge mapped/ -o final.sam
```

**How it works:**
1. Splits genome into chunks (< 2GB each) by accession/chromosome boundaries
2. Builds FM-index for each chunk in parallel
3. Maps reads against all indexes in parallel
4. Merges SAM results (deduplicates by read name)

**Options:**
- `--max-size 2000` - max chunk size in MB (default: 2000)
- `--threads 8` - parallel threads (default: 4)
- `--no-cleanup` - keep intermediate files

---

## Development Roadmap

### ✅ Completed

- **Phase 0**: Critical bug fixes (rarity calculation, TLEN, BWT serialization, panic fixes)
- **Phase 1.1**: Top-N rarest k-mer anchors (97.9% → 99.3% mapping rate)
- **Phase 1.2**: Reverse complement support with SAM FLAG 0x10
- **Phase 1.3**: Paired-end support with full SAM compliance
- **Phase 1.4**: Myers edit distance (23-54x faster than Smith-Waterman)
- **Phase 2.1**: Memory-mapped FASTA (`--mmap`)
- **Phase 2.2**: Parallel index build (rayon)
- **Phase 3.1**: Progress reporting (CLI progress bars)
- **Phase 4**: CAMI Low Complexity benchmark (61 genomes, ~1M reads, k10-k22 sweep, consensus, EM, 92.29% peak accuracy)
- **Phase 6**: NCBI E-utilities integration (search, fetch, update commands)
- **Phase 7**: Large genome workaround (`bitpop-workflow.py`)
- **UX**: `run` command with auto-index caching and smart defaults
- **Tests**: Integration test suite (5 tests)
- **Paired-end**: Discordant pair reconciliation — resolves R1/R2 cross-genome conflicts via top-N candidate overlap
- **BAM output**: Native binary alignment map format with BGZF compression (`--bam` flag)
- **Taxonomic classification**: NCBI taxonomy tree with LCA algorithm (`bit-pop tax` command)
- **Gaussian insert size model**: Probabilistic paired-end classification using normal distribution of observed insert sizes
- **EM post-processing**: Expectation-Maximization for multi-candidate refinement (+0.29% on CAMI)

### 🔧 In Progress

- **Strain resolution**: Investigating approaches for >99.9% identical genomes

### 📋 Planned

- **Phase 2**: SA compression, streaming input, SIMD acceleration (AVX2)
- **Phase 3**: CIGAR accuracy improvements, quality filter enhancements
- **Phase 5**: Read caching, enhanced statistics, API documentation (docs.rs)
- **Multi-index**: Unified FM-index with automatic splitting (>2GB genomes)
- **Strain resolution**: Multi-k consensus, long-read support (PacBio/ONT), known SNP (VCF) integration

### 📊 Expand Benchmarks

- 100+ genomes and eukaryotic genomes
- Direct comparison with Bowtie2, BWA-MEM on multi-genome tasks
- CAMI Low Complexity: completed (61 genomes, ~1M reads, 89.85% accuracy)

## Getting Help

- **Documentation**: This README and [docs/paper.pdf](docs/paper.pdf)
- **Issues**: [GitHub Issues](https://github.com/mladenpop-oss/bit-pop/issues) — bug reports and feature requests
- **Discussions**: [GitHub Discussions](https://github.com/mladenpop-oss/bit-pop/discussions) — questions and feature ideas
- **Citation**: See [CITATION.cff](docs/CITATION.cff) or the DOI below

## Paper

[Read the full paper (PDF)](docs/paper.pdf)

## Availability

Source code available under the MIT License.

## Citation

```bibtex
@software{popovic_2026_bitpop,
  author = {Popovi{\'c}, Mladen},
  title = {Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification},
  year = {2026},
  doi = {10.5281/zenodo.20043593},
  url = {https://github.com/mladenpop-oss/bit-pop}
}
```

Or plain text:

> Popović, M. (2026). Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification. https://doi.org/10.5281/zenodo.20043593

## License

MIT License
