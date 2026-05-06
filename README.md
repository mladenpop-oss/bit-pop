# Bit-Pop: Multi-Genome DNA Read Classification

A proof-of-concept genomic tool for simultaneous multi-genome DNA read classification. While existing aligners such as Bowtie2, BWA, and minimap2 map reads to single reference genomes, Bit-Pop identifies which genome in a collection best matches each read.

## Features

- **Multi-genome indexing**: All reference genomes indexed in a single FM-index structure
- **Speed via bit-level operations**: 2-bit XOR alignment achieving ~2.3 ns per 31-base XOR chunk operation
- **Quality-aware refinement**: Smith-Waterman local alignment with Phred-scaled quality penalties
- **Combined ranking**: Formula balancing alignment score (85%) and k-mer rarity (15%)
- **Top-N rarest k-mer anchors**: Fallback to 2nd/3rd rarest k-mers for improved mapping rate
- **Paired-end support**: Full SAM specification compliance with proper FLAG handling
- **Parallel mapping**: Work-stealing scheduler using rayon for multi-core speedup

## Pipeline

1. **FM-index** (radix-sort suffix arrays) for efficient k-mer lookup
2. **Anchor-based k-mer filtering** (top-N rarest k-mer selection with fallback)
3. **2-bit XOR alignment** (~2.3 ns per 31-base chunk for exact/near-exact matches)
4. **Smith-Waterman refinement** for lower confidence scores (<0.9)
5. **Multi-genome ranking** with combined scoring formula

## Installation

### Prerequisites

- Rust toolchain (2021 edition)
- Python 3.x with Biopython (for read simulation)

### Build

```bash
cargo build --release
```

### Dependencies

```bash
pip install biopython
```

## Quick Start

Build an index from the included genomes and map simulated reads:

```bash
# 1. Build index (k=10, 4 threads)
cargo run --release -- build \
  -f Ecoli_K12_MG1655.fna \
  -f CP029198.1.fasta \
  -f Sac_cerevisiae_complete.fasta \
  -o multi3.bitpop \
  -k 10 -t 4

# 2. Map single-end reads (top_n=2 recommended for balance)
cargo run --release -- map \
  -i multi3.bitpop \
  -r simulated_ecoli_10k_new.fastq \
  -o ecoli_mappings.sam \
  -a hybrid -t 4 --top-n 2

# 3. Map paired-end reads
cargo run --release -- map \
  -i multi3.bitpop \
  --reads-1 SRR1060563_1.fastq \
  --reads-2 SRR1060563_2.fastq \
  -o salmonella_paired.sam \
  -a hybrid -t 4

# 4. View index statistics
cargo run --release -- stats -i multi3.bitpop
```

## Usage

### Build Index

```bash
./target/release/bit-pop build \
  -f genome1.fasta -f genome2.fasta -f genome3.fasta \
  -o index.bitpop \
  -k 10 \
  -t 4
```

### Map Reads

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
  -a hybrid \
  -t 4
```

### Show Index Statistics

```bash
./target/release/bit-pop stats -i index.bitpop
```

### Add Genomes to Existing Index

```bash
./target/release/bit-pop load \
  -i existing.bitpop \
  -f new_genome.fasta \
  -o updated.bitpop
```

### Map Reads Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --index` | Input index path | (required) |
| `-r, --reads` | Input FASTQ/FASTA file | (required for single-end) |
| `-1, --reads-1` | R1 FASTQ for paired-end | (required with -2) |
| `-2, --reads-2` | R2 FASTQ for paired-end | (required with -1) |
| `-o, --output` | Output SAM file | (required) |
| `-a, --align-mode` | Alignment mode: xor, sw, hybrid | xor |
| `-m, --min-score` | Minimum alignment score (0.0-1.0) | 0.7 |
| `-q, --min-quality` | Minimum Phred quality (0 = no filter) | 0 |
| `-t, --reads-threads` | Number of threads | 1 |
| `--top-n` | Number of rarest k-mers to try (1 = single anchor) | 1 |

### Align Modes

- `xor`: Fast 2-bit XOR alignment only
- `sw`: Smith-Waterman refinement for all reads
- `hybrid`: XOR first, SW only when confidence < 0.9

## Benchmark Results

### Setup

- **Genomes**: E. coli K-12 MG1655 (4.6 Mb), S. aureus (2.9 Mb), S. cerevisiae (12.2 Mb)
- **Reads**: 20,000 simulated reads (100 bp, 0.1% error rate, Q30-Q40)
- **k-mer size**: k=10
- **Hardware**: Standard desktop CPU (Windows, PowerShell)

### Results (k=10, top_n=1)

| Genome | Size | Mapped | Mapping Rate | Accuracy |
|--------|------|--------|--------------|----------|
| E. coli | 4.6 Mb | 9,755/10,000 | 97.6% | 99.9% |
| S. aureus | 2.9 Mb | 4,910/5,000 | 98.2% | 99.9% |
| S. cerevisiae | 12.2 Mb | 4,905/5,000 | 98.1% | 100.0% |
| **Total** | **19.7 Mb** | **19,570/20,000** | **97.9%** | **99.9%** |

### Results (k=10, top_n=3)

| Genome | Size | Mapped | Mapping Rate | Accuracy |
|--------|------|--------|--------------|----------|
| E. coli | 4.6 Mb | 9,924/10,000 | 99.2% | 99.9% |
| S. aureus | 2.9 Mb | 4,968/5,000 | 99.4% | 99.9% |
| S. cerevisiae | 12.2 Mb | 4,970/5,000 | 99.4% | 100.0% |
| **Total** | **19.7 Mb** | **19,862/20,000** | **99.3%** | **99.9%** |

**Performance trade-off**: top_n=3 is ~3x slower than top_n=1 (2.8s vs 0.9s for E. coli). Recommended: `--top-n 2` for balance between speed and accuracy.

**Throughput**: ~1,500 reads/second (top_n=1)

## Project Structure

```
├── src/                    # Rust source code
├── benches/                # Benchmark scripts
├── Cargo.toml              # Rust project configuration
├── paper.tex               # Academic paper
├── references.bib          # Bibliography
├── simulate_reads.py       # Read simulation tool
├── analyze_benchmark_new.ps1 # Benchmark analysis
└── README.md               # This file
```

### Data Files

**Genomes:**
- `Ecoli_K12_MG1655.fna` - E. coli K-12 MG1655 (4.6 Mb)
- `CP029198.1.fasta` - Staphylococcus aureus (2.9 Mb)
- `Sac_cerevisiae_complete.fasta` - S. cerevisiae S288C (12.2 Mb)

**Simulated Reads:**
- `simulated_ecoli_10k_new.fastq` - 10,000 E. coli reads
- `simulated_aureus_5k_new.fastq` - 5,000 S. aureus reads
- `simulated_cerevisiae_5k_new.fastq` - 5,000 S. cerevisiae reads

**Real Reads (examples):**
- `SRR1060563_1.fastq` / `_2.fastq` - Salmonella paired-end (~1.1M pairs)

## Limitations

- Proof-of-concept stage; not validated on large-scale real datasets
- No reverse complement support yet (Phase 1.2 planned)
- No clinical validation; academic research tool only
- Large index file sizes (~152MB for 19.7Mb genome, delta compression planned)

## Future Work

- **Phase 1**: Reverse complement support, entropy-adaptive k-mer size, seed-and-extend
- **Phase 2**: Parallel build, SA compression, streaming input, SIMD acceleration
- **Phase 3**: Progress reporting, CIGAR accuracy, quality filter improvements
- **Phase 4**: Integration tests, CAMI benchmark, performance regression tests
- **Phase 6**: NCBI E-utilities integration, local reference cache, batch download
- Expand benchmarks to 100+ genomes and eukaryotic genomes
- Direct comparison with Bowtie2, BWA-MEM on multi-genome tasks

## Availability

Source code available under the MIT License.

## Citation

```
Popović, M. (2025). Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification.
```

## License

MIT License
