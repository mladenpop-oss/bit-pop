# Bit-Pop: Multi-Genome DNA Read Classification

[![CI](https://github.com/mladenpop-oss/bit-pop/actions/workflows/ci.yml/badge.svg)](https://github.com/mladenpop-oss/bit-pop/actions/workflows/ci.yml)
[![Docker](https://github.com/mladenpop-oss/bit-pop/actions/workflows/docker.yml/badge.svg)](https://github.com/mladenpop-oss/bit-pop/pkgs/container/bit-pop)
[![Tests](https://img.shields.io/badge/tests-312%2B%20unit%2C%205%20integration-blue)](https://github.com/mladenpop-oss/bit-pop)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20043593.svg)](https://doi.org/10.5281/zenodo.20043593)

> **Ultra-fast multi-genome DNA read classification.** Maps reads across dozens of genomes using a compact FM-index built in Rust with bit-level parallelism.

While existing aligners (Bowtie2, BWA, minimap2) map reads to a single reference genome, Bit-Pop identifies **which genome** in a collection best matches each read — ideal for metagenomic classification, outbreak detection, and clinical microbiology.

## Benchmarks

| Dataset | Genomes | Reads | Mapping Rate | Accuracy |
|---------|---------|-------|-------------|----------|
| Quick (3 genomes, 19.7 Mb) | 3 | 10k | **99.3%** | **99.9%** |
| CAMI Low (k13 + EM) | 61 | ~1M | 70.0% | **92.29% strain / ~100% species** |
| CAMI Low (k12-k15 consensus + EM) | 61 | ~1M | 91.1% | **90.07%** |
| CAMI Low (k13+k22 + EM) | 61 | ~1M | **99.48%** | 89.86% |
| PacBio HiFi (realistic simulation) | 69 | 86k | **99.0%** | **95.2%** |
| **Ebola strains (Nanopore R9.4, 15% errors)** | **3** | **10k/strain** | **68-72%** | **95-97%** |
| **Ebola strains (Nanopore R10.4, 5-8% errors)** | **3** | **10k/strain** | **80-81%** | **97%** |
| **Ebola mixed sample (3 strains, R9.4)** | **3** | **~3.3k/strain** | **56%** | **intra-clade only** |
| **Ebola + human contamination (R9.4)** | **4** | **5k each** | **59%** | **0% cross-contamination** |
| **Ebola outbreak detection (unknown variant)** | **2** | **27.5k** | **24%** | **76% unmapped → alert** |
| **Ebola 2026 DRC (3 refs, chunk-pct 3%)** | **3** | **27.5k** | **93.5%** | **Bundibugyo identified** |

> **Species-level accuracy is ~100% across all benchmarks.** Misclassifications occur only within clades (sibling strains), never between species. Human reads never map to viral genomes and vice versa.

See [docs/BENCHMARKS.md](docs/BENCHMARKS.md) for full benchmark details and methodology.

## Installation

```bash
# From source (requires Rust: https://rustup.rs)
git clone https://github.com/mladenpop-oss/bit-pop.git
cd bit-pop
cargo build --release

# Docker
docker pull ghcr.io/mladenpop-oss/bit-pop:latest

# Bioconda
conda install -c bioconda bit-pop
```

## Quick Start

```bash
# One-command workflow: build index + map reads
./target/release/bit-pop run \
  data/genomes/Ecoli_K12_MG1655.fna \
  data/reads/reads.fastq

# Paired-end
./target/release/bit-pop run \
  data/genomes/ \
  -1 R1.fastq -2 R2.fastq

# Download from NCBI and map
./target/release/bit-pop run \
  --ncbi "Escherichia coli" \
  reads.fastq

# Multi-k consensus (best accuracy/coverage balance)
bit-pop consensus \
  -i index_k13.bitpop index_k14.bitpop \
  -r reads.fastq -o output.sam

# EM post-processing (+0.38% accuracy)
bit-pop em -i mapped.sam -o mapped_em.sam --temperature 1.0
```

See [docs/USAGE.md](docs/USAGE.md) for full CLI reference.

## GUI

Desktop application (Tauri + Svelte) for users who prefer a graphical interface:

```bash
cd gui
npm install
npm run tauri build
```

**Build** — Create index from genome folder (k-mer size, threads, CAMI support)  
**Load/Run** — Map reads (single index or multi-index consensus mode, chunk-pct control)  
**Results** — Load SAM file, view statistics, filter and sort results  
**Help** — Step-by-step guide and GitHub link

## Outbreak Detection Mode

Bit-Pop works as both a classifier **and** an outbreak detector. The mapping rate itself is a diagnostic signal:

| Mapping Rate | Interpretation | Action |
|---|---|---|
| **>70%** | Known strain in index | Classify normally |
| **30-70%** | Related but divergent | Review top genomes, consider adding references |
| **<30%** | Unknown/novel variant | **ALERT** — possible new outbreak, escalate for sequencing |

**How it works:** When reads don't match any genome in the index, they appear as unmapped. A high unmapped rate signals the presence of a novel variant not in the reference database.

**Real-world example — Ebola 2026 DRC Outbreak:**
- **Index:** 2 reference strains (Zaire, Sudan) — no Bundibugyo
- **Reads:** 27,541 ONT reads from 6 new Bundibugyo outbreak strains
- **Result:** 24% mapped, **76% unmapped** → Unknown variant alert
- **Follow-up:** Adding old Bundibugyo reference → 86% mapped as Bundibugyo (93.5% with chunk-pct 3%)

**Best practice:** Use single index for clearest unmapped signal. Consensus mode may produce more false mappings.

## Key Features

- **FM-index + 2-bit XOR alignment** — ~2.3 ns per 31-base chunk
- **Multi-genome classification** — single index for all reference genomes
- **Multi-k consensus** — combine indexes at different k for accuracy/coverage trade-off
- **Dynamic chunking** — adaptive chunk size for long reads with high error rates (Nanopore)
- **EM post-processing** — soft-assignment refinement for ambiguous reads
- **PacBio HiFi support** — automatic long-read mode (no chunking for reads >1kb)
- **SAM/BAM output** — full spec compliance with CIGAR, NM, MAPQ tags
- **NCBI integration** — download and index genomes directly
- **Paired-end support** — discordant pair reconciliation
- **Taxonomic classification** — LCA algorithm (`bit-pop tax`)
- **312+ unit tests**, 5 integration tests, 17 benchmark groups
- **Android JNI** — native Android app with identical algorithm to CLI (99.95% accuracy parity)
- **Outbreak-ready** — 95-97% strain accuracy on Nanopore R9.4 (Ebola, 15% errors); 97% on R10.4
- **Human contamination safe** — 0% cross-mapping between human and viral reads

## Bit-Pop vs Kraken2

| | Bit-Pop | Kraken2 |
|---|---------|---------|
| Database size | **MB** (your genomes only) | 50-100 GB (full NCBI) |
| Build time | **17 seconds** | Hours to days |
| Custom genomes | **Trivial** | Requires full taxonomy dump |
| Offline use | ✅ | ✅ |
| "Unknown unknown" detection | ✅ High unmapped rate = novel variant | ✅ Full NCBI coverage |

**Use Bit-Pop when** you know which organisms to look for. **Use Kraken2 when** you need broad discovery against all of NCBI.

## Limitations

- **Strain-level resolution**: Genomes >99.9% identical share most k-mers — misclassification between sibling strains is expected and represents a fundamental information-theoretic limit. Species-level accuracy is ~100%.
- **Index size**: ~152 MB per 19.7 Mb genome
- **Max index size**: ~2 GB per index (libsais limit) — use `scripts/bitpop-workflow.py` for larger genomes

See [docs/ADVANCED.md](docs/ADVANCED.md) for experimental features and large genome support.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — Full CLI reference
- [docs/BENCHMARKS.md](docs/BENCHMARKS.md) — Benchmark methodology and results
- [docs/ADVANCED.md](docs/ADVANCED.md) — Experimental features, large genome support
- [docs/paper.pdf](docs/paper.pdf) — Academic paper (in preparation)

## Citation

```bibtex
@software{popovic_2026_bitpop,
  author = {Popović, Mladen},
  title = {Bit-Pop: Multi-Genome DNA Read Classification},
  year = {2026},
  doi = {10.5281/zenodo.20043593},
  url = {https://github.com/mladenpop-oss/bit-pop}
}
```

## License

MIT License
