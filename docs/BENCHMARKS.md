# Bit-Pop Benchmarks

## Summary

| Dataset | Genomes | Reads | Mapping Rate | Accuracy |
|---------|---------|-------|-------------|----------|
| Quick (3 genomes) | 3 | 10k | 99.3% | 99.9% |
| CAMI Low — k13 + EM | 61 | ~1M | 70.0% | 92.29% |
| CAMI Low — k12-k15 + EM | 61 | ~1M | 91.1% | 90.07% |
| CAMI Low — k13+k22 + EM | 61 | ~1M | 99.48% | 89.86% |
| PacBio HiFi (simulated) | 69 | 86k | 99.0% | 95.2% |

> **Species-level accuracy is ~100% across all benchmarks.** Misclassifications occur exclusively within clades (sibling strains), never between unrelated species.

---

## Quick Benchmark

**Dataset:** 3 genomes (E. coli K-12, S. aureus, S. cerevisiae), 19.7 Mb total, 20,000 simulated Illumina reads.

**Results:**
- Mapping rate: 99.3%
- Accuracy: 99.9%
- Speed: 0.9s per 10k reads
- Hardware: Intel i5, no GPU

---

## CAMI Low Complexity Benchmark

**Dataset:** 61 microbial genomes (28 parent strains, 13 evo_* strains, 20 Sample genomes), 157 Mb total, ~1,000,000 simulated Illumina short reads (150bp).

**Hardware:** Intel i5-14400, 16GB RAM, no GPU.

### Configuration Sweep

| Config | Mapping Rate | Strain Accuracy | Notes |
|--------|-------------|-----------------|-------|
| k=10, tn=4 | 86.2% | 86.52% | High coverage, lower precision |
| k=13, tn=4 | 70.0% | 91.91% | Best single-index accuracy |
| k=13, tn=4 + EM t=1.0 | 70.0% | **92.29%** | **Recommended: precision mode** |
| k=12+k=13+k=14+k=15, tn=2 + EM | 91.1% | 90.07% | **Recommended: balanced mode** |
| k=13+k=22, tn=2 + EM | 99.48% | 89.86% | **Recommended: clinical/coverage mode** |
| k=18+, tn=4 | 99.6% | 86.79% | Maximum mapping rate |

### Species vs Strain Accuracy

| Level | Accuracy |
|-------|----------|
| Species-level | ~100% |
| Strain-level (overall) | 89-92% |
| Strain-level (evo_* clades) | 60-92% |

**Key finding:** All misclassifications are intra-clade. A read from `evo_1049056.015` may be assigned to `evo_1049056.011` (sibling strain), but never to an unrelated species. This confirms correct species-level taxonomy in all cases.

### evo_* Strain Results

evo_* genomes are computationally simulated strains (sgEvolver) differing by 1-2 SNPs per Mb from their parent genome. With 150bp reads, P(read spans SNP) ≈ 0.015% — representing a fundamental information-theoretic limit, not an algorithmic issue.

| Genome | Accuracy |
|--------|----------|
| evo_1286_AP.008 | 92.18% |
| evo_1035930.011 | 91.11% |
| evo_1286_AP.026 | 68.95% |
| evo_1049056.031 | 72.96% |
| evo_1049056.039 | 47.23% |

### Runtime

- Index load: ~4s (1.4GB index)
- Mapping 1M reads: ~132s
- EM post-processing: ~0.2s
- Hardware: Intel i5-14400, 16 threads

---

## PacBio HiFi Benchmark

**Dataset:** 69 bacterial genomes (51 unique species), 285 Mb total, 86,248 simulated HiFi reads (8,000-20,000 bp).

**Simulation parameters (realistic error profile):**
- Base error rate: 0.1% (HiFi accuracy)
- Homopolymer error rate: 2% in homopolymer regions (≥4 same bases)
- Coverage variation: Gaussian ±50% (clamped 0.3x-2.0x)
- Chimeric reads: 1% (two fragments joined)
- Read length: 8,000-20,000 bp variable

**Results (k=70, no chunking):**
- Mapping rate: 99.0% (85,351/86,248)
- Accuracy: 95.2% (81,279/85,351)
- Runtime: ~8 minutes
- Hardware: Intel i5, no GPU

**Key finding:** For reads >1,000 bp, chunking is automatically disabled. Direct FM-index alignment outperforms chunk-based voting by +12.6pp accuracy for long reads. The 897 unmapped reads correspond to chimeric reads (~862 expected at 1% chimera rate) — correctly rejected rather than misassigned.

### Short vs Long Read Comparison

| Read Type | Length | Accuracy | Mapping Rate |
|-----------|--------|----------|-------------|
| Illumina short reads | 150 bp | 92.29% | 70-99% |
| PacBio HiFi | 8-20 kb | **95.2%** | **99.0%** |

Longer reads improve strain disambiguation because P(read spans SNP) scales with read length: 150 bp → 0.015%, 10,000 bp → 1.0%.

---

## Methodology Notes

- All benchmarks run on consumer hardware (Intel i5-14400, no GPU)
- Ground truth known for all benchmarks (simulated reads)
- CAMI Low dataset: standard community benchmark for metagenomic classifiers
- PacBio simulation: custom Python script with homopolymer, chimera, and coverage variation modeling (see `scripts/simulate_reads.py`)
- Accuracy metric: exact genome name match (strain-level); species-level computed separately
- EM temperature=1.0 recommended (temperature=0.1 over-concentrates probability mass)
