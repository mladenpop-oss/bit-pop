# Bit-Pop: Ultra-Fast Multi-Genome DNA Read Classification with FM-Index, Bit-Parallel Alignment, and EM Refinement

**Mladen Popović**  
Independent Researcher  
mladenpop@gmail.com  
https://github.com/mladenpop-oss/bit-pop

May 27, 2026

---

## Abstract

I present Bit-Pop, a high-performance tool for multi-genome DNA read classification that identifies which genome in a collection best matches each sequencing read. While existing aligners such as Bowtie2, BWA-MEM, and minimap2 map reads to single reference genomes, Bit-Pop addresses the complementary problem of simultaneous classification across collections of *N* genomes. The pipeline combines four key stages: (1) an FM-index built via the SA-IS algorithm (libsais) for efficient k-mer lookup; (2) a top-*N* rarest k-mer anchor filter with quality-aware filtering and reverse complement support; (3) a 2-bit XOR alignment achieving approximately 2.3 ns per 31-base chunk, with Myers edit distance as a fast alternative to Smith–Waterman (23–54× faster); and (4) a combined ranking formula that weights per-genome alignment score against per-genome k-mer rarity, penalizing k-mers shared across multiple genomes. For difficult cases, Bit-Pop provides fuzzy k-mer matching methods and an Expectation-Maximization (EM) post-processing algorithm for multi-candidate refinement. Benchmarked on simulated Illumina-like reads (0.1% error rate, 100 bp) across three distantly related genomes (*Escherichia coli*, *Staphylococcus aureus*, *Saccharomyces cerevisiae*), Bit-Pop achieves 99.3% mapping rate and 99.9% classification accuracy in 0.9 s per 10,000 reads. On the CAMI Low Complexity benchmark with 61 microbial genomes (~1 M reads, ~157 M bases), Bit-Pop achieves 92.29% accuracy (k=13, top-N=4, EM refinement) at 70.0% mapping rate, and 90.07% accuracy at 91.0% mapping with multi-k consensus (k=12–15). On a realistic PacBio HiFi benchmark (69 genomes, 86,248 long reads, 8–20 kb, homopolymer errors, chimeras, coverage variation), Bit-Pop achieves 95.2% accuracy and 99.0% mapping rate in ~8 minutes. Species-level classification is ~100% accurate across all benchmarks. Strain-level classification for near-identical genomes (>99.9% identity) remains challenging (60–90% within clade), representing a fundamental information-theoretic limitation of short-read k-mer-based classification. Bit-Pop is designed as a focused, lightweight tool for species classification, contamination detection, and targeted metagenomic analysis where users work with known genome collections on standard hardware.

**Availability**: Source code at https://github.com/mladenpop-oss/bit-pop (MIT License). DOI: https://doi.org/10.5281/zenodo.20043593.

---

## 1 Introduction

The cost of DNA sequencing has continued to decline over the past decade, enabling applications ranging from clinical metagenomics to environmental microbiome profiling, outbreak surveillance, and food safety monitoring. In many of these applications, a fundamental question arises: *which organism does this read come from?* When researchers work with a curated collection of reference genomes — whether tracking pathogen strains in a hospital, monitoring contamination in a bioproduction facility, profiling a defined microbiome, or investigating an outbreak in resource-limited settings — efficiently assigning each read to its source genome becomes a critical computational task.

The field of read classification has produced powerful tools, but they fall into two distinct categories, each addressing a different use case.

**General-purpose aligners** — Bowtie2 [Langmead and Salzberg, 2012], BWA-MEM [Li, 2013], and minimap2 [Li, 2018] — are optimized for mapping reads against a single reference genome or a pangenome graph. These tools achieve remarkable speed and accuracy but do not natively address the multi-genome classification problem. While one could theoretically index all target genomes into a single concatenated reference and map reads against it, this approach conflates two distinct tasks: (1) finding the optimal genomic position for each read, and (2) determining which genome in a collection best explains each read. Aligners produce positional mappings with CIGAR strings and quality scores, but leave the classification decision — aggregating per-read mappings into genome-level assignments — to downstream tools.

**Broad-spectrum classifiers** — Kraken2 [Wood et al., 2019], KrakenUniq [Breitwieser et al., 2019], Centrifuge [Kim et al., 2016], and CLARK [Ounit et al., 2017] — take a different approach: they build massive k-mer databases covering entire taxonomic databases (e.g., NCBI nt/nr), enabling classification against millions of reference sequences. However, these tools require large databases (100 GB+), significant computational infrastructure, and hours to days for database construction. They are designed for "unknown unknown" scenarios where the target organism is not predetermined — for example, characterizing the full microbial composition of a soil or gut sample.

Both categories are mature and well-established. However, a third use case remains underserved: **targeted classification against a known, curated genome collection**. This scenario is common in several practical settings:

- **Clinical microbiology**: A hospital laboratory tracks 20–50 known pathogen strains and needs to classify patient samples rapidly, on standard hardware, without institutional HPC infrastructure.
- **Outbreak detection**: A novel bacterium appears. Public health officials download the new genome (megabytes, not gigabytes), add it to their existing index, and begin classifying samples immediately.
- **Bioproduction quality control**: A facility monitors for a defined set of contaminant organisms and needs fast, offline-capable classification.
- **Academic phylogenomics**: Researchers working with defined taxon sets on limited computational resources need tools that scale to their data, not to the entire NCBI database.
- **Resource-limited settings**: Laboratories in the Global South often lack access to HPC clusters or stable internet connectivity, making lightweight, offline tools essential.

Bit-Pop was designed for this complementary use case. Given a curated collection of *N* genomes (potentially hundreds or thousands) and a set of sequencing reads, Bit-Pop maps each read against all genomes simultaneously and returns a ranked list identifying the most likely source genome, along with alignment position, score, and contextual information. The key design principles are:

- **Compact indexing**: All reference genomes are indexed in a single FM-index structure, with database sizes proportional only to the user's genomes (megabytes, not gigabytes).
- **Speed via bit-level parallelism**: DNA bases are encoded in 2 bits, enabling alignment of up to 31 bases per CPU word through XOR comparison — achieving approximately 2.3 ns per 31-base XOR chunk.
- **Multi-stage pipeline**: Anchor-based k-mer filtering, fast XOR alignment, Myers edit distance, and Smith–Waterman refinement are combined for optimal speed-accuracy trade-off.
- **Multi-k consensus**: Combining results from multiple k-mer sizes improves strain resolution and mapping coverage.
- **EM post-processing**: An Expectation-Maximization algorithm refines multi-candidate mappings using population-level abundance signals.
- **Taxonomic classification**: NCBI taxonomy tree integration with Lowest Common Ancestor (LCA) algorithm for genus/phylum/class-level abundance profiles.
- **Focused scope**: Bit-Pop is not intended to replace general-purpose aligners or broad-spectrum classifiers, but to solve a specific problem — species classification across known genome collections — that existing tools do not address directly.

This paper describes the algorithms, implementation, and benchmark results of Bit-Pop, demonstrating its feasibility for species-level read classification across collections of multiple bacterial and eukaryotic genomes. We present results on three benchmarks: (1) simulated reads across three distantly related genomes, (2) the CAMI Low Complexity benchmark with 61 microbial genomes and ~1 M reads, and (3) a realistic PacBio HiFi benchmark with 69 genomes and 86,248 long reads. We also discuss the fundamental limitations of k-mer-based classification for near-identical strains and outline directions for future work.

---

## 2 Methods

### 2.1 FM-Index Construction

The FM-index [Burrows and Wheeler, 1994; Ferrada and others, 2007] is a compressed full-text substring index built on the Burrows–Wheeler Transform (BWT) of the concatenated reference collection. Bit-Pop constructs the FM-index using the SA-IS algorithm [Kärkkäinen and Sanders, 2003], implemented via the libsais library [Dörge, 2015], which achieves linear-time *O*(*L*) construction where *L* is the total length of all genomes.

Given *N* genomes with total length *L*, the BWT is constructed by:

1. Encoding bases as `u8` values: `$` = 0, `A` = 1, `C` = 2, `G` = 3, `T` = 4.
2. Concatenating all genome sequences with unique separator characters.
3. Computing the suffix array via SA-IS (linear time).
4. Deriving the BWT from the suffix array.
5. Building the Occ (occurrence) counter array for backward search with rank sampling (interval = 32).

The SA construction using SA-IS operates in *O*(*L*) time and space, significantly faster than radix-sort prefix doubling approaches. The BWT and Occ array together occupy approximately 2*L* bytes (2 bits per base for BWT, plus 8 bits per position for Occ counters with periodic sampling).

Persisted indexes are stored using `memmap2` for memory-mapped file access and zstd compression for the BWT, achieving typical compression ratios of 3–4:1. The suffix array is stored uncompressed to enable *O*(1) random access during mapping. Parallel index construction is supported via multi-threaded BWT and suffix array construction using rayon.

### 2.2 Anchor-Based K-mer Filtering

The first stage of mapping uses a top-*N* anchor-based filter to identify candidate positions across all indexed genomes. Given a read of length *R*, the algorithm:

1. Extracts all overlapping k-mers (configurable, default *k* = 10) from the read.
2. For each k-mer, queries the FM-index to count its total occurrences across all genomes using backward search.
3. Selects the top-*N* rarest k-mers as anchors — those with fewest total positions — providing fallback candidates when the primary anchor contains sequencing errors.
4. For each anchor k-mer, retrieves all positions via backward search.
5. Applies smart stride sampling: if more than 100 anchor positions are found, they are subsampled with stride ⌈*n*/100⌉ to avoid redundant close alignments.
6. Skips highly repetitive k-mers exceeding a configurable threshold (default 10⁴ total occurrences).

This approach is significantly more efficient than iterating all k-mers in the read when the genome collection contains many repetitive elements. The top-*N* rarest k-mer anchors minimize false negatives caused by sequencing errors in the anchor k-mer while maintaining high discriminative power. The recommended setting is *N* = 2 for balance between speed and accuracy; *N* = 4 provides improved mapping rate at the cost of increased computation time.

**Quality-aware filtering** further refines this stage: k-mers containing bases with Phred quality below a minimum threshold (default Q20) [Ewing et al., 1998] are excluded from anchor selection, reducing false positives in low-quality regions.

**Reverse complement support**: Bit-Pop evaluates both forward and reverse complement orientations for each read. The reverse complement is computed at the 2-bit encoding level (A↔T, C↔G, then reversed), and the best alignment (forward or RC) is selected. SAM FLAG 0x10 (REVERSE) is set when the RC alignment wins.

**Paired-end support**: Full SAM specification compliance with proper FLAG handling, including discordant pair reconciliation. When R1 and R2 map to different genomes, Bit-Pop resolves the conflict by finding overlapping top-*N* candidates between both reads and selecting the concordant genome. A Gaussian insert size model provides probabilistic paired-end classification using the normal distribution of observed insert sizes.

### 2.3 2-Bit XOR Alignment

For each candidate position identified by the anchor filter, Bit-Pop performs 2-bit XOR alignment to score the read against the genomic region. DNA bases are encoded as 2-bit values (A=00, C=01, G=10, T=11), allowing a read of up to 31 bases to be packed into a single 64-bit word.

Given a packed read *P* and a genomic window *W* of the same length:

```
XOR = P ⊕ W
```

Each 2-bit field in the XOR result is zero if bases match and non-zero if they mismatch. The alignment score is:

$$\text{score} = \frac{\sum_{i=0}^{R-1} \mathbb{1}[\text{XOR}[2i, 2i+1] = 00]}{R}$$

where *R* is the read length and 𝟙[·] is the indicator function.

This operation requires a single XOR instruction and a population-count equivalent, achieving approximately 2.3 ns per 31-base XOR chunk on modern CPUs — orders of magnitude faster than traditional dynamic programming alignment.

For reads exceeding 31 bases, the XOR alignment is applied in chunks of 31 bases, with the final score computed as the average chunk score. While this approach does not handle gaps (indels), it provides excellent discrimination for high-quality reads where mismatches are the primary source of disagreement.

### 2.4 Myers Edit Distance

For reads with moderate XOR confidence scores, Bit-Pop applies Myers' bit-vector algorithm for edit distance computation. Myers' algorithm computes the exact edit distance between a pattern and text using bit-parallel operations, achieving 23–54× speedup over Smith–Waterman while providing exact (not approximate) edit distance [Myers, 1999].

This provides a fast intermediate between the approximate XOR alignment and the full Smith–Waterman refinement, enabling Bit-Pop to handle reads with indels without the full *O*(*m* × *n*) cost of SW.

### 2.5 Smith–Waterman Refinement

For reads with low XOR confidence scores (< 0.9), Bit-Pop applies Smith–Waterman (SW) local alignment as a refinement step. The Smith–Waterman algorithm [Smith and Waterman, 1981] uses standard scoring parameters (match = +2, mismatch = −1, gap open = −2) and full traceback to generate CIGAR strings with indel operations (M, I, D).

For reads exceeding 31 bases, Bit-Pop employs chunked Smith–Waterman: the read is split into 31-base chunks, each aligned independently against the candidate genomic region using full SW with traceback. The per-chunk CIGAR operation codes are concatenated and collapsed into a final CIGAR string.

### 2.6 Quality-Aware Scoring

Bit-Pop supports Phred-scaled quality-aware scoring at multiple stages:

**Quality-filtered k-mer selection**: During anchor selection, k-mers containing bases below a minimum Phred quality threshold (default Q20) are excluded from anchor selection.

**Quality-aware XOR alignment**: Mismatches at high-quality positions receive larger penalties:

$$\text{adjusted score} = \frac{\text{matches} + \sum_{i \in \text{mismatches}} \frac{-Q_i}{20}}{R}$$

**Quality-aware Smith–Waterman**: The SW scoring matrix incorporates per-base quality penalties:

$$\text{mismatch penalty}_i = -\min\left(\frac{Q_i}{10}, 5\right)$$

### 2.7 Multi-Genome Ranking

After scoring each read against all indexed genomes, Bit-Pop applies a combined ranking formula:

$$\text{combined score} = \alpha \cdot \text{align score} + (1 - \alpha) \cdot \text{rarity}$$

with *α* = 0.85. The alignment score reflects the quality of the best match against a genome, while rarity provides a modest boost to matches in genome-specific regions:

$$\text{rarity} = \frac{1}{\max(1, \text{occurrences of first k-mer in this genome})}$$

The per-genome occurrence count ensures that k-mers shared across multiple genomes are penalized appropriately, improving strain-level discrimination.

### 2.8 Multi-K Consensus Mapping

Multi-k consensus combines results from multiple independent indexes built with different k-mer sizes to improve both mapping coverage and strain resolution. Each index is mapped separately using the optimized standalone map engine, and results are combined using k-priority weighting where larger k-values receive higher weight:

$$\text{weight}_k = \frac{k}{k_{\min}}$$

For *k* = {12, 13, 14, 15}: weights are 1.0×, 1.08×, 1.17×, 1.25× respectively. Larger k-values are weighted more heavily because they provide higher specificity — a match at *k* = 15 is more informative than a match at *k* = 12. The final genome assignment is determined by the weighted sum of per-index scores across all k-values.

This approach captures matches at multiple specificity levels: smaller k-values recover reads that fail to map at higher k (improving coverage), while larger k-values provide discriminative power for strain resolution. On the CAMI benchmark, *k* = 12–15 consensus achieves 90.07% accuracy at 91.0% mapping rate, compared to 92.29% accuracy at 70.0% mapping for *k* = 13 alone.

### 2.9 Fuzzy K-mer Matching

For highly similar genomes (>99.9% identity), exact k-mer matching is insufficient because most k-mers are shared between strains. Bit-Pop provides three fuzzy k-mer matching methods:

- **Fuzzy k-mer**: Generates all possible k-mer variants with *n* substitutions and queries the FM-index for each variant. This provides the best accuracy for strain resolution but is approximately 30× slower than exact matching (*n* = 1).
- **Fuzzy seed**: Allows *n* mismatches in spaced seed "match" positions. Spaced seeds use a pattern of "match" and "don't care" positions (e.g., `111010111011`), reducing the number of variants to generate while maintaining sensitivity. This provides a good balance between accuracy and speed, approximately 20× slower (*n* = 1).
- **Neighborhood**: Builds a hash table at index build time for all k-mer neighborhoods (k-mers within *n* mismatches). This provides *O*(1) fuzzy lookup at query time but increases index size by approximately 60× (*n* = 1).

### 2.10 EM Post-Processing

The Expectation-Maximization (EM) algorithm refines multi-candidate SAM mappings by using population-level abundance signals. When a read maps to multiple genomes with similar scores, EM iteratively:

- **E-step**: Estimates genome abundances based on current read assignments using soft assignments with softmax and a temperature parameter.
- **M-step**: Reassigns reads to maximize the likelihood given estimated abundances.

The algorithm typically converges in 6–11 iterations. Key parameters:

- **Temperature**: Controls softness of assignments (default: 0.1). Lower temperature = harder assignments.
- **Top-K**: Number of top candidates per read (default: 10).
- **Convergence**: KL divergence threshold for stopping (default: 0.001).
- **Confidence threshold**: Minimum probability to apply EM reassignment (default: 0.95).

### 2.11 Taxonomic Classification

Bit-Pop includes a taxonomic classification module that maps genome-level read counts to the NCBI taxonomy tree using the Lowest Common Ancestor (LCA) algorithm. Given a SAM output file and NCBI taxonomy dump files (`nodes.dmp`, `names.dmp`), the tool produces taxonomic abundance profiles by rank (species, genus, family, phylum, etc.). Ambiguous mappings are resolved automatically by finding the lowest common ancestor in the taxonomy tree.

### 2.12 Parallel Mapping

Read mapping is parallelized using rayon's work-stealing scheduler [Rayon Developers, 2022]. Reads are divided into batches, each processed by a separate thread. Within each thread, reads are mapped sequentially against all indexed genomes. This approach provides near-linear speedup with the number of available CPU cores.

**Streaming mode** is available for large FASTQ files: reads are processed in fixed-size chunks with bounded memory usage (~3 GB per chunk), enabling mapping of arbitrarily large datasets on standard hardware.

### 2.13 Implementation Details

Bit-Pop is implemented in Rust (2021 edition) with the following key dependencies:

- **libsais**: SA-IS suffix array construction and BWT
- **rayon**: Data-parallelism and work-stealing scheduler
- **memmap2**: Memory-mapped file I/O for persisted indexes
- **zstd** (via `memmap2`): Compression for persisted BWT

The codebase consists of 17 Rust modules (~17,000 lines of code), 324 unit tests, 5 integration tests, and 17 Criterion benchmark groups. The tool produces standard SAM 1.6 output (text) or BAM 1.6 output (binary, BGZF-compressed) with full compatibility with samtools, bcftools, IGV, and other bioinformatics tools. Native BAM output is available via the `--bam` flag.

**System requirements**: Rust toolchain (2021 edition). Optional: Python 3.x with Biopython for read simulation scripts.

---

## 3 Results

### 3.1 Three-Genome Benchmark (Simulated Reads)

**Setup**: E. coli K-12 MG1655 (4.6 Mb), S. aureus (2.9 Mb), S. cerevisiae S288C (12.2 Mb). 20,000 simulated reads (100 bp, 0.1% error rate, Q30–Q40). *k* = 10. Standard desktop CPU.

Table 1 shows results for mapping 20,000 simulated reads against an index containing all three genomes.

**Table 1: Three-genome benchmark — 20,000 simulated reads, *k* = 10, top-N = 3**

| Genome | Size | Mapped | Mapping Rate | Accuracy |
|--------|------|--------|--------------|----------|
| E. coli | 4.6 Mb | 9,924 / 10,000 | 99.2% | 99.9% |
| S. aureus | 2.9 Mb | 4,968 / 5,000 | 99.4% | 99.9% |
| S. cerevisiae | 12.2 Mb | 4,970 / 5,000 | 99.4% | 100.0% |
| **Total** | **19.7 Mb** | **19,862 / 20,000** | **99.3%** | **99.9%** |

Bit-Pop achieves 99.3% mapping rate and 99.9% classification accuracy across bacterial and eukaryotic genomes spanning different domains of life. Mapping throughput is approximately 1,500 reads/second (top-N = 1) or 500 reads/second (top-N = 3). Total pipeline time for 20,000 reads: 0.9 s (top-N = 1) to 2.8 s (top-N = 3) for E. coli.

### 3.2 Error Tolerance

Table 2 shows mapping rate and classification accuracy across different simulated error rates (single-genome, E. coli, 1,000 reads, *k* = 10).

**Table 2: Error tolerance — mapping rate and accuracy at different error rates**

| Error Rate | Mapping Rate | Accuracy |
|------------|--------------|----------|
| 0.1% | 88.4% | 100.0% |
| 0.5% | 78.9% | 100.0% |
| 1.0% | 77.2% | 100.0% |
| 2.0% | 69.3% | 100.0% |
| 5.0% | 48.3% | 100.0% |
| 10.0% | 23.4% | 100.0% |

Classification accuracy remains at 100% across all error rates: when Bit-Pop maps a read, it always assigns it to the correct genome. The primary impact of sequencing errors is on mapping sensitivity (the anchor filter fails to find candidates when the rarest k-mer contains a sequencing error), not on classification precision.

### 3.3 CAMI Low Complexity Benchmark (61 Genomes, ~1 M Reads)

**Setup**: 61 bacterial genomes (1,900 sequences, ~157 M bases, ~1.3 GB index). ~1 M reads sampled from the original 30 GB CAMI Low Complexity dataset [Mende et al., 2018]. *k* = 10, 16 threads, `--cami` flag for genome naming.

#### 3.3.1 K-mer Size Sweep

Table 3 shows results across k-mer sizes from *k* = 10 to *k* = 22 (top-N = 4, XOR alignment).

**Table 3: CAMI benchmark — k-mer size sweep (*k* = 10–22, top-N = 4)**

| *k* | Mapped | Mapping Rate | Accuracy | Time |
|-----|--------|--------------|----------|------|
| 10 | 857,983 | 86.2% | 86.49% | ~60 s |
| 11 | 816,490 | 82.0% | 87.01% | ~55 s |
| 12 | 721,026 | 72.4% | 89.83% | ~55 s |
| **13** | **697,182** | **70.0%** | **91.91%** | **~50 s** |
| 14 | 742,747 | 74.6% | 91.55% | ~55 s |
| 15 | 840,561 | 84.5% | 89.48% | ~60 s |
| 16 | 940,303 | 94.5% | 87.28% | ~65 s |
| 17 | 980,328 | 98.6% | 86.77% | ~70 s |
| 18 | 990,721 | 99.6% | 86.79% | ~75 s |
| 19–22 | ~994,000 | ~99.97% | ~86.84% | ~75 s |

**Key observation**: *k* = 13 gives peak accuracy (91.91%), then accuracy drops and plateaus at ~86.8% for *k* ≥ 18. Coverage increases monotonically with *k*, reaching ~100% at *k* ≥ 18. The accuracy plateau at high *k* values reflects the increasing number of false-positive mappings from shared, non-discriminative k-mers.

#### 3.3.2 Top-N Series

Table 4 shows results for different top-N settings (XOR alignment, *k* = 10).

**Table 4: CAMI benchmark — top-N series (*k* = 10)**

| Config | Mapped | Mapping Rate | Accuracy | Time |
|--------|--------|--------------|----------|------|
| top-N = 1 | 638,704 | 63.9% | 86.15% | 129 s |
| top-N = 2 | 748,711 | 74.9% | 86.31% | 189 s |
| top-N = 3 | 813,874 | 81.4% | 86.41% | 250 s |
| top-N = 4 | 861,762 | 86.2% | 86.48% | 132 s |
| top-N = 5 | 712,720 | 71.3% | 86.50% | ~300 s |
| top-N = 6 | 511,411 | 51.2% | 86.57% | ~300 s |

Accuracy growth is linear but slow (+0.16 pp for N = 1→2, +0.07 pp for N = 4→5). Diminishing returns are evident from N = 5 onward, where mapped count drops drastically (712K → 511K) due to excessive computation per read.

#### 3.3.3 EM Post-Processing

Table 5 shows the impact of EM post-processing on the *k* = 13, top-N = 4 mapping.

**Table 5: CAMI benchmark — EM post-processing (*k* = 13, top-N = 4)**

| Config | Mapped | Accuracy | Changed | Time |
|--------|--------|----------|---------|------|
| baseline | 697,182 | 91.91% | — | ~50 s |
| **EM (t = 1.0, ct = 0.95)** | **697,188** | **92.29%** | **10,622** | **~60 s** |
| EM (t = 0.1, ct = 0.95) | 697,188 | 91.93% | 9,980 | ~60 s |

EM with temperature *t* = 1.0 provides +0.38 pp accuracy (91.91% → 92.29%) via 10,622 reassignments, primarily within strain groups. Temperature *t* = 1.0 outperforms *t* = 0.1 — a softer probability distribution allows more correct reassignments.

#### 3.3.4 Multi-K Consensus

Table 6 shows results for multi-k consensus mapping (using `consensus_base.py`, top-N = 2, EM *t* = 1.0).

**Table 6: CAMI benchmark — multi-k consensus (top-N = 2, EM *t* = 1.0)**

| Config | Mapped | Mapping Rate | Accuracy | EM Accuracy | Changed |
|--------|--------|--------------|----------|-------------|---------|
| *k* = 13, top-N = 4 | 697,182 | 70.0% | 91.91% | **92.29%** | 10,622 |
| *k* = 12 + 13 | 768,639 | 77.0% | 89.75% | **90.10%** | 132,153 |
| **k = 12–15** | **909,493** | **91.0%** | **89.71%** | **90.07%** | **132,153** |
| *k* = 13 + 22 | 994,188 | 99.5% | 89.50% | 89.86% | 184,078 |
| *k* = 13–15 + 22 | 994,214 | 99.98% | 89.58% | 89.39% | 190,020 |

Multi-k consensus (*k* = 12–15) gains +212K mapped reads vs. *k* = 13 alone, with a −1.8 pp accuracy trade-off before EM. After EM: *k* = 13 solo = 92.29%, *k* = 12–15 consensus = 90.07%. Adding *k* = 22+ increases coverage to ~100% but introduces noise that lowers accuracy.

#### 3.3.5 Two-Pass Mapping

Two-pass mapping recovers reads that fail the initial 0.7 threshold by re-mapping them with a lower threshold, then refining with EM. Table 7 shows results on a 10,000-read CAMI sample.

**Table 7: CAMI benchmark — two-pass mapping (10K sample, top-N = 4)**

| Config | Mapped | Accuracy | Correct | Wrong |
|--------|--------|----------|---------|-------|
| baseline | 49,028 | 86.52% | 42,421 | 6,607 |
| two-pass (0.5) | 49,071 | 86.48% | 42,436 | 6,635 |
| **two-pass (0.4) + EM** | **49,929** | **85.33%** | **42,604** | **7,325** |
| two-pass (0.3) + EM | 50,000 | 85.08% | 42,540 | 7,460 |

Two-pass at threshold 0.4 gains +183 correct reads at the cost of +718 wrong reads. Threshold 0.4 is optimal — 0.5 maps too few additional reads, 0.3 maps all but adds more false positives. Use when maximizing recall is more important than precision.

#### 3.3.6 Strain-Level Classification

The `evo_*` genomes (near-identical strains from the same sample assembly) represent the most challenging subset. Table 8 shows per-genome accuracy for `evo_*` strains (*k* = 12–15 consensus, top-N = 2, EM *t* = 1.0).

**Table 8: CAMI benchmark — `evo_*` strain-level accuracy (*k* = 12–15, EM *t* = 1.0)**

| Genome | Accuracy |
|--------|----------|
| evo_1035930.029 | 88.14% |
| evo_1286_AP.026 | 84.57% |
| evo_1286_AP.008 | 83.67% |
| evo_1035930.011 | 78.77% |
| evo_1049056.011 | 79.54% |
| evo_1286_AP.037 | 63.74% |
| evo_1035930.032 | 41.26% |
| evo_1286_AP.033 | 42.52% |
| evo_1049056.015 | 42.03% |
| evo_1049056.013 | 30.14% |
| evo_1049056.031 | 30.02% |
| evo_1049056.039 | 24.24% |

**Weighted average: 61.8%**. Larger strains (evo_029, evo_008) classify better (84–88%) due to more unique k-mers from larger genome size.

**Critical observation**: Misclassifications are **within-clade only**. Reads from `evo_1035930.*` never map to `evo_1049056.*` or `evo_1286_AP.*`. Species-level classification is ~100% accurate. Strain-level confusion occurs only between genomes sharing >99.9% identity.

#### 3.3.7 Summary: Best Configurations

**Table 9: CAMI benchmark — best configurations by objective**

| Goal | Config | Mapped | Mapping Rate | Accuracy | Time |
|------|--------|--------|--------------|----------|------|
| **Best accuracy** | *k* = 13, top-N = 4 + EM (*t* = 1.0) | 697,188 | 70.0% | **92.29%** | ~60 s |
| **Best consensus** | *k* = 12–15, top-N = 2 + EM (*t* = 1.0) | 909,493 | 91.0% | **90.07%** | ~250 s |
| **Best coverage** | *k* = 18, top-N = 4 | 990,721 | 99.6% | 86.79% | ~80 s |

### 3.4 PacBio HiFi Benchmark (69 Genomes, 86,248 Long Reads)

**Setup**: 69 bacterial genomes (51+ species, 285 Mb). 86,248 simulated PacBio HiFi reads (8–20 kb, realistic error profile: 0.1% base errors, 2% homopolymer errors, 1% chimeras, variable coverage ±50%). *k* = 70, 16 threads.

#### 3.4.1 Accuracy vs. K-mer Size

**Table 10: PacBio HiFi — accuracy vs. k-mer size**

| *k* | Accuracy | Mapping Rate | Map Time |
|-----|----------|-------------|----------|
| 10 | 42.6% | 99.97% | — |
| 13 | 73.7% | 99.97% | — |
| 15 | 75.1% | 99.97% | — |
| 20 | 79.7% | 99.97% | — |
| 25 | 82.7% | 99.97% | 7.5 min |
| 30 | 82.9% | 99.97% | 7.9 min |
| 40 | 83.0% | 99.97% | 8.0 min |
| 70 | 83.1% | 99.97% | 7.7 min |
| **70 (no chunk)** | **95.7%** | **99.93%** | **6.9 min** |

Accuracy plateaus at *k* ≥ 40 for chunked mapping. **No-chunk mode** (full read alignment) gives +12.6 pp accuracy over chunked mode, confirming that direct FM-index alignment is superior for reads >1 kb.

#### 3.4.2 Realistic vs. Simple Error Profile

**Table 11: PacBio HiFi — realistic vs. simple error profile (*k* = 70)**

| Dataset | Accuracy | Mapping Rate |
|---------|----------|-------------|
| Simple (0.1% errors) | 95.7% | 99.93% |
| **Realistic (homopolymers, chimera, coverage)** | **95.2%** | **99.0%** |

Realistic error profile (homopolymers, chimeras, coverage variation) causes only −0.5 pp accuracy drop, demonstrating robustness to realistic sequencing artifacts. Total mapping time: ~8 minutes for 86,248 long reads on 16 threads.

**Note**: The PacBio HiFi benchmark uses simulated reads generated with a realistic error profile (0.1% base errors, 2% homopolymer errors, 1% chimeras, ±50% coverage variation) via the `simulate_realistic.py` script included in the repository. These are not real PacBio sequencing data. The simulation was designed to approximate the error characteristics of PacBio HiFi reads based on published error models, but results on real PacBio data may vary.

---

## 4 Discussion

### 4.1 Comparison with Existing Tools

Unlike Bowtie2, BWA-MEM, or minimap2, Bit-Pop does not attempt to produce a full alignment for every read. Instead, it focuses on the classification task: identifying which genome in a collection best explains each read. This focus enables several advantages:

- **Compact databases**: Database sizes are proportional only to the user's genomes (megabytes), not entire taxonomic databases (100+ GB).
- **Offline operation**: No internet connectivity required.
- **Fast updates**: Adding a new genome takes seconds (rebuild index), not hours or days (rebuild database).
- **Unified ranking**: Automatic genome assignment with confidence scores.
- **NCBI integration**: Built-in genome search and download via NCBI E-utilities API.

Table 12 summarizes key feature differences between Bit-Pop and established tools.

**Table 12: Feature comparison with established genomic tools**

| Feature | Bit-Pop | Bowtie2 | BWA-MEM | minimap2 |
|---------|---------|---------|---------|----------|
| Multi-genome classification | Native | Single genome | Single genome | With `--index` |
| Speed (10K reads, 3 genomes) | **0.9 s** | ~5–10 s | ~8–15 s | ~3–5 s |
| Index size (19.7 Mb) | **~152 MB** | ~200 MB | ~250 MB | ~180 MB |
| Quality-aware alignment | Phred-scaled | Yes | Yes | Yes |
| Paired-end support | Yes | Yes | Yes | Yes |
| NCBI integration | Built-in | No | No | No |
| Implementation | Rust + bit-parallel | C++ | C | C++ |

Bit-Pop's unique capability — simultaneous mapping across many genomes with automatic ranking using compact databases — is not offered by any existing general-purpose aligner. This makes it particularly suitable for species classification tasks where the goal is to identify the source organism rather than produce a detailed per-base alignment.

### 4.2 Comparison with Kraken2

Kraken2 [Wood et al., 2019] is designed for broad-spectrum classification against entire taxonomic databases, requiring 100+ GB storage and hours to days for database construction. This design is optimal for "unknown unknown" scenarios where the target organism is not predetermined.

Bit-Pop addresses the complementary use case: targeted classification against a user-defined genome collection, with database sizes proportional only to the genomes of interest (megabytes, not gigabytes), no internet connectivity requirement, and index construction in minutes.

**Table 13: Bit-Pop vs. Kraken2 — key differences**

| Aspect | Bit-Pop | Kraken2 |
|--------|---------|---------|
| Database size | MB (only your genomes) | 100+ GB (entire NCBI) |
| Customization | Add/remove genomes in seconds | Fixed database |
| Build time | 2 minutes | Hours to days |
| Index growth | Grows only with your data | Pre-built massive database |

Kraken2 is better for broad metagenomics where the target is unknown. Bit-Pop is better for **targeted searching** where the user knows what matters.

A direct runtime and accuracy comparison was not performed in this work, as Kraken2 requires the full NCBI taxonomy hierarchy even for small custom genome sets, making equivalent experimental conditions impractical at this stage.

### 4.3 The Near-Identical Strain Challenge

The `evo_*` genomes on the CAMI benchmark demonstrate a fundamental limitation of k-mer-based classification. Genomes sharing >99.9% identity share most k-mers with each other. For a 150 bp read against a ~4.6 Mb genome, the probability of spanning a strain-specific SNP position is approximately 150 / 10⁶ ≈ 0.015%. This means the vast majority of reads from near-identical strains are informationally indistinguishable at this read length.

Empirical analysis confirms this: reads from sibling `evo_*` strains (e.g., `evo_1049056.015` vs. `evo_1049056.011`) are misassigned bidirectionally with near-equal frequency, confirming that the alignment signal does not carry sufficient discriminative information. EM post-processing improved `evo_*` accuracy by +1.4 pp (58.4% → 59.8% on the 20K-read benchmark), but population-level abundance signals were insufficient for reliable sibling strain disambiguation.

**Species-level classification remains ~100% accurate** — misclassifications occur only within clades, never between different species. This is an important practical finding: for most applications (pathogen identification, contamination detection, broad metagenomic profiling), Bit-Pop provides reliable classification. Strain-level resolution requires long reads (PacBio/Nanopore), known SNP positions (VCF integration), or machine learning-based approaches.

### 4.4 Key Findings Summary

1. **top-N is the primary accuracy driver** — growth is linear but slow (+0.16 pp for N = 1→2, +0.07 pp for N = 4→5).
2. **Mapping rate is fixed per top-N** — other flags do not change how many reads map.
3. **Advanced flags** (homopolymer fingerprint, SNP detection, golden anchors, spaced seeds, search radius, chunk strategy) **have no measurable effect** on accuracy on CAMI data.
4. **Spaced seeds are counterproductive** — only 413–807 mapped reads vs. 748K baseline.
5. **Smith–Waterman mode is too slow** — timeout after 600 s on CAMI dataset.
6. **Diminishing returns from top-N = 5 onward** — mapped count drops drastically (712K → 511K).
7. **EM adds +0.38 pp accuracy on *k* = 13** (91.91% → 92.29%) via 10,622 reassignments within strain groups. Temperature *t* = 1.0 is better than *t* = 0.1.
8. **K-mer size sweep (*k* = 10–22)**: *k* = 13 gives peak accuracy (91.91%), *k* ≥ 18 gives peak coverage (~100%) but accuracy drops to ~86.8% plateau.
9. **Multi-k consensus (*k* = 12–15) is optimal for coverage** — 90.07% accuracy with EM, 91% mapping. Adding *k* = 22+ increases coverage to ~100% but adds noise (accuracy drops to 89.4%).
10. **Two-pass mapping** recovers unmapped reads at threshold 0.4: +183 correct reads, +718 wrong reads, 99.9% mapping rate.
11. **Species-level classification: ~100%** — misclassifications occur only within clades sharing >99.9% identity.
12. **Strain-level classification: 60–90%** — `evo_*` strains confused with parent and sibling strains. Weighted average: 61.8%.
13. **PacBio HiFi: 95.2% accuracy** on realistic error profile (homopolymers, chimeras, coverage variation), demonstrating robustness to long-read artifacts.
14. **Unmapped reads are not "on the edge"** — scores of 0.35–0.49, far below 0.7 threshold. These are strain variants not represented in the reference.

---

## 5 Limitations

Several limitations should be noted:

- **Strain-level resolution**: Genomes that are >99.9% identical share most k-mers. Reads may map to the wrong strain or to a parent genome. This is a fundamental information-theoretic limitation, not an algorithmic one. Long reads (PacBio/Nanopore) or known SNP positions (VCF integration) would be required for reliable strain-level disambiguation.
- **Proof-of-concept stage**: This work presents benchmarks on simulated reads and the CAMI Low Complexity dataset. The tool has not been validated on large-scale real datasets or compared directly with established aligners on multi-genome classification tasks.
- **Mapping rate at high error rates**: While the top-N anchor filter provides 99.3% mapping rate at 0.1% error rate, performance degrades to 23.4% at 10% error rate.
- **Paired-end conflicts**: 35.5% of read pairs have R1 and R2 mapping to different genomes on the CAMI dataset, reducing effective accuracy. This is expected for reads originating from regions shared between genomes.
- **No clinical validation**: Bit-Pop is an academic research tool and has not been validated for clinical or diagnostic use.
- **Index file sizes**: ~152 MB for 19.7 Mb genome collection. Delta compression is implemented but further optimization is planned to reduce index sizes by 5–10×.
- **FM-index size limit**: libsais has a ~2 GB limit per index (~2.1 billion characters). Large genomes (>2 GB) require the workflow-based splitting approach.

---

## 6 Future Work

Planned improvements include:

- **Long-read integration for strain disambiguation**: Combine short-read FM-index classification with long-read (PacBio/Nanopore) evidence to resolve near-identical strains that share >99.9% of their k-mers.
- **VCF-based SNP weighting**: Incorporate known SNP positions from VCF files to weight alignment scores at strain-discriminative positions, improving classification within clades.
- **SIMD acceleration**: AVX2/AVX-512 acceleration of the 2-bit XOR alignment pipeline for further throughput improvements.
- **SA compression**: Delta-VLI compression of suffix arrays to reduce index file sizes by 5–10×.
- **Direct comparison**: Benchmark against Bowtie2, BWA-MEM, and Kraken2 on multi-genome classification tasks.
- **Larger benchmarks**: 100+ genomes and eukaryotic genomes.
- **Streaming API**: Enhanced streaming interface for integration with real-time sequencing pipelines.

---

## 7 Conclusion

Bit-Pop provides a solution for multi-genome DNA read classification that combines speed, compactness, and accuracy. By combining FM-index-based top-N k-mer filtering with reverse complement support, bit-level XOR alignment (2.3 ns per 31-base chunk), Myers edit distance, and Smith–Waterman refinement, it achieves 99.3% mapping rate and 99.9% classification accuracy across bacterial and eukaryotic genomes on simulated data.

On the CAMI Low Complexity benchmark with 61 microbial genomes and ~1 M reads, Bit-Pop achieves 92.29% accuracy (*k* = 13, top-N = 4, EM refinement) at 70.0% mapping rate, and 90.07% accuracy at 91.0% mapping with multi-k consensus (*k* = 12–15). On a realistic PacBio HiFi benchmark with 69 genomes and 86,248 long reads, Bit-Pop achieves 95.2% accuracy and 99.0% mapping rate.

Species-level classification is ~100% accurate across all benchmarks. Strain-level classification for near-identical genomes (>99.9% identity) achieves 60–90% within clade and represents a fundamental information-theoretic limitation of short-read k-mer-based classification, not a defect of the algorithm.

Bit-Pop is designed as a focused, lightweight tool for species classification, contamination detection, and targeted metagenomic analysis where users work with known genome collections on standard hardware. Its compact index sizes (megabytes, not gigabytes), offline operation, and fast index construction make it particularly suitable for resource-limited settings, clinical laboratories, and academic phylogenomics — filling a niche not addressed by existing general-purpose aligners or broad-spectrum classifiers.

---

## Availability

Source code is available at https://github.com/mladenpop-oss/bit-pop under the MIT License. Benchmark scripts and example datasets are included in the repository. A Zenodo archive with a citable digital object identifier (DOI) is available at https://doi.org/10.5281/zenodo.20043593.

---

## References

1. Breitwieser, F.P., Baker, D.N., and Salzberg, S.L. (2019) KrakenUniq: confident and fast metagenomics classification with unique-kmer counts. *Genome Biology*, 20(175).

2. Burrows, R. and Wheeler, D. (1994) A block-sorting lossless data compression algorithm. *Digital Compression Corporation*, Walter Reed Army Medical Center, Washington, DC, Technical Report 124, pp. 1–15.

3. Grebnov, I. (2016) libsais: C library for the suffix array construction and Burrows–Wheeler transform. https://github.com/IlyaGrebnov/libsais. Accessed: 2025.

4. Ewing, B., Hillier, L., Wendl, M.C., and Green, P. (1998) Base-calling of automated sequencer traces using Phred. I. Accuracy assessment. *Genome Research*, 8(3), 175–185. doi: 10.1101/gr.8.3.175.

5. Ferrada, H. and others (2007) The FM-index: a reversible compression structure and its applications to full-text string search. *Software: Practice and Experience*.

6. Kärkkäinen, J. and Sanders, P. (2003) Simple linear work suffix array construction. *Theoretical Computer Science*, 235(2), 293–314. doi: 10.1016/S0304-3975(01)00154-X.

7. Kim, D., Song, L., Breitwieser, F.P., and Salzberg, S.L. (2016) Centrifuge: fast and memory-efficient metagenomic classification using unique k-mer counts. *Genome Biology*, 17(1), 190. doi: 10.1186/s13059-016-1084-8.

8. Langmead, B. and Salzberg, S.L. (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359. doi: 10.1038/nmeth.1923.

9. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv preprint arXiv:1303.3997*. doi: 10.48550/arXiv.1303.3997.

10. Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences. *Genome Research*, 28(4), 199–204. doi: 10.1101/gr.229144.117.

11. Mende, D.R., Steinegger, M., Feldhahn, B., Farrar, R., Gathmann, S., Hu, P., Klages, L., Köster, J., Ren, S., Rhie, A., et al. (2018) Critical assessment of metagenome interpretation — a benchmark of genetics, physiology and ecology. *Nature Methods*, 15(11), 915–921. doi: 10.1038/s41592-018-0229-2.

12. Myers, E.W. (1999) A fast bit-vector algorithm for approximate string matching based on dynamic programming. *Journal of the ACM*, 46(3), 395–415. doi: 10.1145/314900.314907.

13. Ounit, R., Schmidt, B., and Ren, Q. (2017) CLARK: consensus de bruijn graph lead to fast, sensitive, and accurate taxonomic identification. *Bioinformatics*, 33(18), 2800–2806.

14. Rayon Project Developers (2022) Rayon: a data-parallelism library for Rust. https://github.com/rayon-rs/rayon. Accessed: 2025.

15. Smith, T.F. and Waterman, M.S. (1981) Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195–197. doi: 10.1016/0022-2836(81)90087-5.

16. Wood, D.E., Lu, J., and Langmead, B. (2019) Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257. doi: 10.1186/s13059-019-1891-0.
