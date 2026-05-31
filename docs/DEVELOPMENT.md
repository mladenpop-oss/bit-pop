# Bit-Pop Development

## Project Structure

```
bit-pop/
├── src/
│   ├── lib.rs          # Core library (MappingResult, pipeline)
│   ├── main.rs         # CLI entry point
│   ├── fm.rs           # FM-index (SA-IS, BWT, backward search)
│   ├── align.rs        # Alignment (XOR, SW, Myers edit distance)
│   ├── rank.rs         # Multi-genome ranking (alignment + k-mer rarity)
│   ├── em.rs           # EM post-processing (soft assignment)
│   ├── sam.rs          # SAM/BAM output (CIGAR, NM, MAPQ tags)
│   ├── fasta.rs        # FASTA parsing
│   ├── fastq.rs        # FASTQ parsing
│   ├── taxonomy.rs     # NCBI taxonomy tree + LCA algorithm
│   ├── ncbi.rs         # NCBI E-utilities integration
│   └── ...
├── bin/
│   └── extract_seqs.rs # Extract genome sequences from .bitpop index
├── benches/            # Criterion benchmarks (17 groups)
├── tests/              # Integration tests (5 tests)
├── scripts/
│   ├── simulate_reads.py        # Realistic read simulation (Biopython)
│   ├── cami_accuracy.py         # SAM accuracy vs CAMI ground truth
│   ├── consensus_base.py        # Multi-k consensus (legacy, use bit-pop consensus)
│   ├── bitpop-workflow.py       # Large genome workflow tool
│   └── em_classifier.py         # Python EM prototype (reference)
├── data/
│   ├── genomes/        # Reference genomes (.fna, .fasta)
│   └── reads/          # Sequencing reads (.fastq)
├── docs/
│   ├── paper.tex       # Academic paper (in preparation)
│   ├── paper.pdf       # Compiled paper
│   ├── references.bib  # Bibliography
│   └── CITATION.cff    # Citation metadata
├── gui/                # Desktop GUI (Svelte + Tauri)
│   ├── src/            # Svelte frontend
│   └── src-tauri/      # Tauri Rust backend
└── Cargo.toml
```

---

## Roadmap

### ✅ Completed

**Core engine:**
- FM-index with SA-IS construction, BWT, backward search
- 2-bit DNA encoding + XOR alignment (~2.3 ns per 31-base chunk)
- Myers edit distance (23-54x faster than Smith-Waterman)
- Multi-genome ranking: alignment score × k-mer rarity
- Reverse complement support (SAM FLAG 0x10)
- Parallel mapping and index build (rayon)

**I/O:**
- SAM output with full CIGAR, NM, MAPQ, MD tags
- BAM output with BGZF compression
- Paired-end support with discordant pair reconciliation
- Gaussian insert size model for paired-end classification
- Memory-mapped FASTA (`--mmap`)

**Features:**
- Top-N rarest k-mer anchors (`--top-n`)
- EM post-processing (soft assignment, temperature=1.0)
- Multi-k consensus (`bit-pop consensus`)
- Auto chunking: ON for reads <1000bp, OFF for long reads
- Two-pass mapping (`--two-pass`)
- NCBI E-utilities integration (`--ncbi`)
- Taxonomic classification with LCA (`bit-pop tax`)
- Auto index caching (reuses `.bitpop` when genomes unchanged)
- Large genome workflow (`scripts/bitpop-workflow.py`)
- Desktop GUI (Tauri + Svelte)
- Docker container
- Bioconda package

**Benchmarks:**
- Quick benchmark: 3 genomes, 99.9% accuracy
- CAMI Low: 61 genomes, ~1M reads, 92.29% strain / ~100% species accuracy
- PacBio HiFi: 69 genomes, 86k long reads, 95.2% accuracy

### 🔧 In Progress

- Paper submission (Oxford Bioinformatics, target: August 2026)
- RNA-Pop: FM-index engine for RNA-seq quantification (early testing)

### 📋 Planned

- SIMD acceleration (AVX2) for XOR alignment
- SA compression for reduced index size
- Streaming input for large FASTQ files
- CIGAR accuracy improvements for chunked reads
- Unified multi-index FM-index (>2GB genomes natively)
- VCF integration for known SNP positions (strain disambiguation)
- ONT long-read benchmark (real data)
- 100+ genome benchmark
- docs.rs API documentation

---

## Testing

```bash
# Run all tests
cargo test

# Integration tests only
cargo test --test integration_tests

# Benchmarks
cargo bench
```

**Coverage:**
- 312+ unit tests (alignment, indexing, serialization, SAM, EM, taxonomy)
- 5 integration tests (build, map, multi-genome, SAM format, cache reuse)
- 17 Criterion benchmark groups

---

## Architecture Notes

### Scoring Formula

```
final_score = alignment_score × (0.5 + 0.5 × relevance)
relevance   = 0.4 × rarity + 0.6 × proximity
rarity      = 1 / genome_kmer_count  (per-genome, not global sum)
proximity   = kmer_count / read_length
```

**Note:** Rarity is computed per-genome (not summed across all genomes). This was a critical bug fix that improved CAMI accuracy from ~89% to ~92%.

### EM Algorithm

- Temperature=1.0 recommended (0.1 over-concentrates, suppresses reassignment)
- Convergence threshold: KL divergence < 0.001
- Max iterations: 50 (typically converges in 10)
- Effect: +0.38pp strain accuracy on CAMI k13

### Long Read Mode

Chunking is automatically disabled for reads >1,000 bp. Direct FM-index alignment for long reads outperforms chunk-based voting by +12.6pp accuracy (PacBio HiFi benchmark).

### Index Format (.bitpop)

Binary serialized FM-index:
- BWT + suffix array (libsais, max ~2GB per index)
- k-mer occurrence table per genome
- Genome name → ID mapping
- Build time: ~17 seconds for 157 Mb (61 genomes)

---

## Known Limitations

- libsais 2GB limit per index — use `scripts/bitpop-workflow.py` for larger genomes
- Strain-level resolution limited to ~60-92% for >99.9% identical genomes (information-theoretic limit)
- Chunked reads use generic CIGAR without per-base mismatch detail
- Not validated for clinical use

---

## Related Projects

- [bit-pep](https://github.com/animesh/bit-pep) — FM-index engine adapted for proteomics (fork by Animesh Sharma, NTNU)
- RNA-Pop — RNA-seq quantification on the same engine (in development)
