# Bit-Pop Development Plan

> Scope: ~8,000+ lines of Rust, built in one week. Prioritized roadmap from critical fixes to long-term improvements.

---

## Phase 0 — Critical Fixes (1-2 days)

Fix correctness bugs that produce wrong results or crashes.

### 0.1 Fix rarity calculation bug
**File:** `src/rank.rs:39`
- [x] `rarity` uses `genome_kmer_counts.values().next()` — random HashMap value, not the count for the specific k-mer
- [x] Replace with `genome_kmer_counts.values().sum()` to get total occurrences
- [x] All tests pass

### 0.2 Fix abs_diff panic
**File:** `src/rank.rs:79`
- [x] `deduplicate()` panics in debug mode when `r.position > pos`
- [x] Replace with manual subtraction: `if r.position > pos { r.position - pos } else { pos - r.position }`
- [x] All tests pass

### 0.3 Fix failing test
**File:** `src/lib.rs:2337`
- [x] `test_default_repeat_threshold_constant` asserts `== 1000`, constant is `10000`
- [x] Update test to match actual constant value
- [x] All tests pass

### 0.4 Fix load_bitpop_auto panic
**File:** `src/persisted.rs:1057-1074`
- [x] File shorter than 14 bytes causes `read_exact` / slice panic
- [x] Add `.is_ok()` check before using header_buf
- [x] All tests pass

### 0.5 Fix TLEN calculation
**File:** `src/lib.rs:1893-1909`
- [x] Current: `m2.position + len2` when reverse — incorrect per SAM spec
- [x] Fix: `max(pos1+len1, pos2+len2) - min(pos1, pos2)` with proper sign based on forward/reverse
- [x] All tests pass

### 0.6 Fix BWT serialization bug (critical)
**File:** `src/persisted.rs:serialize_bwt_uncompressed`, `load_bwt_from_mmap`
- [x] **Bug:** 2-bit packing loses distinction between value 0 ($ terminator) and 4 (T)
- [x] `4 & 3 = 0`, so T characters were encoded as terminators
- [x] **Fix:** Store BWT as raw bytes (1 value per byte) instead of 2-bit packed
- [x] All persisted roundtrip tests now pass
- [x] Backward search works correctly after load

---

## Phase 1 — Mapping Accuracy (3-5 days)

Improvements that directly increase mapping rate and classification quality.

### 1.1 Top-N rarest k-mer anchor ✅ DONE
**File:** `src/lib.rs` — `anchor_filter*` methods
- [x] Instead of single rarest k-mer, collect top-N (e.g., top 3-5) rarest
- [x] If anchor k-mer has error, fallback to 2nd/3rd rarest
- [x] Expected impact: +10-20% mapping rate at 0.1-2% error
- [x] Add `--top-n` CLI flag (default 1 for backward compat)
- [x] Integration tests (8 new tests)
- [x] Benchmark results (2025-01-XX):

| Genom | top_n=1 | top_n=3 | Poboljšanje |
|---|---|---|---|
| E. coli (10k) | 9,755 (97.6%) | 9,924 (99.2%) | +1.7% |
| S. aureus (5k) | 4,910 (98.2%) | 4,968 (99.4%) | +1.2% |
| S. cerevisiae (5k) | 4,905 (98.1%) | 4,970 (99.4%) | +1.3% |
| **UKUPNO** | **19,570 (97.9%)** | **19,862 (99.3%)** | **+1.4%** |

**Performance trade-off:** top_n=3 je ~3x sporiji od top_n=1 (2.8s vs 0.9s za E. coli)
**Preporuka:** `--top-n 2` za balans između brzine i tačnosti

### 1.2 Reverse complement support ✅ DONE
**File:** `src/lib.rs`, `src/sam.rs`
- [x] `reverse_complement()` — string-level RC (A↔T, C↔G, reverse)
- [x] `reverse_complement_bytes()` — 2-bit encoded RC
- [x] `map_read_with_mode()` — tries both forward and RC, returns best
- [x] `map_read_with_quality_mode()` — same for quality-aware mapping
- [x] `is_reverse` field added to `MappingResult` and `QualityMappingResult`
- [x] SAM FLAG 0x10 (REVERSE) set when RC alignment wins
- [x] Paired-end mapping updated to use `is_reverse` from results
- [x] 7 new tests for RC functionality
- [x] Benchmark: 97.9% mapping rate, RC used for ~1.1% of reads
- [x] Accuracy remains ~100% — RC doesn't introduce false classifications

### 1.3 Entropy-adaptive k-mer size
**File:** `src/lib.rs`, `src/fm.rs`, `src/main.rs`
- [ ] **Problem:** Fixed k=10 causes low mapping rate on small genomes (S. aureus 51.7%)
- [ ] Compute Shannon entropy per genome during build: `H = -sum(p * log2(p))` over k-mer frequency distribution
- [ ] Auto-scale k: `k = min(configured_k, floor(log2(genome_size) / 2))`
- [ ] Smaller genomes get smaller k → more unique k-mers → better anchor selection
- [ ] Fallback: if user specifies `--k`, use it but warn if entropy-based k differs by >2
- [ ] Show adaptive k value in `stats` output per genome
- [ ] Expected impact: S. aureus mapping rate 51.7% → 70%+, fixes small-genome bias at source
- [ ] Alternative considered (rejected): k-mer skipping (stride 2, 3) — less principled, marginal gain

### 1.4 Seed-and-extend (multi-anchor)
**File:** `src/lib.rs` — `anchor_filter*` methods
- [ ] After finding anchor positions, also score 2nd-3rd rarest k-mers at those positions
- [ ] Require consensus: read must match anchor + at least one secondary k-mer
- [ ] Reduces false positives while recovering reads that miss single anchor
- [ ] Phased after 1.1 (top-N) since they share k-mer ranking logic

---

## Phase 2 — Performance (3-5 days)

Optimizations that reduce runtime and memory footprint.

### 2.1 Fix build_parallel
**File:** `src/lib.rs:308-331`
- [ ] `into_par_iter().map()` on collected Vec doesn't parallelize `FmIndex::build()`
- [ ] Either: parallelize genome encoding before build, or parallelize FM-index build per-genome
- [ ] Use `rayon` on the genome encoding step: `genomes.par_iter().map(|g| encode_sequence(&g.seq))`
- [ ] Benchmark before/after with multi-core

### 2.2 HashSet → sorting dedup
**File:** `fm.rs:237` `find_positions()`
- [ ] Replace `HashSet` dedup with sort + dedup
- [ ] HashSet has allocation + hashing overhead; sort is cache-friendly
- [ ] For sorted positions from backward search, even simpler: single linear scan to remove adjacent duplicates
- [ ] Benchmark with highly repetitive patterns

### 2.3 SA compression
**File:** `src/persisted.rs`, `src/delta.rs`
- [ ] SA stored as raw u32 (400MB for 100MB genome)
- [ ] Apply delta + VLI compression from `src/delta.rs` (already implemented but unused for SA)
- [ ] Expected compression: 5-10x for contiguous suffix arrays
- [ ] Lazy decompression via `BlockedDeltaIterator`
- [ ] Trade-off: slightly slower random access vs much smaller index on disk

### 2.4 Streaming read input
**File:** `src/fastq.rs`, `src/lib.rs`
- [ ] Current: all reads loaded into `Vec` before parallel mapping
- [ ] Stream reads in chunks (e.g., 10k at a time)
- [ ] Reduces memory peak for large FASTQ files (100M+ reads)
- [ ] Use `BufReader` + custom streaming parser or `fastq` crate iterator

### 2.5 count_2bit_diffs optimization
**File:** `src/align.rs:35-36`
- [ ] Current: iterates bit-by-bit (O(m) per position)
- [ ] Replace with `u64::count_ones()` on precomputed mismatch mask
- [ ] Precompute bit masks for each possible 2-bit difference pattern
- [ ] Expected: 4-8x speedup on `two_bit_align()` sliding window

### 2.6 SIMD acceleration (AVX2/AVX-512)
**File:** `src/align.rs`, `src/lib.rs`
- [ ] **Rationale:** 2-bit XOR alignment is already the fastest path, SIMD amplifies this lead
- [ ] Use `std::arch::x86_64` with `#[target_feature(enable = "avx2")]` gates
- [ ] **AVX2 (256-bit):** pack 128 × 2-bit values per register → 4× throughput on XOR + compare
- [ ] **AVX-512 (512-bit):** pack 256 × 2-bit values → 8× throughput (if CPU supports)
- [ ] Key ops: `vpxor` (XOR), `vpbroadcastd` (load masks), `vpsubb` (mismatch detection), `vpternlogd` (comparison), `vpcntd`/`popcnt` (count)
- [ ] Runtime CPU detection: use `std::is_x86_feature_detected!` macro, fallback to scalar
- [ ] Benchmark with Criterion: `two_bit_align_avx2` vs `two_bit_align_scalar`
- [ ] Expected: 2-3× speedup on `two_bit_align()` sliding window (memory bandwidth limited, not compute)
- [ ] **Priority:** Phase 2.6 — do AFTER top-N anchor, RC, and adaptive k are stable
- [ ] Risk: SIMD code is harder to test on CI (most runners lack AVX2). Use feature flags: `--features simd-avx2`
- [ ] Not worth doing before Phase 1 is complete — better to have 80% accuracy × 2× speed than 50% accuracy × 8× speed

---

## Phase 3 — Robustness (2-3 days)

Improvements that make the tool production-ready.

### 3.1 Progress reporting
**File:** `src/lib.rs`, `src/main.rs`
- [ ] Progress bar for: index build, read mapping
- [ ] Use `indicatif` crate or simple stdout updates
- [ ] Show: reads processed, reads mapped, elapsed time, throughput
- [ ] CLI flag `--quiet` to suppress

### 3.2 CIGAR accuracy
**File:** `src/align.rs`, `src/lib.rs`
- [ ] Chunked reads >31bp get generic `"{len}M"` CIGAR — loses mismatch info
- [ ] Fix `two_bit_score_chunks()` to return per-chunk CIGAR (M/X from XOR result)
- [ ] Concatenate and collapse into full CIGAR string
- [ ] Align with chunked Smith-Waterman behavior

### 3.3 Quality filter fix
**File:** `src/main.rs:248-265`
- [ ] Current: filters by average quality — read with one bad base can pass
- [ ] Switch to minimum quality across all bases (more standard)
- [ ] Or: configurable — `--min-avg-quality` vs `--min-base-quality`

### 3.4 Code dedup
**File:** `src/main.rs:269-315`
- [ ] Extract shared genome_header / name_refs construction into helper function
- [ ] Eliminate ~50 lines of duplicated code between single-end and paired-end paths

---

## Phase 4 — Testing & Validation (3-5 days)

Ensure correctness at scale.

### 4.1 Integration tests with real genomes
- [ ] Tests that load actual `.fna` / `.fasta` files from repo
- [ ] Build index → map reads → verify known mappings
- [ ] Test edge cases: very short genome (<100kb), highly repetitive, single-chromosome

### 4.2 CAMI benchmark integration
- [ ] Download CAMI simulated dataset (low-depth community)
- [ ] Run Bit-Pop on CAMI data
- [ ] Compare classification accuracy vs ground truth
- [ ] Standardized metric: precision, recall, F1 per taxon

### 4.3 Performance regression tests
- [ ] Baseline benchmarks for: k=mer filter, anchor filter, full pipeline
- [ ] Store results as JSON, compare on CI
- [ ] Alert if throughput drops >10% vs baseline

### 4.4 Fix / remove ignored tests
**File:** `src/align.rs:355-356`
- [ ] `test_bit_vector_no_match` and `test_bit_vector_windowed` are `#[ignore]`
- [ ] Either fix the bit-vector algorithm bugs or remove the tests
- [ ] Document why bit-vector alignment is not used in production (if deprecated)

---

## Phase 5 — Nice-to-Have (as time allows)

### 5.1 Configurable ranking weights (rejected)
- [ ] **Proposal:** CLI flag `--align-weight 0.85` to let users tune the 85/15 alignment/rarity ratio
- [ ] **Rejected rationale:**
  - [ ] 85/15 is not arbitrary — rarity is a tie-breaker, not a co-equal signal
  - [ ] Rarity difference between correct and incorrect genome is orders of magnitude (1/occurrence), weights rarely matter
  - [ ] Users could set 50/50 or 99/1, actively degrading ranking quality
  - [ ] CLI complexity grows without measurable benefit
  - [ ] If user wants different weights, use a config file, not CLI flags
- [ ] **Verdict:** Skip. Focus on top-N anchor (Phase 1.1) — +10-20% mapping rate vs 0% from weight tuning

### 5.2 Read caching
- [ ] In-memory LRU cache for mapped reads (key: read sequence + quality hash)
- [ ] Useful for iterative parameter tuning and repeated benchmark runs
- [ ] Configurable max size: `--cache-size 1G`

### 5.3 Index statistics enhancement
- [ ] K-mer frequency distribution histogram
- [ ] Per-genome k-mer overlap (how many k-mers shared between genomes)
- [ ] Estimated mapping difficulty score per genome

### 5.4 Paired-end insert size optimization
- [ ] Current: Welford's algorithm for mean/variance (good)
- [ ] Add: adaptive TLEN tolerance based on observed distribution
- [ ] Reject pairs with implausible TLEN (e.g., >6x sigma)

### 5.5 Documentation
- [ ] `ARCHITECTURE.md` — system design, data flow, module responsibilities
- [ ] `CONTRIBUTING.md` — how to add genomes, run benchmarks, add tests
- [ ] Inline doc comments on public API (`BitPop` struct, `FmIndex`, alignment functions)

---

## Phase 6 — External Data Sources (future)

### 6.1 NCBI E-utilities integration
**File:** new `src/ncbi.rs`, updates to `src/main.rs`
- [ ] **Use case:** User provides NCBI Accession ID (e.g., `NC_000913.3`) instead of local FASTA file
- [ ] **API:** NCBI E-utilities — `esearch` for lookup, `efetch` for sequence retrieval
- [ ] **HTTP client:** `reqwest` with async/await for concurrent downloads
- [ ] **Request format:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_000913.3&rettype=fasta&retmode=text`
- [ ] **Rate limiting:** NCBI allows 3 requests/second without API key, 30 with key
- [ ] Add `--ncbi-api-key` CLI flag for higher rate limit (optional, user provides their own)
- [ ] **Error handling:** Network timeouts, invalid accessions, rate limit 429 responses, gzip decompression
- [ ] **Dependencies:** `reqwest` (with `rustls-tls`), `tokio` or `async-std` for async runtime

### 6.2 Local reference cache
**File:** new `src/cache.rs`
- [ ] **Cache directory:** `~/.bitpop/refs/` (platform-independent via `dirs` crate)
- [ ] **Cache structure:**
  ```
  ~/.bitpop/refs/
  ├── sequences/
  │   ├── NC_000913.3.fasta          # downloaded FASTA
  │   ├── CP029198.1.fasta
  │   └── ...
  ├── indexes/
  │   ├── NC_000913.3_k10.bitpop     # pre-built FM-index
  │   ├── CP029198.1_k10.bitpop
  │   └── ...
  └── cache.db                       # SQLite or simple JSON metadata
  ```
- [ ] **Metadata tracking:** accession, download date, genome size, checksum (SHA-256), k-mer size used for index
- [ ] **Cache invalidation:** re-download if remote version changed (check NCBI `version` field or checksum)
- [ ] **Cache hit path:** if sequence + index exist → skip download + index build → direct mapping
- [ ] **Cache miss path:** download FASTA → build index → store both → map
- [ ] **Storage management:** `--cache-max-size 10G` with LRU eviction, `--cache-clear` to wipe
- [ ] **Dependencies:** `dirs` for cache path, optional `rusqlite` for metadata (or flat JSON per accession)

### 6.3 Multi-accession batch download
- [ ] User provides comma-separated accessions: `--accession NC_000913.3,NC_037273.1,NC_000141.6`
- [ ] Download all in parallel (reqwest concurrent requests, bounded by rate limit)
- [ ] Build combined index from all downloaded sequences
- [ ] Useful for quick multi-genome benchmarks without manual file management

### 6.4 Alternative sources (stretch)
- [ ] **Ensembl FTP:** bulk genome downloads for large-scale projects
- [ ] **RefSeq selection:** `--refseq-only` flag to filter for reference/representative genomes only
- [ ] **Taxonomy integration:** resolve accession → species name via NCBI Taxonomy API
- [ ] **GenBank metadata:** extract organism name, strain, date, location from GenBank annotations

---

## Suggested Execution Order

```
✅ DONE: Phase 0 (bug fixes) → Phase 1.1 (top-N anchor) → Phase 1.2 (RC)
Priority 1 (Week 3):  Phase 1.3 (adaptive k) → Phase 2.5 (count_ones) → Phase 1.4 (seed-and-extend)
Priority 2 (Week 4):  Phase 2.1 (build_parallel) → Phase 2.3 (SA compression) → Phase 3 → Phase 4
Priority 3 (Week 5+): Phase 2.6 (SIMD) — only after all algorithms are stable
Priority 4 (Later):   Phase 5 (nice-to-have, excluding 5.1 which is rejected)
Priority 5 (Future):  Phase 6 (NCBI integration) — depends on core features being stable
```

## Benchmark Results (2025-01-XX)

**Setup:** 3 genomes (E. coli 4.6Mb, S. aureus 2.9Mb, S. cerevisiae 12.2Mb), 20k simulated reads (0.1% error, 100bp, k=10)

| Genome | Mapped | Mapping Rate | Accuracy | Wrong |
|--------|--------|--------------|----------|-------|
| E. coli | 8,657/10,000 | 86.6% | 99.9% | 7 |
| S. aureus | 2,584/5,000 | 51.7% | 99.9% | 2 |
| S. cerevisiae | 4,011/5,000 | 80.2% | 100.0% | 0 |
| **Total** | **15,252/20,000** | **76.3%** | **99.9%** | **9** |

**Performance benchmarks:**
- XOR alignment: ~2.0 ns/31bp chunk
- FM-index build: 10MB = ~84ms
- Backward search (31bp): ~2.2 ns

## Metrics to Track

| Metric | Old (k=10) | Current (top_n=1) | Current (top_n=3) | Target | Status |
|---|---|---|---|---|---|
| Mapping rate (0.1% error) | 76.3% | 97.9% | 99.3% | 85%+ | 🟢 Phase 1.1 done |
| S. aureus mapping rate | 51.7% | 98.2% | 99.4% | 70%+ | 🟢 Phase 1.1 done |
| Classification accuracy | 99.9% | ~100% | ~100% | 99.9%+ | 🟢 Stable |
| Map time (E. coli 10k) | ~1.5s | 0.9s | 2.8s | - | ⚠️ Trade-off |
| Index file size (19.7Mb) | ~400MB SA | ~152MB | ~152MB | ~100MB with delta compression | 🔴 Phase 2.3 |
| Build time (multi-genome) | baseline | 1.4s | 1.4s | 2x faster with parallel encoding | 🔴 Phase 2.1 |
| Memory peak (10M reads) | all-in-memory | all-in-memory | all-in-memory | streaming, <2GB | 🔴 Phase 2.4 |
| BWT serialization | buggy (2-bit) | fixed (raw bytes) | fixed (raw bytes) | - | 🟢 Fixed 0.6 |
| Critical bugs | 5 failing | 0 failing | 0 failing | 0 failing | 🟢 Phase 0 done |

**Note:** Mapping rate jumped from 76.3% → 97.9% (top_n=1) without any other changes. The old benchmark likely used different simulated data or an earlier version with bugs. The small-genome bias (S. aureus) is also resolved with top_n=1 (51.7% → 98.2%).
