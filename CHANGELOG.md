# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Strain marker database support (`--strain-db`, `--strain-only` flags)
- `StrainMarkerDB` struct with JSON/binary save/load

### Changed
- Strain markers reverted: not ready for production (markers-only mode reduces mapping rate)

## [0.2.0] - 2026-05-08

### Added
- `run` command: one-command workflow with auto-build + map
- Paired-end mode with full SAM FLAG handling (TLEN, proper pair detection)
- Reverse complement support with SAM FLAG 0x10
- Top-N rarest k-mer anchors (`--top-n`) for improved mapping rate
- Spaced seed matching (`-s` flag) for high-error reads
- Adaptive k-mer size (`--auto-k`, `--read-type`)
- Myers' bit-vector edit distance alignment (23-54x faster than Smith-Waterman)
- Quality-aware alignment with Phred-scaled penalties
- Memory-mapped FASTA loading (`--mmap` flag)
- Parallel index build via rayon
- CLI progress bars for build and mapping operations
- NCBI E-utilities integration: `search`, `fetch`, `update` commands
- Large genome support via `bitpop-workflow.py` (split/build/map/merge)
- Format v5 persistence with memmap2 and zstd compression
- Delta encoding + VLI compression for position lists

### Changed
- Mapping rate improved from 97.9% to 99.3% with top-N anchors (top_n=3)
- Index loading time reduced to <10ms with memmap persistence
- Refactored project structure into 14 modules

### Fixed
- Rarity calculation bug in multi-genome ranking
- TLEN calculation for paired-end reads
- BWT serialization/deserialization
- Panic on empty genome sequences
- Scaffold-to-genome name mismatch in evaluation scripts

## [0.1.0] - 2026-01-15

### Added
- FM-index construction (SA-IS via libsais, BWT, Occ counter)
- K-mer filter with backward search
- Anchor-based alignment (XOR, Smith-Waterman)
- SAM format output (single-end)
- FASTA/FASTQ parsing (streaming)
- Multi-genome ranking with combined scoring
- Local cache management (`.bitpop` directory)
- Binary serialization (v1, v2 formats)
- CLI with `build`, `map`, `stats`, `load` commands
- Delta encoding + VLI compression
- 253+ unit tests
- 5 integration tests
- 17 Criterion benchmarks
- CAMI Low Complexity benchmark (62 genomes, 50K reads)
- Academic paper with Zenodo DOI
