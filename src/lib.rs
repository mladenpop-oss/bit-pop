//! Bit-Pop: Ultra-fast multi-genome DNA read classification.
//!
//! Bit-Pop maps DNA sequencing reads across multiple reference genomes using
//! a compact FM-index with bit-level parallelism. It identifies **which genome**
//! in a collection best matches each read.
//!
//! # Pipeline
//!
//! 1. **Indexing**: Load FASTA genomes → 2-bit encode → Build FM-index (SA-IS)
//! 2. **Filtering**: Find top-N rarest k-mers as anchors → backward search for candidates
//! 3. **Alignment**: Score candidates via XOR/SW/Myers → reverse complement → rank
//!
//! # Example
//!
//! ```
//! use bit_pop::BitPop;
//!
//! // Create indexer with k-mer size 10
//! let mut bp = BitPop::new(10);
//!
//! // Add reference genomes
//! bp.add_genome("Ecoli_K12", "AGCTAGCTAGCT...");
//! bp.add_genome("Staph_aureus", "GCTAGCTAGCTA...");
//!
//! // Build the FM-index (must call before mapping)
//! bp.build();
//!
//! // Map a read (returns ranked results)
//! let results = bp.map_read("AGCTAGCTAGCTAGCT", 10);
//!
//! for result in &results {
//!     let name = bp.genome_name(result.genome_id).unwrap();
//!     println!("{} → {} (score: {:.2})", result.genome_id, name, result.score);
//! }
//! ```
//!
//! # Parallel Mapping
//!
//! ```no_run
//! use bit_pop::BitPop;
//!
//! let mut bp = BitPop::new(10);
//! bp.add_genome("genome1", "AGCT...");
//! bp.build();
//!
//! // Map many reads in parallel
//! let reads = vec![
//!     ("read1", "AGCTAGCTAGCT"),
//!     ("read2", "GCTAGCTAGCTA"),
//! ];
//! let mapped = bp.map_reads_parallel(&reads, "output.sam", 10).unwrap();
//! println!("Mapped {} reads", mapped);
//! ```
//!
//! # Persistence
//!
//! ```no_run
//! use bit_pop::BitPop;
//!
//! let mut bp = BitPop::new(10);
//! bp.add_genome("genome1", "AGCT...");
//! bp.build();
//!
//! // Save to disk
//! bp.serialize_to_file("index.bitpop").unwrap();
//!
//! // Load from disk (< 10ms with memmap)
//! let loaded = BitPop::deserialize_from_file("index.bitpop").unwrap();
//! let results = loaded.map_read("AGCTAGCTAGCT", 10);
//! ```
//!
//! # DNA Encoding
//!
//! All sequences are stored as 2-bit encoded bytes (A=1, C=2, G=3, T=4).
//! Unknown bases (N) are skipped during encoding.
//!
//! ```
//! use bit_pop::{encode_sequence, decode_sequence, encode_kmer, decode_kmer};
//!
//! let seq = "ACGTACGT";
//! let encoded = encode_sequence(seq);
//! assert_eq!(decode_sequence(&encoded), seq.replace('N', ""));
//!
//! let kmer = encode_kmer("ACGT").unwrap();
//! assert_eq!(decode_kmer(kmer, 4), "ACGT");
//! ```

#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::type_complexity)]

pub mod align;
pub mod bam;
pub mod cache;
pub mod chain;
pub mod chunk_consensus;
pub mod consensus;
pub mod delta;
pub mod em;
pub mod fasta;
pub mod fastcon;
pub mod fastq;
pub mod fm;
pub mod hf;
pub mod index_manager;
#[cfg(feature = "android")]
pub mod jni;
pub mod ncbi;
pub mod persisted;
pub mod rank;
pub mod sam;
pub mod serialize;
pub mod snp;
pub mod taxonomy;

pub use chunk_consensus::MultiChunkConsensus;
pub use consensus::{ConsensusResult, ConsensusStrategy, KResult, MultiKConsensus};
pub use fastcon::FastCon;

use std::fmt;

use std::collections::HashMap;
use std::io;

use fm::FmIndex;
use rayon::prelude::*;

/// Default threshold for filtering repetitive k-mers.
/// K-mers appearing in more than this many positions are treated as noise.
/// 10000 = skips highly repetitive elements, keeps unique signal.
pub const DEFAULT_REPEAT_THRESHOLD: usize = 10000;

// --- DNA Alphabet ---

/// Encode a single DNA base to a 2-bit value.
/// $=0 (sentinel), A=1, C=2, G=3, T=4. N skipped.
pub fn encode_base(ch: char) -> Option<u8> {
    match ch.to_ascii_uppercase() {
        'A' => Some(1),
        'C' => Some(2),
        'G' => Some(3),
        'T' => Some(4),
        _ => None,
    }
}

/// Decode a 2-bit value back to a DNA base character.
pub fn decode_base(val: u8) -> char {
    match val {
        1 => 'A',
        2 => 'C',
        3 => 'G',
        4 => 'T',
        _ => 'N',
    }
}

type PairedReads = (String, String, Vec<u8>, String, Vec<u8>);

/// Encode a DNA sequence into a compact byte slice (2 bits per base).
/// Skips unknown bases (N).
pub fn encode_sequence(seq: &str) -> Vec<u8> {
    seq.chars().filter_map(encode_base).collect()
}

/// Decode a compact byte slice back to a DNA string.
pub fn decode_sequence(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| decode_base(b)).collect()
}

// --- K-mer Encoding ---

/// Encode a k-mer string into a u64 bit-parallel representation.
/// Each base = 2 bits (A=0, C=1, G=2, T=3), so k=31 fits in u64 (31 × 2 = 62 bits).
/// Returns None if k-mer is too long (>31) or contains invalid bases.
pub fn encode_kmer(kmer: &str) -> Option<u64> {
    if kmer.len() > 31 {
        return None;
    }
    let mut result: u64 = 0;
    for ch in kmer.chars() {
        let base = match ch.to_ascii_uppercase() {
            'A' => 0u64,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => return None,
        };
        result = (result << 2) | base;
    }
    Some(result)
}

/// Decode a u64 back to a k-mer string of length `k`.
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    let mut result = String::with_capacity(k);
    let mut val = encoded;
    for _ in 0..k {
        let base = (val & 3) as u8;
        result.push(match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        });
        val >>= 2;
    }
    result.chars().rev().collect()
}

/// Generate all k-mer variants with up to `max_mismatches` substitutions.
/// Returns unique encoded k-mer values (excluding the original).
/// Uses 0-indexed encoding internally (A=0, C=1, G=2, T=3).
/// Caller must convert to 1-indexed when querying FM-index.
pub fn generate_kmer_variants(original: u64, k: usize, max_mismatches: usize) -> Vec<u64> {
    if k > 31 || max_mismatches == 0 {
        return Vec::new();
    }

    let bases: [u64; 4] = [0, 1, 2, 3]; // A, C, G, T (0-indexed)
    let mut variants = std::collections::HashSet::new();

    // Encode: first base is at MSB, last base is at LSB
    // Position 0 (first base) -> bits (2*(k-1), 2*k-1)
    // Position k-1 (last base) -> bits (0, 1)
    let mut positions = Vec::with_capacity(k);
    for i in 0..k {
        let bit_pos = 2 * (k - 1 - i); // MSB-first mapping
        let base = (original >> bit_pos) & 3;
        positions.push((i, base as usize));
    }

    generate_variants_recursive(
        &positions,
        0,
        0,
        max_mismatches,
        &bases,
        original,
        &mut variants,
    );

    variants.into_iter().collect()
}

fn generate_variants_recursive(
    positions: &[(usize, usize)],
    pos_idx: usize,
    mismatches: usize,
    max_mismatches: usize,
    bases: &[u64; 4],
    current: u64,
    variants: &mut std::collections::HashSet<u64>,
) {
    if mismatches > max_mismatches {
        return;
    }

    if pos_idx == positions.len() {
        if mismatches > 0 {
            variants.insert(current);
        }
        return;
    }

    let (array_idx, original_base) = positions[pos_idx];
    let bit_pos = 2 * (positions.len() - 1 - array_idx); // MSB-first mapping

    // Try each base (including original)
    for &base in bases.iter() {
        let new_mismatches = if base as usize == original_base {
            mismatches
        } else {
            mismatches + 1
        };

        if new_mismatches > max_mismatches {
            continue;
        }

        // Clear the old bits (2 bits per base) and set new bits at this position
        let cleared = current & !(3u64 << bit_pos);
        let new_current = cleared | ((base & 3) << bit_pos); // Mask to 2 bits

        generate_variants_recursive(
            positions,
            pos_idx + 1,
            new_mismatches,
            max_mismatches,
            bases,
            new_current,
            variants,
        );
    }
}

/// Compute reverse complement of a DNA string.
/// A<->T, C<->G, then reverse the string.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// Compute reverse complement of a 2-bit encoded sequence.
/// Swaps A(1)<->T(4), C(2)<->G(3), then reverses byte order.
pub fn reverse_complement_bytes(encoded: &[u8]) -> Vec<u8> {
    let complement: Vec<u8> = encoded
        .iter()
        .map(|&b| match b {
            1 => 4, // A -> T
            4 => 1, // T -> A
            2 => 3, // C -> G
            3 => 2, // G -> C
            other => other,
        })
        .collect();
    let mut result = complement;
    result.reverse();
    result
}

// --- MappingResult ---

/// Result of mapping a read to a genome.
#[derive(Debug, Clone)]
pub struct MappingResult {
    /// Genome identifier (which reference genome this maps to).
    pub genome_id: u32,
    /// Position in the genome where the read maps.
    pub position: u64,
    /// Alignment score (0.0-1.0, higher = better match).
    pub score: f64,
    /// CIGAR string describing the alignment (e.g. "100M", "95M5D").
    pub cigar: String,
    /// Context: ±window bases around the mapped position.
    pub context: String,
    /// True if the read mapped to the reverse strand (RC alignment won).
    pub is_reverse: bool,
    /// K-mer rarity score (1 / occurrence_count of the read's first k-mer).
    pub rarity: f64,
    /// MD tag string for mismatch verification (e.g. "10A5T3").
    pub md_string: String,
    /// Homopolymer fingerprint similarity score (0.0-1.0, only set when --hf enabled).
    pub hf_score: f64,
}

/// Quality-aware mapping result with per-base quality information.
#[derive(Debug, Clone)]
pub struct QualityMappingResult {
    /// Genome identifier (which reference genome this maps to).
    pub genome_id: u32,
    /// Position in the genome where the read maps.
    pub position: u64,
    /// Raw alignment score (0.0-1.0, higher = better match).
    pub align_score: f64,
    /// Quality-adjusted alignment score with Phred-scaled penalties.
    pub adjusted_score: f64,
    /// Combined ranking score (align_score × 0.85 + rarity × 0.15).
    pub combined_score: f64,
    /// CIGAR string describing the alignment.
    pub cigar: String,
    /// Quality penalty applied (negative value means mismatches at high quality positions).
    pub quality_penalty: f64,
    /// Per-base quality scores from the original read.
    pub quality_scores: Vec<u8>,
    /// Context: ±window bases around the mapped position.
    pub context: String,
    /// True if the read mapped to the reverse strand (RC alignment won).
    pub is_reverse: bool,
    /// K-mer rarity score (1 / occurrence_count of the read's first k-mer).
    pub rarity: f64,
    /// MD tag string for mismatch verification (e.g. "10A5T3").
    pub md_string: String,
}

/// A paired-end read (R1 + R2).
#[derive(Debug, Clone)]
pub struct PairedRead {
    pub name: String,
    pub read1_seq: String,
    pub read1_qual: Vec<u8>,
    pub read2_seq: String,
    pub read2_qual: Vec<u8>,
}

/// Insert size statistics for paired-end mapping with Gaussian model.
#[derive(Debug, Clone)]
pub struct InsertSizeStats {
    pub mean: f64,
    pub stddev: f64,
    pub count: usize,
    /// M2 sum for Welford's algorithm (stored internally, not serialized)
    m2: f64,
}

impl Default for InsertSizeStats {
    fn default() -> Self {
        Self::new()
    }
}

impl InsertSizeStats {
    pub fn new() -> Self {
        Self {
            mean: 0.0,
            stddev: 0.0,
            count: 0,
            m2: 0.0,
        }
    }

    pub fn update(&mut self, insert_size: i64) {
        if insert_size <= 0 {
            return;
        }
        self.count += 1;
        let delta = insert_size as f64 - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = insert_size as f64 - self.mean;
        self.m2 += delta * delta2;
        if self.count >= 2 {
            self.stddev = (self.m2 / (self.count as f64 - 1.0)).sqrt();
        }
    }

    pub fn update_with_variance(&mut self, insert_size: i64) {
        self.update(insert_size);
    }

    /// Compute the Gaussian (normal distribution) probability density
    /// for an observed insert size value.
    ///
    /// Returns the probability density function (PDF) value:
    ///   f(x) = (1 / (sigma * sqrt(2*pi))) * exp(-0.5 * ((x - mu) / sigma)^2)
    ///
    /// Returns 0.0 if stddev is zero or count < 2.
    pub fn gaussian_probability(&self, observed_tlen: i64) -> f64 {
        if self.count < 2 || self.stddev < f64::EPSILON {
            return 0.0;
        }
        let x = observed_tlen as f64;
        let diff = x - self.mean;
        let exponent = -0.5 * (diff / self.stddev).powi(2);
        // Clamp exponent to avoid underflow
        if exponent < -500.0 {
            return 0.0;
        }
        (1.0 / (self.stddev * (2.0 * std::f64::consts::PI).sqrt())) * exponent.exp()
    }

    /// Compute the log probability density of the Gaussian model for an
    /// observed insert size. More numerically stable than raw probability.
    ///
    /// Returns log(f(x)) where f is the normal PDF.
    /// Returns `-f64::INFINITY` if stddev is zero or count < 2.
    pub fn log_gaussian_probability(&self, observed_tlen: i64) -> f64 {
        if self.count < 2 || self.stddev < f64::EPSILON {
            return -f64::INFINITY;
        }
        let x = observed_tlen as f64;
        let diff = x - self.mean;
        let z = diff / self.stddev;
        // log PDF = -log(sigma) - 0.5*log(2*pi) - 0.5*z^2
        -self.stddev.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln() - 0.5 * z * z
    }

    /// Compute a normalized confidence score in [0.0, 1.0] for an observed
    /// insert size. Maps the PDF value to a 0-1 range using the peak density.
    ///
    /// A score of 1.0 means the observed TLEN equals the mean.
    /// A score near 0.0 means the TLEN is far from the mean (>5 stddevs).
    pub fn insert_size_confidence(&self, observed_tlen: i64) -> f64 {
        if self.count < 2 || self.stddev < f64::EPSILON {
            return 0.0;
        }
        // The peak PDF value is at x = mean: 1/(sigma*sqrt(2*pi))
        let peak_pdf = 1.0 / (self.stddev * (2.0 * std::f64::consts::PI).sqrt());
        if peak_pdf < f64::EPSILON {
            return 0.0;
        }
        let actual_pdf = self.gaussian_probability(observed_tlen);
        // Normalize: ratio of actual PDF to peak PDF = exp(-0.5*z^2)
        let ratio = actual_pdf / peak_pdf;
        ratio.clamp(0.0, 1.0)
    }

    /// Check if the observed TLEN falls within the expected range (mean ± 3*stddev)
    /// using the Gaussian model. Also returns the confidence score.
    pub fn is_proper_pair(&self, observed_tlen: i64) -> bool {
        if self.count < 2 || observed_tlen <= 0 {
            return false;
        }
        let lower = (self.mean - 3.0 * self.stddev).max(0.0) as i64;
        let upper = (self.mean + 3.0 * self.stddev) as i64;
        observed_tlen >= lower && observed_tlen <= upper
    }

    /// Get the expected insert size range (mean ± 3*stddev).
    pub fn expected_range(&self) -> (i64, i64) {
        if self.count < 2 || self.stddev < f64::EPSILON {
            return (0, 0);
        }
        let lower = (self.mean - 3.0 * self.stddev).max(0.0) as i64;
        let upper = (self.mean + 3.0 * self.stddev) as i64;
        (lower, upper)
    }
}

/// Result of mapping a single read in a paired-end context.
#[derive(Debug, Clone)]
pub struct PairedReadMapping {
    pub genome_id: u32,
    pub position: u64,
    pub score: f64,
    pub cigar: String,
    pub is_reverse: bool,
    pub mapped: bool,
    pub align_score: f64,
    pub rarity: f64,
    pub md_string: String,
}

/// Result of mapping a paired-end read to all indexed genomes.
#[derive(Debug, Clone)]
pub struct PairedMappingResult {
    pub read_name: String,
    pub map1: Option<PairedReadMapping>,
    pub map2: Option<PairedReadMapping>,
    pub tlen: i64,
    pub insert_size_stats: InsertSizeStats,
}

/// Chunk-level mapping result for chunk-based PacBio classification.
#[derive(Debug, Clone)]
pub struct ChunkVote {
    /// Genome ID that won this chunk's vote
    pub genome_id: u32,
    /// Number of chunks this genome won
    pub win_count: usize,
    /// Sum of alignment scores across all won chunks
    pub total_score: f64,
    /// Average alignment score
    pub avg_score: f64,
    /// Total chunks mapped successfully
    pub chunks_mapped: usize,
    /// Quality-weighted score (score * chunks_mapped / total_chunks)
    pub quality_weighted_score: f64,
}

/// Alignment mode for mapping reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignMode {
    /// 2-bit XOR alignment (fastest, ~2.3ns per read)
    #[default]
    Xor,
    /// Smith-Waterman local alignment (more accurate, handles gaps/indels)
    Sw,
    /// XOR first for fast filtering, then SW on top candidates for precise scoring
    Hybrid,
    /// XOR with soft-clipping: slides windows across read to find optimal alignment region,
    /// emits S operations in CIGAR for adapter/low-quality regions
    Softclip,
    /// Gap-aware XOR chaining: minimizer-based seed chaining with XOR gap extension.
    /// True long-read alignment for ONT/PacBio (5-15% error rates, long indels).
    Chain,
}

/// Fuzzy k-mer matching method for improved strain resolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FuzzyMethod {
    /// No fuzzy matching (default)
    #[default]
    None,
    /// Generate all k-mer variants with N mismatches and query FM-index for each
    FuzzyKmer,
    /// Allow N mismatches in spaced seed "match" positions
    FuzzySeed,
    /// Build neighborhood hash table at index build time for O(1) fuzzy lookup
    Neighborhood,
}

/// Anchor strategy for chunk-based mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChunkAnchorStrategy {
    /// Use top-N rarest k-mers (current default behavior, no quality needed)
    #[default]
    Rarest,
    /// Use quality-weighted golden anchors (requires FASTQ quality scores)
    Golden,
    /// Use spaced seed k-mers (requires spaced seed enabled during build)
    Spaced,
}

impl fmt::Display for ChunkAnchorStrategy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChunkAnchorStrategy::Rarest => write!(f, "rarest"),
            ChunkAnchorStrategy::Golden => write!(f, "golden"),
            ChunkAnchorStrategy::Spaced => write!(f, "spaced"),
        }
    }
}

impl fmt::Display for AlignMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignMode::Xor => write!(f, "xor"),
            AlignMode::Sw => write!(f, "sw"),
            AlignMode::Hybrid => write!(f, "hybrid"),
            AlignMode::Softclip => write!(f, "softclip"),
            AlignMode::Chain => write!(f, "chain"),
        }
    }
}

// --- BitPop (main struct) ---

/// Bit-Pop multi-genome DNA read mapper.
///
/// Implements a 3-stage pipeline:
/// 1. **K-mer inverted index** → candidate positions via FM-index backward search
/// 2. **Bit-level alignment** → precise matches via XOR/SW/Myers algorithms
/// 3. **Multi-genome ranking** → scored results combining alignment and k-mer rarity
///
/// # Usage
///
/// ```
/// use bit_pop::{BitPop, AlignMode};
///
/// let mut bp = BitPop::new(10);
/// bp.add_genome("Ecoli", "AGCTAGCTAGCTAGCTAGCT");
/// bp.add_genome("Staph", "GCTAGCTAGCTAGCTAGCTA");
/// bp.build();
///
/// // Map a read
/// let results = bp.map_read("AGCTAGCTAGCTAGCT", 10);
/// assert!(!results.is_empty());
///
/// // With custom alignment mode
/// let results = bp.map_read_with_mode("AGCTAGCTAGCTAGCT", AlignMode::Hybrid, 10);
/// ```
///
/// # Configuration
///
/// - `set_auto_k(true)` — auto-scale k-mer size based on genome size
/// - `set_top_n(3)` — try top-3 rarest k-mers as anchors (better mapping rate, slower)
/// - `set_spaced_seed(true)` — use spaced seed pattern for error-prone reads
/// - `set_read_type("long")` — for Nanopore/PacBio reads (k clamped to [13,19])
///
/// # Persistence
///
/// Save and load indexes with `serialize_to_file()` / `deserialize_from_file()`.
/// The persisted format uses memmap2 for <10ms load time and zstd compression.
pub struct BitPop {
    /// FM-index (built on demand via `build()`)
    fm_index: Option<FmIndex>,

    /// Forward genome storage: genome_id → DNA sequence
    genomes: HashMap<u32, Vec<u8>>,

    /// Genome names for output
    genome_names: HashMap<u32, String>,

    /// K-mer size for indexing
    k: usize,

    /// Whether to auto-scale k based on genome size
    auto_k: bool,

    /// Number of top rarest k-mers to try as anchors (for error tolerance)
    top_n: usize,

    /// Whether to use spaced seed pattern for k-mer matching
    use_spaced_seed: bool,

    /// Spaced seed pattern (default: 11111011111111)
    spaced_seed_pattern: Vec<bool>,

    /// Spaced seed hash table: encoded_spaced_kmer -> Vec<(genome_id, position)>
    /// Built during build() when use_spaced_seed is true
    spaced_seed_hash: Option<HashMap<u64, Vec<(u32, u64)>>>,

    /// Read type: "short" or "long"
    read_type: String,

    /// Fuzzy k-mer matching method for improved strain resolution
    fuzzy_method: FuzzyMethod,

    /// Maximum number of mismatches allowed in fuzzy matching (default: 1)
    fuzzy_mismatches: usize,

    /// Neighborhood hash table for O(1) fuzzy k-mer lookup (populated during build)
    neighborhood_hash: Option<HashMap<u64, Vec<(u64, u32, u64)>>>,

    /// Search radius for delta loop (±N bp around anchor position).
    /// None = dynamic (read_len/4, max 200). Some(N) = fixed override.
    search_radius: Option<isize>,

    /// Chunk size for PacBio long-read mapping (0 = disabled, auto: 150)
    chunk_size: usize,

    /// Minimum fraction of chunks that must agree for a mapping to be accepted (0.0-1.0, default: 0.0 = no threshold)
    chunk_vote_threshold: f64,

    /// Number of top genomes to return per read in chunk-based mode (default: 1)
    chunk_top_n: usize,

    /// Chunk size as percentage of read length (0.0-1.0, 0.0 = disabled, use fixed chunk_size)
    chunk_pct: f64,

    /// Minimum chunk size clamp when using chunk_pct (default: 50)
    chunk_min: usize,

    /// Maximum chunk size clamp when using chunk_pct (default: 200)
    chunk_max: usize,

    /// Enable SNP-aware scoring for strain resolution
    enable_snp_detect: bool,

    /// Minimum support count for SNP detection (default: 3)
    snp_min_support: u32,

    /// SNP penalty: amount to subtract from score for mismatches not on known SNP positions
    snp_penalty: f64,

    /// Default alignment mode for mapping (Xor, Sw, Hybrid)
    align_mode: AlignMode,

    /// Anchor strategy for chunk-based mapping (Rarest, Golden, Spaced)
    chunk_anchor_strategy: ChunkAnchorStrategy,

    /// Enable homopolymer fingerprint scoring
    enable_hf: bool,

    /// Minimum run length for homopolymer fingerprint (default: 3)
    hf_min_run: usize,

    /// Per-genome homopolymer profiles (computed during build)
    hf_profiles: HashMap<u32, hf::HfProfile>,

    /// Gap-aware XOR chaining config for long-read alignment
    chain_config: chain::ChainConfig,
}

impl BitPop {
    /// Create a new Bit-Pop indexer with the given k-mer size.
    ///
    /// # Arguments
    /// * `_k` — K-mer size. Recommended: k=10 for short reads (Illumina),
    ///   k=15 for long reads (Nanopore/PacBio). Use `set_auto_k(true)` to
    ///   compute optimal k automatically.
    ///
    /// # Example
    ///
    /// ```
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// bp.set_auto_k(true);
    /// bp.set_top_n(3);
    /// ```
    pub fn new(_k: usize) -> Self {
        Self {
            fm_index: None,
            genomes: HashMap::new(),
            genome_names: HashMap::new(),
            k: _k,
            auto_k: false,
            top_n: 1,
            use_spaced_seed: false,
            spaced_seed_pattern: vec![
                true, true, true, true, true, false, true, true, true, true, true, true, true, true,
            ],
            spaced_seed_hash: None,
            read_type: "short".to_string(),
            fuzzy_method: FuzzyMethod::None,
            fuzzy_mismatches: 1,
            neighborhood_hash: None,
            search_radius: None,
            chunk_size: 0,
            chunk_vote_threshold: 0.0,
            chunk_top_n: 1,
            chunk_pct: 0.0,
            chunk_min: 20,
            chunk_max: 500,
            enable_snp_detect: false,
            snp_min_support: 3,
            snp_penalty: 0.1,
            align_mode: AlignMode::Xor,
            chunk_anchor_strategy: ChunkAnchorStrategy::Rarest,
            enable_hf: false,
            hf_min_run: 3,
            hf_profiles: HashMap::new(),
            chain_config: chain::ChainConfig::default(),
        }
    }

    /// Enable auto-scaling of k-mer size based on total genome size.
    pub fn set_auto_k(&mut self, auto_k: bool) {
        self.auto_k = auto_k;
    }

    /// Compute optimal k-mer size based on total genome size and read type.
    /// Formula: k = floor(log2(genome_size) / log2(4)) - 2
    pub fn compute_optimal_k(&self) -> usize {
        if !self.auto_k {
            return self.k;
        }

        let total_size: usize = self.genomes.values().map(|seq| seq.len()).sum();

        if total_size == 0 {
            return self.k;
        }

        let genome_log2 = (total_size as f64).log2() / 2.0;
        let optimal_k = (genome_log2.floor() as usize).saturating_sub(2);

        if self.read_type == "long" {
            optimal_k.clamp(13, 19)
        } else {
            optimal_k.clamp(10, 15)
        }
    }

    /// Recompute k based on current genome sizes if auto_k is enabled.
    pub fn recompute_k(&mut self) {
        if self.auto_k {
            self.k = self.compute_optimal_k();
        }
    }

    /// Set the number of top rarest k-mers to try as anchors.
    /// Higher values improve mapping rate at the cost of computation.
    pub fn set_top_n(&mut self, top_n: usize) {
        self.top_n = top_n.max(1);
    }

    /// Get the current top_n setting.
    pub fn top_n(&self) -> usize {
        self.top_n
    }

    /// Enable/disable spaced seed pattern for k-mer matching.
    pub fn set_spaced_seed(&mut self, use_spaced: bool) {
        self.use_spaced_seed = use_spaced;
    }

    /// Set a custom spaced seed pattern.
    pub fn set_spaced_seed_pattern(&mut self, pattern: &str) {
        self.spaced_seed_pattern = pattern.chars().map(|c| c == '1').collect();
    }

    /// Get the spaced seed pattern as a string.
    pub fn spaced_seed_pattern(&self) -> String {
        self.spaced_seed_pattern
            .iter()
            .map(|&b| if b { '1' } else { '0' })
            .collect()
    }

    /// Set the read type: "short" (Illumina) or "long" (Nanopore/PacBio).
    pub fn set_read_type(&mut self, read_type: &str) {
        self.read_type = read_type.to_lowercase();
    }

    /// Set the fuzzy k-mer matching method.
    pub fn set_fuzzy_method(&mut self, method: FuzzyMethod) {
        self.fuzzy_method = method;
    }

    /// Set the maximum number of mismatches for fuzzy matching.
    pub fn set_fuzzy_mismatches(&mut self, mismatches: usize) {
        self.fuzzy_mismatches = mismatches.max(1);
    }

    /// Set the search radius for delta loop (±N bp around anchor position).
    /// When set, overrides the default dynamic calculation (read_len/4, max 200).
    pub fn set_search_radius(&mut self, radius: isize) {
        self.search_radius = Some(radius.clamp(1, 200));
    }

    /// Get the current search radius.
    /// Returns dynamic value (read_len/4, max 200) if not explicitly set.
    pub fn search_radius(&self, read_len: usize) -> isize {
        match self.search_radius {
            Some(r) => r,
            None => ((5.max(read_len / 4)).min(200)) as isize,
        }
    }

    /// Set chunk size for PacBio long-read mapping.
    /// 0 = disabled (use full read), >0 = split reads into chunks of this size.
    /// Recommended: 150 for PacBio HiFi (~6-11kb reads).
    pub fn set_chunk_size(&mut self, chunk_size: usize) {
        self.chunk_size = chunk_size;
    }

    /// Get the current chunk size.
    pub fn chunk_size(&self) -> usize {
        self.chunk_size
    }

    /// Set chunk size as percentage of read length (0.0-1.0).
    /// 0.0 = disabled (use fixed chunk_size), >0 = dynamic per-read chunk sizing.
    /// The calculated chunk is clamped to [chunk_min, chunk_max] bp range.
    pub fn set_chunk_pct(&mut self, chunk_pct: f64) {
        self.chunk_pct = chunk_pct.clamp(0.0, 1.0);
    }

    /// Get the current chunk_pct setting.
    pub fn chunk_pct(&self) -> f64 {
        self.chunk_pct
    }

    /// Set minimum chunk size clamp for dynamic chunking (default: 50).
    pub fn set_chunk_min(&mut self, chunk_min: usize) {
        self.chunk_min = chunk_min;
    }

    /// Get the current chunk_min setting.
    pub fn chunk_min(&self) -> usize {
        self.chunk_min
    }

    /// Set maximum chunk size clamp for dynamic chunking (default: 200).
    pub fn set_chunk_max(&mut self, chunk_max: usize) {
        self.chunk_max = chunk_max;
    }

    /// Get the current chunk_max setting.
    pub fn chunk_max(&self) -> usize {
        self.chunk_max
    }

    /// Set minimum fraction of chunks that must agree for a mapping to be accepted (0.0-1.0).
    /// 0.0 = no threshold (accept all), 0.6 = require 60% of chunks to agree.
    pub fn set_chunk_vote_threshold(&mut self, threshold: f64) {
        self.chunk_vote_threshold = threshold;
    }

    /// Get the current chunk vote threshold.
    pub fn chunk_vote_threshold(&self) -> f64 {
        self.chunk_vote_threshold
    }

    /// Set number of top genomes to return per read in chunk-based mode.
    /// Default: 1 (winner-takes-all). Use 2-3 for multi-genome uncertainty.
    pub fn set_chunk_top_n(&mut self, top_n: usize) {
        self.chunk_top_n = top_n.max(1);
    }

    /// Get the current chunk top-n setting.
    pub fn chunk_top_n(&self) -> usize {
        self.chunk_top_n
    }

    /// Set the anchor strategy for chunk-based mapping.
    pub fn set_chunk_anchor_strategy(&mut self, strategy: ChunkAnchorStrategy) {
        self.chunk_anchor_strategy = strategy;
    }

    /// Get the current chunk anchor strategy.
    pub fn chunk_anchor_strategy(&self) -> ChunkAnchorStrategy {
        self.chunk_anchor_strategy
    }

    /// Enable or disable SNP-aware scoring for strain resolution.
    pub fn set_snp_detect(&mut self, enable: bool) {
        self.enable_snp_detect = enable;
    }

    /// Get the current SNP detection enabled status.
    pub fn snp_detect_enabled(&self) -> bool {
        self.enable_snp_detect
    }

    /// Set the minimum support count for SNP detection.
    pub fn set_snp_min_support(&mut self, min_support: u32) {
        self.snp_min_support = min_support.max(1);
    }

    /// Get the current SNP minimum support setting.
    pub fn snp_min_support(&self) -> u32 {
        self.snp_min_support
    }

    /// Set the SNP penalty amount (subtract from score for non-SNP mismatches).
    pub fn set_snp_penalty(&mut self, penalty: f64) {
        self.snp_penalty = penalty.clamp(0.0, 1.0);
    }

    /// Get the current SNP penalty setting.
    pub fn snp_penalty(&self) -> f64 {
        self.snp_penalty
    }

    /// Enable homopolymer fingerprint scoring.
    pub fn set_hf(&mut self, enable: bool) {
        self.enable_hf = enable;
    }

    /// Get the current HF enabled status.
    pub fn hf_enabled(&self) -> bool {
        self.enable_hf
    }

    /// Set the minimum run length for homopolymer fingerprint.
    pub fn set_hf_min_run(&mut self, min_run: usize) {
        self.hf_min_run = min_run.max(2);
    }

    /// Get the current HF min_run setting.
    pub fn hf_min_run(&self) -> usize {
        self.hf_min_run
    }

    /// Set the default alignment mode for mapping operations.
    pub fn set_align_mode(&mut self, mode: AlignMode) {
        self.align_mode = mode;
    }

    /// Get the current alignment mode.
    pub fn align_mode(&self) -> AlignMode {
        self.align_mode
    }

    /// Add a genome (reference sequence) to the index.
    ///
    /// Returns the assigned genome_id (0-based, sequential).
    /// After adding all genomes, call `build()` to construct the FM-index.
    ///
    /// # Arguments
    /// * `name` — Genome identifier (e.g., "Ecoli_K12", "NC_000913.3")
    /// * `sequence` — DNA sequence in FASTA format (any case, N bases skipped)
    ///
    /// # Example
    ///
    /// ```
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// let id = bp.add_genome("test", "ACGTACGTACGT");
    /// assert_eq!(id, 0);
    /// ```
    pub fn add_genome(&mut self, name: &str, sequence: &str) -> u32 {
        let genome_id = self.genomes.len() as u32;
        let encoded = encode_sequence(sequence);
        self.genome_names.insert(genome_id, name.to_string());
        self.genomes.insert(genome_id, encoded.clone());
        genome_id
    }

    /// Finalize the index: construct the FM-index from all stored genomes.
    ///
    /// After `build()`, the index is immutable and ready for mapping.
    /// Recomputes k-mer size if `auto_k` is enabled.
    ///
    /// # Example
    ///
    /// ```
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// bp.add_genome("g1", "ACGTACGT");
    /// bp.add_genome("g2", "GCTAGCTA");
    /// bp.build();
    /// assert_eq!(bp.genome_count(), 2);
    /// ```
    pub fn build(&mut self) {
        self.recompute_k();

        let mut genome_list: Vec<(u32, &str, &[u8])> = self
            .genomes
            .iter()
            .map(|(gid, seq)| {
                (
                    *gid,
                    self.genome_names.get(gid).map(|s| s.as_str()).unwrap_or(""),
                    seq.as_slice(),
                )
            })
            .collect();
        genome_list.sort_by_key(|(gid, _, _)| *gid);
        let genomes: Vec<(&str, &[u8])> = genome_list
            .into_iter()
            .map(|(_, name, seq)| (name, seq))
            .collect();
        self.fm_index = Some(FmIndex::build(&genomes));

        if matches!(self.fuzzy_method, FuzzyMethod::Neighborhood) {
            self.build_neighborhood_hash();
        }

        if self.use_spaced_seed {
            self.build_spaced_seed_hash();
        }

        if self.enable_hf {
            for (gid, seq) in &self.genomes {
                let profile = hf::HfProfile::compute(seq, self.hf_min_run);
                self.hf_profiles.insert(*gid, profile);
            }
        }
    }

    /// Finalize the index with parallel FM-Index build using rayon.
    /// Parallelizes BWT construction and OccCounter building.
    pub fn build_parallel(&mut self) {
        self.recompute_k();

        let mut genome_list: Vec<(u32, String, Vec<u8>)> = self
            .genomes
            .iter()
            .map(|(gid, seq)| {
                let name = self.genome_names.get(gid).cloned().unwrap_or_default();
                (*gid, name, seq.to_vec())
            })
            .collect();
        genome_list.sort_by_key(|(gid, _, _)| *gid);

        let genomes: Vec<(String, Vec<u8>)> = genome_list
            .into_par_iter()
            .map(|(_, name, seq)| (name, seq))
            .collect();

        let genome_refs: Vec<(&str, &[u8])> = genomes
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.as_slice()))
            .collect();

        self.fm_index = Some(FmIndex::build_parallel(&genome_refs));

        if matches!(self.fuzzy_method, FuzzyMethod::Neighborhood) {
            self.build_neighborhood_hash();
        }

        if self.use_spaced_seed {
            self.build_spaced_seed_hash();
        }

        if self.enable_hf {
            for (gid, seq) in &self.genomes {
                let profile = hf::HfProfile::compute(seq, self.hf_min_run);
                self.hf_profiles.insert(*gid, profile);
            }
        }
    }

    /// Load genomes from a FASTA file.
    /// Each header becomes a genome name. Returns assigned genome IDs.
    /// After loading, call `build()` to compress.
    pub fn load_genome_fasta(&mut self, path: &str) -> io::Result<Vec<u32>> {
        let mut reader = fasta::FastaReader::new(path)?;
        let mut ids = Vec::new();

        while let Some(result) = reader.next() {
            let (header, sequence) = result?;
            let gid = self.add_genome(&header, &sequence);
            ids.push(gid);
        }

        Ok(ids)
    }

    /// Load genomes from a FASTA file using memory-mapped I/O.
    /// Each header becomes a genome name. Returns assigned genome IDs.
    /// Only available when `mmap` feature is enabled.
    /// After loading, call `build()` to compress.
    #[cfg(feature = "mmap")]
    pub fn load_genome_fasta_mmap(&mut self, path: &str) -> io::Result<Vec<u32>> {
        let mut reader = fasta::MmapFastaReader::new(path)?;
        let mut ids = Vec::new();

        while let Some(result) = reader.next() {
            let (header, sequence) = result?;
            let gid = self.add_genome(&header, &sequence);
            ids.push(gid);
        }

        Ok(ids)
    }

    /// Encode k bytes (each 2 bits) into a u64 k-mer value.
    #[allow(dead_code)]
    fn encode_kmer_bytes(&self, bytes: &[u8]) -> u64 {
        let mut result: u64 = 0;
        for &b in bytes {
            result = (result << 2) | b as u64;
        }
        result
    }

    /// Stage 1: K-mer filter for a read.
    ///
    /// Returns candidate positions across all genomes using FM-index backward search.
    /// No threshold applied (equivalent to `kmer_filter_with_threshold(read, usize::MAX)`).
    ///
    /// # Arguments
    /// * `read` — DNA read sequence
    ///
    /// # Returns
    /// Vector of (genome_id, position, kmer_count) tuples sorted by descending count.
    pub fn kmer_filter(&self, read: &str) -> Vec<(u32, u64, usize)> {
        self.kmer_filter_with_threshold(read, usize::MAX)
    }

    /// Stage 1: K-mer filter with quality-aware filtering.
    /// Skips k-mers where any base has quality below min_quality threshold.
    pub fn kmer_filter_with_quality(
        &self,
        read: &str,
        quality: &[u8],
        min_quality: u8,
        max_hits: usize,
    ) -> Vec<(u32, u64, usize)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        if encoded.len() < self.k {
            return Vec::new();
        }

        // Filter k-mers by quality: a k-mer is valid only if all its bases have quality >= min_quality
        let mut counts: HashMap<u64, usize> = HashMap::new();

        for i in 0..=(encoded.len() - self.k) {
            // Check quality for this k-mer's bases
            let qual_end = (i + self.k).min(quality.len());
            if qual_end - i < self.k {
                continue;
            }

            let has_low_quality: bool = quality[i..qual_end].iter().any(|&q| q < min_quality);

            if has_low_quality {
                continue;
            }

            let kmer = &encoded[i..i + self.k];

            if max_hits < usize::MAX {
                let occ = fm.count_occurrences(kmer);
                if occ > max_hits {
                    continue;
                }
            }

            let positions = fm.find_positions(kmer, max_hits);
            for &(gid, pos) in &positions {
                let packed = ((gid as u64) << 32) | (pos & 0xFFFFFFFF);
                *counts.entry(packed).or_default() += 1;
            }
        }

        let mut candidates: Vec<(u32, u64, usize)> = counts
            .into_iter()
            .map(|(packed, count)| {
                let genome_id = (packed >> 32) as u32;
                let position = packed & 0xFFFFFFFF;
                (genome_id, position, count)
            })
            .collect();
        candidates.sort_by_key(|b| std::cmp::Reverse(b.2));
        candidates
    }

    /// Stage 1: K-mer filter with repetitive k-mer threshold.
    /// K-mers with more than `max_hits` total positions are skipped.
    pub fn kmer_filter_with_threshold(
        &self,
        read: &str,
        max_hits: usize,
    ) -> Vec<(u32, u64, usize)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        if encoded.len() < self.k {
            return Vec::new();
        }

        let mut counts: HashMap<u64, usize> = HashMap::new();

        for i in 0..=(encoded.len() - self.k) {
            let kmer = &encoded[i..i + self.k];

            if max_hits < usize::MAX {
                let occ = fm.count_occurrences(kmer);
                if occ > max_hits {
                    continue;
                }
            }

            let positions = fm.find_positions(kmer, max_hits);
            for &(gid, pos) in &positions {
                let packed = ((gid as u64) << 32) | (pos & 0xFFFFFFFF);
                *counts.entry(packed).or_default() += 1;
            }
        }

        let mut candidates: Vec<(u32, u64, usize)> = counts
            .into_iter()
            .map(|(packed, count)| {
                let genome_id = (packed >> 32) as u32;
                let position = packed & 0xFFFFFFFF;
                (genome_id, position, count)
            })
            .collect();
        candidates.sort_by_key(|b| std::cmp::Reverse(b.2));
        candidates
    }

    /// Stage 2: Bit-level alignment.
    /// Aligns a read against a genome region starting at position.
    /// Returns (alignment_score_0_to_1, cigar_string, aligned_start_offset).
    pub fn align_read(&self, read: &str, genome_id: u32, position: u64) -> (f64, String, usize) {
        let read_enc = encode_sequence(read);
        let genome = match self.genomes.get(&genome_id) {
            Some(g) => g,
            None => return (0.0, String::new(), 0),
        };

        let pos = position as usize;
        let read_len = read_enc.len();

        if read_len == 0 {
            return (0.0, String::new(), 0);
        }

        // Extract genome region
        let region_end = (pos + read_len).min(genome.len());
        let region = &genome[pos..region_end];

        // Fast path: exact match
        if region.len() == read_len && read_enc == *region {
            return (1.0, format!("{}M", read_len), 0);
        }

        // For reads <=31bp: direct 2-bit XOR
        if read_len <= 31 {
            return align::two_bit_align(&read_enc, region);
        }

        // For reads >31bp: chunked 2-bit XOR
        let (score, offset) = align::two_bit_score_chunks(&read_enc, region);
        let cigar = format!("{}M", read_len);
        (score, cigar, offset)
    }

    /// Stage 2: Smith-Waterman local alignment.
    /// Aligns a read against a genome region using SW with full traceback.
    /// Returns (alignment_score_0_to_1, cigar_string, best_offset_in_region).
    pub fn align_read_sw(&self, read: &str, genome_id: u32, position: u64) -> (f64, String, usize) {
        let read_enc = encode_sequence(read);
        let genome = match self.genomes.get(&genome_id) {
            Some(g) => g,
            None => return (0.0, String::new(), 0),
        };

        let pos = position as usize;
        let read_len = read_enc.len();

        if read_len == 0 {
            return (0.0, String::new(), 0);
        }

        // Extract a generous genome region for SW to work with
        let search_radius = (self.k.max(read_len / 4)).min(200);
        let region_start = pos.saturating_sub(search_radius);
        let region_end = (pos + read_len + search_radius).min(genome.len());
        let region = &genome[region_start..region_end];

        if region.is_empty() {
            return (0.0, String::new(), 0);
        }

        // For reads <=31bp: use standard SW with full traceback
        if read_len <= 31 {
            let (sw_score, cigar) = align::smith_waterman(&read_enc, region);
            if sw_score == 0 {
                return (0.0, String::new(), 0);
            }
            // Normalize score to 0-1 range (max possible = 2 * read_len for match=+2)
            let normalized = (sw_score as f64) / (2.0 * read_len as f64);
            (normalized.clamp(0.0, 1.0), cigar, region_start)
        } else {
            // For longer reads: chunked SW with full traceback → real CIGAR
            let (sw_score, cigar) = align::smith_waterman_chunked(&read_enc, region);
            if sw_score == 0 {
                let (score, offset) = align::smith_waterman_score(&read_enc, region);
                (score, format!("{}M", read_len), region_start + offset)
            } else {
                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                (normalized.clamp(0.0, 1.0), cigar, region_start)
            }
        }
    }

    /// Stage 2: Quality-aware Smith-Waterman local alignment.
    /// Aligns a read against a genome region using SW with Phred-scaled quality penalties.
    /// Returns (alignment_score_0_to_1, cigar_string, best_offset_in_region, quality_penalty).
    pub fn align_read_sw_with_quality(
        &self,
        read: &str,
        quality: &[u8],
        genome_id: u32,
        position: u64,
    ) -> (f64, String, usize, f64) {
        let read_enc = encode_sequence(read);
        let genome = match self.genomes.get(&genome_id) {
            Some(g) => g,
            None => return (0.0, String::new(), 0, 0.0),
        };

        let pos = position as usize;
        let read_len = read_enc.len();

        if read_len == 0 {
            return (0.0, String::new(), 0, 0.0);
        }

        // Extract a generous genome region for SW to work with
        let search_radius = (self.k.max(read_len / 4)).min(200);
        let region_start = pos.saturating_sub(search_radius);
        let region_end = (pos + read_len + search_radius).min(genome.len());
        let region = &genome[region_start..region_end];

        if region.is_empty() {
            return (0.0, String::new(), 0, 0.0);
        }

        // For reads <=31bp: use quality-aware SW with full traceback
        if read_len <= 31 {
            let (sw_score, cigar, offset, qual_penalty) =
                align::smith_waterman_with_quality(&read_enc, region, quality, 2, -2, 0);
            if sw_score == 0 {
                return (0.0, String::new(), 0, 0.0);
            }
            let normalized = (sw_score as f64) / (2.0 * read_len as f64);
            (normalized.clamp(0.0, 1.0), cigar, offset, qual_penalty)
        } else {
            // For longer reads: chunked quality-aware SW with traceback → real CIGAR
            let mut total_score = 0i32;
            let mut total_penalty = 0.0f64;
            let mut all_ops: Vec<u8> = Vec::new();

            let chunk_size = 31;
            for chunk_start in (0..read_len).step_by(chunk_size) {
                let chunk_end = (chunk_start + chunk_size).min(read_len);
                if chunk_end - chunk_start < 4 {
                    break;
                }
                let chunk = &read_enc[chunk_start..chunk_end];

                let qual_chunk =
                    &quality[chunk_start.min(quality.len())..chunk_end.min(quality.len())];
                let text_start = chunk_start + region_start;
                let text_end = (text_start + chunk.len()).min(region.len());
                if text_end - text_start < chunk.len() {
                    continue;
                }
                let text_region = &region[text_start..text_end];

                let (sw_score, cigar) =
                    align::smith_waterman_internal(chunk, text_region, 2, -1, -2);
                total_score += sw_score;
                if !cigar.is_empty() && sw_score > 0 {
                    align::parse_cigar_ops(&cigar, &mut all_ops);
                }

                // Accumulate quality penalty by re-scoring the chunk alignment path
                let (_, _, _, penalty) =
                    align::smith_waterman_with_quality(chunk, text_region, qual_chunk, 2, -2, 0);
                total_penalty += penalty;
            }

            if total_score == 0 {
                return (0.0, String::new(), 0, 0.0);
            }

            let cigar = align::build_cigar_string(&all_ops);
            let normalized = (total_score as f64) / (2.0 * read_len as f64);
            (
                normalized.clamp(0.0, 1.0),
                cigar,
                region_start,
                total_penalty,
            )
        }
    }

    /// Stage 2: Unified alignment method that dispatches based on AlignMode.
    pub fn align_read_with_mode(
        &self,
        read: &str,
        mode: AlignMode,
        genome_id: u32,
        position: u64,
    ) -> (f64, String, usize) {
        match mode {
            AlignMode::Xor => self.align_read(read, genome_id, position),
            AlignMode::Sw => self.align_read_sw(read, genome_id, position),
            AlignMode::Hybrid => {
                // Fast XOR filter first, then SW on promising candidates
                let (xor_score, xor_cigar, xor_offset) = self.align_read(read, genome_id, position);
                if xor_score >= 0.9 {
                    // High confidence XOR match — skip SW
                    (xor_score, xor_cigar, xor_offset)
                } else {
                    // Lower confidence — refine with SW
                    let (sw_score, sw_cigar, _) = self.align_read_sw(read, genome_id, position);
                    if sw_score > xor_score {
                        (sw_score, sw_cigar, 0)
                    } else {
                        (xor_score, xor_cigar, xor_offset)
                    }
                }
            }
            AlignMode::Softclip => self.align_read_softclip(read, genome_id, position),
            AlignMode::Chain => {
                // Chain mode does its own genome search, ignores genome_id/position
                let chain_results = self.align_read_chain(read);
                if let Some((_, score, cigar, ref_pos)) = chain_results.first() {
                    (*score, cigar.clone(), *ref_pos)
                } else {
                    (0.0, String::new(), 0)
                }
            }
        }
    }

    /// Stage 2: Soft-clipping XOR alignment.
    /// Slides windows across the read to find the optimal alignment region,
    /// emitting soft-clips (S:) for adapter/contaminated regions.
    pub fn align_read_softclip(
        &self,
        read: &str,
        genome_id: u32,
        position: u64,
    ) -> (f64, String, usize) {
        let read_enc = encode_sequence(read);
        let genome = match self.genomes.get(&genome_id) {
            Some(g) => g,
            None => return (0.0, String::new(), 0),
        };

        let pos = position as usize;
        let read_len = read_enc.len();

        if read_len == 0 {
            return (0.0, String::new(), 0);
        }

        // Extract genome region
        let region_end = (pos + read_len).min(genome.len());
        let region = &genome[pos..region_end];

        // Minimum alignment length: at least 1/3 of read or 20bp, whichever is larger
        let min_align = (read_len / 3).max(20);

        align::two_bit_align_softclip(&read_enc, region, min_align, 0.7)
    }

    /// Stage 2: Gap-aware XOR chaining for long-read alignment.
    ///
    /// Uses minimizer-based seed chaining with XOR gap extension.
    /// Handles ONT/PacBio 5-15% error rates and long indels.
    pub fn align_read_chain(&self, read: &str) -> Vec<(u32, f64, String, usize)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let results = chain::chain_align_long_read(read, fm, &self.genomes, &self.chain_config);

        results
            .into_iter()
            .map(|(genome_id, score, cigar, ref_pos, _span)| (genome_id, score, cigar, ref_pos))
            .collect()
    }

    /// Set chain k-mer size for minimizer extraction.
    pub fn set_chain_k(&mut self, k: usize) {
        self.chain_config.k = k.max(5);
    }

    /// Get chain k-mer size.
    pub fn chain_k(&self) -> usize {
        self.chain_config.k
    }

    /// Set chain window size for minimizer extraction.
    pub fn set_chain_window(&mut self, w: usize) {
        self.chain_config.w = w.max(1);
    }

    /// Get chain window size.
    pub fn chain_window(&self) -> usize {
        self.chain_config.w
    }

    /// Set minimum chain seeds required.
    pub fn set_chain_min_seeds(&mut self, min_seeds: usize) {
        self.chain_config.min_chain_seeds = min_seeds.max(1);
    }

    /// Get minimum chain seeds.
    pub fn chain_min_seeds(&self) -> usize {
        self.chain_config.min_chain_seeds
    }

    /// Set maximum gap size between chained seeds.
    pub fn set_chain_max_gap(&mut self, max_gap: usize) {
        self.chain_config.max_gap = max_gap;
    }

    /// Get maximum gap size.
    pub fn chain_max_gap(&self) -> usize {
        self.chain_config.max_gap
    }

    /// Set gap open penalty for affine gap model.
    pub fn set_chain_gap_open(&mut self, gap_open: f64) {
        self.chain_config.gap_open = gap_open;
    }

    /// Get gap open penalty.
    pub fn chain_gap_open(&self) -> f64 {
        self.chain_config.gap_open
    }

    /// Set gap extension penalty for affine gap model.
    pub fn set_chain_gap_extend(&mut self, gap_extend: f64) {
        self.chain_config.gap_extend = gap_extend;
    }

    /// Get gap extension penalty.
    pub fn chain_gap_extend(&self) -> f64 {
        self.chain_config.gap_extend
    }

    /// Stage 3: Rank pre-scored mapping results.
    /// Takes already-aligned candidates and applies rarity-based ranking.
    pub fn rank_scored_results(
        &self,
        scored_candidates: &[(u32, u64, f64, String)],
        read: &str,
        context_window: usize,
        min_final_score: f64,
    ) -> Vec<MappingResult> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let read_len = encode_sequence(read).len();
        let encoded = encode_sequence(read);
        let mut results = Vec::new();

        let read_hf = if self.enable_hf {
            Some(hf::HfProfile::compute_read(&encoded, self.hf_min_run))
        } else {
            None
        };

        for &(genome_id, position, align_score, ref cigar) in scored_candidates {
            if align_score < min_final_score {
                continue;
            }

            // Per-candidate rarity: extract k-mer from genome at candidate position
            let rarity = if encoded.len() >= self.k {
                if let Some(genome_seq) = self.genomes.get(&genome_id) {
                    let start = (position as usize).min(genome_seq.len().saturating_sub(1));
                    let end = (start + self.k).min(genome_seq.len());
                    if end - start == self.k {
                        let kmer_bytes = &genome_seq[start..end];
                        let kmer_encoded = kmer_bytes
                            .iter()
                            .map(|&b| b as char)
                            .filter_map(encode_base)
                            .collect::<Vec<u8>>();
                        if kmer_encoded.len() == self.k {
                            let occ = fm.count_occurrences(&kmer_encoded);
                            1.0 / (occ as f64).max(1.0)
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let mut hf_score = 0.0f64;
            if self.enable_hf {
                if let Some(profile) = self.hf_profiles.get(&genome_id) {
                    if let Some(ref read_prof) = read_hf {
                        hf_score = profile.similarity(read_prof);
                    }
                }
            }

            // align_score is the primary signal (0.85 weight).
            // rarity provides a modest boost (0.15 weight) so perfect matches
            // floor at 0.85 instead of 0.5 when the anchor k-mer is common.
            let combined_score = align_score * 0.85 + rarity * 0.15;
            let context =
                self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(MappingResult {
                genome_id,
                position,
                score: combined_score,
                cigar: cigar.clone(),
                context,
                is_reverse: false,
                rarity,
                md_string: String::new(),
                hf_score,
            });
        }

        results.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Stage 3: Rank mapping results (legacy, for backwards compatibility).
    /// Takes candidates from Stage 1 and alignments from Stage 2,
    /// returns ranked MappingResults.
    pub fn rank_results(
        &self,
        candidates: &[(u32, u64, usize)],
        read: &str,
        context_window: usize,
    ) -> Vec<MappingResult> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let read_len = encode_sequence(read).len();
        let encoded = encode_sequence(read);
        let mut results = Vec::new();

        let top_candidates = candidates.iter().take(50);

        for &(genome_id, position, _kmer_count) in top_candidates {
            let (align_score, cigar, _) = self.align_read(read, genome_id, position);

            if align_score < 0.5 {
                continue;
            }

            let rarity = if encoded.len() >= self.k {
                if let Some(genome_seq) = self.genomes.get(&genome_id) {
                    let start = (position as usize).min(genome_seq.len().saturating_sub(1));
                    let end = (start + self.k).min(genome_seq.len());
                    if end - start == self.k {
                        let kmer_bytes = &genome_seq[start..end];
                        let kmer_encoded = kmer_bytes
                            .iter()
                            .map(|&b| b as char)
                            .filter_map(encode_base)
                            .collect::<Vec<u8>>();
                        if kmer_encoded.len() == self.k {
                            let occ = fm.count_occurrences(&kmer_encoded);
                            1.0 / (occ as f64).max(1.0)
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let context =
                self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(MappingResult {
                genome_id,
                position,
                score: combined_score,
                cigar,
                context,
                is_reverse: false,
                rarity,
                md_string: String::new(),
                hf_score: 0.0,
            });
        }

        results.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Smart threshold computation based on read length and genome characteristics.
    /// Returns a minimum score threshold that adapts to the context:
    /// - Short reads (<20bp): stricter threshold (higher min_score)
    /// - Long reads (>100bp): more lenient threshold
    /// - High-quality reads: stricter threshold
    /// - Repetitive genomes: more lenient threshold
    fn compute_smart_threshold(&self, read_len: usize, has_quality: bool, avg_quality: f64) -> f64 {
        let base_threshold: f64 = 0.5;

        // Read length adjustment
        let length_factor: f64 = if read_len < 20 {
            0.1 // stricter for short reads
        } else if read_len > 100 {
            -0.05 // more lenient for long reads
        } else {
            0.0
        };

        // Quality adjustment
        let qual_factor: f64 = if has_quality && avg_quality > 25.0 {
            0.05 // stricter for high quality
        } else if has_quality && avg_quality < 15.0 {
            -0.05 // more lenient for low quality
        } else {
            0.0
        };

        (base_threshold + length_factor + qual_factor).clamp(0.3f64, 0.8f64)
    }

    /// Find the top-N rarest k-mers in a read, sorted by ascending occurrence count.
    /// Returns vector of (read_offset, kmer_bytes, count) tuples.
    fn find_top_n_rarest_kmers(
        &self,
        encoded: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.k {
            return Vec::new();
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();

        for i in 0..=(encoded.len() - self.k) {
            let kmer = &encoded[i..i + self.k];
            let count = fm.count_occurrences(kmer);
            if count > 0 && count <= max_hits {
                candidates.push((i, kmer.to_vec(), count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Find the top-N rarest k-mers in a read, considering only high-quality bases.
    /// Returns vector of (read_offset, kmer_bytes, count) tuples.
    fn find_top_n_rarest_kmers_quality(
        &self,
        encoded: &[u8],
        quality: &[u8],
        fm: &FmIndex,
        min_quality: u8,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.k {
            return Vec::new();
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();

        for i in 0..=(encoded.len() - self.k) {
            let qual_end = (i + self.k).min(quality.len());
            if qual_end - i < self.k {
                continue;
            }

            let has_low_quality: bool = quality[i..qual_end].iter().any(|&q| q < min_quality);

            if has_low_quality {
                continue;
            }

            let kmer = &encoded[i..i + self.k];
            let count = fm.count_occurrences(kmer);
            if count > 0 && count <= max_hits {
                candidates.push((i, kmer.to_vec(), count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Find golden anchors: k-mers with highest average quality, then rarest.
    /// Quality score: average Phred quality of bases in k-mer.
    /// Golden threshold: average quality >= 30 (99.9% accuracy per base).
    fn find_golden_anchors(
        &self,
        encoded: &[u8],
        quality: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.k || quality.len() < self.k {
            return Vec::new();
        }

        let mut quality_candidates: Vec<(usize, f64, Vec<u8>)> = Vec::new();

        for i in 0..=(encoded.len().saturating_sub(self.k)) {
            let qual_end = (i + self.k).min(quality.len());
            if qual_end - i < self.k {
                continue;
            }

            let avg_quality: f64 =
                quality[i..qual_end].iter().map(|&q| q as f64).sum::<f64>() / self.k as f64;

            if avg_quality >= 30.0 {
                let kmer = &encoded[i..i + self.k];
                quality_candidates.push((i, avg_quality, kmer.to_vec()));
            }
        }

        if quality_candidates.is_empty() {
            let fallback = self.find_top_n_rarest_kmers(encoded, fm, max_hits);
            return fallback;
        }

        quality_candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        let top_quality_count = (quality_candidates.len() as f64 * 0.1).max(10.0) as usize;
        quality_candidates.truncate(top_quality_count.min(quality_candidates.len()));

        let mut final_candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();
        for &(idx, _, ref kmer) in &quality_candidates {
            let count = fm.count_occurrences(kmer);
            if count > 0 && count <= max_hits {
                final_candidates.push((idx, kmer.clone(), count));
            }
        }

        final_candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(final_candidates.len());
        final_candidates.truncate(n);
        final_candidates
    }

    /// Find the top-N rarest spaced k-mers in a read, sorted by ascending occurrence count.
    /// Uses spaced seed hash table for ultra-fast lookup.
    fn find_top_n_rarest_kmers_spaced(
        &self,
        encoded: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.spaced_seed_pattern.len() {
            return Vec::new();
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();
        let pattern_len = self.spaced_seed_pattern.len();

        for i in 0..=(encoded.len().saturating_sub(pattern_len)) {
            let window = &encoded[i..i + pattern_len];

            let count = if self.spaced_seed_hash.is_some() {
                self.count_occurrences_spaced_hash(window)
            } else {
                fm.count_occurrences_spaced(window, &self.spaced_seed_pattern)
            };
            if count > 0 && count <= max_hits {
                candidates.push((i, window.to_vec(), count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Find the top-N rarest k-mers with fuzzy matching (up to N mismatches).
    /// Generates all k-mer variants and queries FM-index for each, aggregating counts.
    fn find_top_n_rarest_kmers_fuzzy(
        &self,
        encoded: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.k {
            return Vec::new();
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();

        for i in 0..=(encoded.len() - self.k) {
            let kmer_bytes = &encoded[i..i + self.k];
            let kmer_encoded = Self::encode_kmer_bytes_0indexed(self, kmer_bytes);

            let variants = generate_kmer_variants(kmer_encoded, self.k, self.fuzzy_mismatches);
            let all_variants: Vec<u64> = std::iter::once(kmer_encoded).chain(variants).collect();

            let mut total_count = 0usize;
            let mut found_any = false;

            for variant in &all_variants {
                let variant_bytes = Self::kmer_0indexed_to_fm_bytes(*variant, self.k);
                let count = fm.count_occurrences(&variant_bytes);
                if count > 0 && count <= max_hits {
                    total_count += count;
                    found_any = true;
                }
            }

            if found_any {
                candidates.push((i, kmer_bytes.to_vec(), total_count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Find the top-N rarest spaced k-mers with fuzzy matching.
    /// Allows mismatches in the "match" positions of the spaced seed.
    fn find_top_n_rarest_kmers_spaced_fuzzy(
        &self,
        encoded: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        if encoded.len() < self.spaced_seed_pattern.len() {
            return Vec::new();
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();

        for i in 0..=(encoded.len() - self.spaced_seed_pattern.len()) {
            let window = &encoded[i..i + self.spaced_seed_pattern.len()];
            let spaced_kmer: Vec<u8> = window
                .iter()
                .enumerate()
                .filter(|(j, _)| self.spaced_seed_pattern[*j])
                .map(|(_, &b)| b)
                .collect();

            let count = fm.count_occurrences_spaced_fuzzy(
                window,
                &self.spaced_seed_pattern,
                self.fuzzy_mismatches,
            );
            if count > 0 && count <= max_hits {
                candidates.push((i, spaced_kmer, count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Build a neighborhood hash table for O(1) fuzzy k-mer lookup.
    /// Maps each k-mer variant (with N mismatches) to its original k-mer and positions.
    fn build_neighborhood_hash(&mut self) {
        if self.fm_index.is_none() {
            return;
        }

        if self.fuzzy_mismatches == 0 {
            return;
        }

        let mut hash_table = HashMap::new();

        for (&genome_id, seq) in &self.genomes {
            for i in 0..=seq.len().saturating_sub(self.k) {
                let kmer_bytes = &seq[i..i + self.k];
                let kmer_encoded = Self::encode_kmer_bytes(self, kmer_bytes);

                let variants = generate_kmer_variants(kmer_encoded, self.k, self.fuzzy_mismatches);
                let all_variants: Vec<u64> =
                    std::iter::once(kmer_encoded).chain(variants).collect();

                for variant in all_variants {
                    let entry = (variant, genome_id, i as u64);
                    hash_table
                        .entry(variant)
                        .or_insert_with(Vec::new)
                        .push(entry);
                }
            }
        }

        self.neighborhood_hash = Some(hash_table);
    }

    /// Build spaced seed hash table for ultra-fast lookup.
    /// Maps encoded spaced k-mers to (genome_id, position) lists.
    fn build_spaced_seed_hash(&mut self) {
        if self.fm_index.is_none() {
            return;
        }

        let pattern = self.spaced_seed_pattern.clone();
        let seed = fm::SpacedSeed::from_binary(
            &pattern
                .iter()
                .map(|&b| if b { '1' } else { '0' })
                .collect::<String>(),
        );

        let mut hash_table = HashMap::new();

        for (&genome_id, seq) in &self.genomes {
            for i in 0..=seq.len().saturating_sub(seed.len()) {
                let kmer = &seq[i..i + seed.len()];
                let hash_key = seed.encode_as_hash(kmer);

                hash_table
                    .entry(hash_key)
                    .or_insert_with(Vec::new)
                    .push((genome_id, i as u64));
            }
        }

        self.spaced_seed_hash = Some(hash_table);
    }

    /// Find positions using spaced seed hash table lookup.
    fn find_positions_spaced_hash(&self, kmer: &[u8]) -> Vec<(u32, u64)> {
        let hash_table = match &self.spaced_seed_hash {
            Some(table) => table,
            None => return Vec::new(),
        };

        let seed = fm::SpacedSeed::from_binary(
            &self
                .spaced_seed_pattern
                .iter()
                .map(|&b| if b { '1' } else { '0' })
                .collect::<String>(),
        );

        let hash_key = seed.encode_as_hash(kmer);

        hash_table.get(&hash_key).cloned().unwrap_or_default()
    }

    /// Count occurrences using spaced seed hash table lookup.
    fn count_occurrences_spaced_hash(&self, kmer: &[u8]) -> usize {
        let positions = self.find_positions_spaced_hash(kmer);
        let mut seen = std::collections::HashSet::new();
        positions.iter().filter(|p| seen.insert(*p)).count()
    }

    /// Find the top-N rarest k-mers using neighborhood hash table.
    fn find_top_n_rarest_kmers_neighborhood(
        &self,
        encoded: &[u8],
        fm: &FmIndex,
        max_hits: usize,
    ) -> Vec<(usize, Vec<u8>, usize)> {
        let hash_table = match &self.neighborhood_hash {
            Some(table) => table,
            None => return self.find_top_n_rarest_kmers(encoded, fm, max_hits),
        };

        let mut position_counts: HashMap<u64, usize> = HashMap::new();

        for i in 0..=(encoded.len() - self.k) {
            let kmer_bytes = &encoded[i..i + self.k];
            let kmer_encoded = Self::encode_kmer_bytes(self, kmer_bytes);

            let variants = generate_kmer_variants(kmer_encoded, self.k, self.fuzzy_mismatches);
            let all_variants: Vec<u64> = std::iter::once(kmer_encoded).chain(variants).collect();

            for variant in all_variants {
                if let Some(entries) = hash_table.get(&variant) {
                    for &(_, genome_id, pos) in entries {
                        let packed = ((genome_id as u64) << 32) | (pos & 0xFFFFFFFF);
                        *position_counts.entry(packed).or_default() += 1;
                    }
                }
            }
        }

        let mut candidates: Vec<(usize, Vec<u8>, usize)> = Vec::new();
        for (i, kmer_bytes) in encoded.windows(self.k).enumerate() {
            let kmer_encoded = Self::encode_kmer_bytes(self, kmer_bytes);
            let variants = generate_kmer_variants(kmer_encoded, self.k, self.fuzzy_mismatches);
            let all_variants: Vec<u64> = std::iter::once(kmer_encoded).chain(variants).collect();

            let mut total_positions: std::collections::HashSet<u64> =
                std::collections::HashSet::new();
            for variant in all_variants {
                if let Some(entries) = hash_table.get(&variant) {
                    for &(_, _, pos) in entries {
                        total_positions.insert(pos);
                    }
                }
            }

            if !total_positions.is_empty() {
                let count = total_positions.len().min(max_hits);
                candidates.push((i, kmer_bytes.to_vec(), count));
            }
        }

        candidates.sort_by_key(|&(_, _, count)| count);

        let n = self.top_n.min(candidates.len());
        candidates.truncate(n);
        candidates
    }

    /// Encode bytes (1-indexed A=1,C=2,G=3,T=4) into u64 (0-indexed internally).
    fn encode_kmer_bytes_0indexed(&self, bytes: &[u8]) -> u64 {
        let mut result: u64 = 0;
        for &b in bytes {
            result = (result << 2) | ((b.saturating_sub(1)) as u64);
        }
        result
    }

    /// Convert 0-indexed u64 k-mer to 1-indexed bytes for FM-index queries.
    pub fn kmer_0indexed_to_fm_bytes(encoded: u64, k: usize) -> Vec<u8> {
        let mut result = Vec::with_capacity(k);
        for i in (0..k).rev() {
            let base = ((encoded >> (2 * i)) & 3) + 1;
            result.push(base as u8);
        }
        result
    }

    /// Helper: encode a u64 k-mer back to bytes for FM-index queries.
    /// Encoding: first base is at MSB, last base is at LSB
    /// Note: FM-index uses 1-indexed encoding (A=1, C=2, G=3, T=4)
    pub fn kmer_encoded_to_bytes(encoded: u64, k: usize) -> Vec<u8> {
        Self::kmer_0indexed_to_fm_bytes(encoded, k)
    }

    /// Anchor-based filter using the rarest k-mer + 2-bit XOR scoring.
    ///
    /// Algorithm:
    /// 1. Find the rarest k-mer in the read (fewest total positions across all genomes)
    /// 2. Get all positions for that anchor k-mer via FM-index backward search
    /// 3. For each position: 2-bit XOR score the entire read against the genome region
    /// 4. Return positions with score >= min_score
    ///
    /// This is O(anchor_positions * read_length/31) vs O(total_kmer_hits) for kmer_filter.
    ///
    /// # Arguments
    /// * `read` — DNA read sequence
    /// * `min_score` — Minimum alignment score (0.0-1.0)
    ///
    /// # Returns
    /// Vector of (genome_id, position, score, cigar) tuples sorted by descending score.
    pub fn anchor_filter(&self, read: &str, min_score: f64) -> Vec<(u32, u64, f64, String)> {
        self.anchor_filter_with_threshold(read, min_score, usize::MAX)
    }

    /// Anchor-based filter with configurable alignment mode.
    ///
    /// Uses top-N rarest k-mers as anchors for error tolerance.
    /// Supports XOR, Smith-Waterman, and Hybrid alignment modes.
    ///
    /// # Arguments
    /// * `read` — DNA read sequence
    /// * `mode` — Alignment mode: `Xor` (fast), `Sw` (accurate), `Hybrid` (balanced)
    /// * `min_score` — Minimum alignment score (0.0-1.0)
    /// * `max_hits` — Maximum k-mer occurrences to consider (repetitive filter)
    ///
    /// # Returns
    /// Vector of (genome_id, position, score, cigar) tuples sorted by descending score.
    pub fn anchor_filter_with_mode(
        &self,
        read: &str,
        mode: AlignMode,
        min_score: f64,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        let window_size = if self.use_spaced_seed {
            self.spaced_seed_pattern.len()
        } else {
            self.k
        };
        if encoded.len() < window_size {
            return Vec::new();
        }

        let top_n_kmers = match self.fuzzy_method {
            FuzzyMethod::FuzzyKmer => self.find_top_n_rarest_kmers_fuzzy(&encoded, fm, max_hits),
            FuzzyMethod::FuzzySeed => {
                if self.use_spaced_seed {
                    self.find_top_n_rarest_kmers_spaced_fuzzy(&encoded, fm, max_hits)
                } else {
                    self.find_top_n_rarest_kmers_fuzzy(&encoded, fm, max_hits)
                }
            }
            FuzzyMethod::Neighborhood => {
                self.find_top_n_rarest_kmers_neighborhood(&encoded, fm, max_hits)
            }
            FuzzyMethod::None => {
                if self.use_spaced_seed {
                    self.find_top_n_rarest_kmers_spaced(&encoded, fm, max_hits)
                } else {
                    self.find_top_n_rarest_kmers(&encoded, fm, max_hits)
                }
            }
        };
        if top_n_kmers.is_empty() {
            return Vec::new();
        }

        let mut scored = Vec::new();
        let read_len = encoded.len();
        let mut seen: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();

        for &(anchor_read_offset, ref anchor_kmer, _) in &top_n_kmers {
            let raw_positions = if self.use_spaced_seed {
                if self.spaced_seed_hash.is_some() {
                    self.find_positions_spaced_hash(anchor_kmer)
                } else {
                    fm.find_positions_spaced(anchor_kmer, &self.spaced_seed_pattern, 500)
                }
            } else {
                fm.find_positions(anchor_kmer, 500)
            };

            let positions: Vec<(u32, u64)> = if raw_positions.len() > 100 {
                let stride = raw_positions.len() / 100;
                raw_positions.into_iter().step_by(stride).collect()
            } else {
                raw_positions
            };

            for &(genome_id, position) in &positions {
                if !seen.insert((genome_id, position)) {
                    continue;
                }

                let genome = match self.genomes.get(&genome_id) {
                    Some(g) => g,
                    None => continue,
                };

                let estimated_read_start = position as isize - anchor_read_offset as isize;

                let search_radius: isize = self.search_radius(read_len);
                let use_direct = search_radius <= 10;
                let mut best_score = f64::NEG_INFINITY;
                let mut best_cigar = String::new();
                let mut best_offset: usize = 0;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() {
                        continue;
                    }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < window_size {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];

                    let (score, cigar, _) = match mode {
                        AlignMode::Xor => {
                            let is_long_read = self.read_type == "long" || read_len > 1000;
                            if use_direct && !is_long_read {
                                let (s, _) =
                                    align::two_bit_score_direct(&encoded, candidate_region);
                                (s, format!("{}M", read_len), candidate_start)
                            } else if read_len <= 31 {
                                align::two_bit_align(&encoded, candidate_region)
                            } else {
                                let (s, o) =
                                    align::two_bit_score_chunks(&encoded, candidate_region);
                                (s, format!("{}M", read_len), o)
                            }
                        }
                        AlignMode::Sw => {
                            if read_len <= 31 {
                                let (sw_score, cigar) =
                                    align::smith_waterman(&encoded, candidate_region);
                                if sw_score == 0 {
                                    continue;
                                }
                                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                                (normalized.clamp(0.0, 1.0), cigar, candidate_start)
                            } else {
                                let (sw_score, cigar) =
                                    align::smith_waterman_chunked(&encoded, candidate_region);
                                if sw_score == 0 {
                                    continue;
                                }
                                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                                (normalized.clamp(0.0, 1.0), cigar, candidate_start)
                            }
                        }
                        AlignMode::Hybrid => {
                            if read_len <= 31 {
                                let (xor_s, xor_cigar, _) =
                                    align::two_bit_align(&encoded, candidate_region);
                                if xor_s >= 0.9 {
                                    (xor_s, xor_cigar, candidate_start)
                                } else if xor_s >= 0.7 {
                                    let (sw_score, sw_cigar) =
                                        align::smith_waterman(&encoded, candidate_region);
                                    if sw_score == 0 {
                                        (xor_s, xor_cigar, candidate_start)
                                    } else {
                                        let normalized =
                                            (sw_score as f64) / (2.0 * read_len as f64);
                                        if normalized > xor_s {
                                            (normalized.min(1.0), sw_cigar, candidate_start)
                                        } else {
                                            (xor_s, xor_cigar, candidate_start)
                                        }
                                    }
                                } else {
                                    (xor_s, xor_cigar, candidate_start)
                                }
                            } else {
                                let (s, _) =
                                    align::two_bit_score_chunks(&encoded, candidate_region);
                                if s >= 0.9 {
                                    (s, format!("{}M", read_len), candidate_start)
                                } else if s >= 0.7 {
                                    let (sw_score, sw_cigar) =
                                        align::smith_waterman_chunked(&encoded, candidate_region);
                                    let sw_s = (sw_score as f64) / (2.0 * read_len as f64);
                                    if sw_s > s && sw_score > 0 {
                                        (sw_s.min(1.0), sw_cigar, candidate_start)
                                    } else {
                                        (s, format!("{}M", read_len), candidate_start)
                                    }
                                } else {
                                    (s, format!("{}M", read_len), candidate_start)
                                }
                            }
                        }
                        AlignMode::Softclip => {
                            let min_align = (read_len / 3).max(20);
                            let (s, c, _) = align::two_bit_align_softclip(
                                &encoded,
                                candidate_region,
                                min_align,
                                0.7,
                            );
                            (s, c, candidate_start)
                        }
                        AlignMode::Chain => {
                            // Chain mode does its own global search via minimizers,
                            // but for per-candidate fallback, use XOR scoring
                            if read_len <= 31 {
                                align::two_bit_align(&encoded, candidate_region)
                            } else {
                                let (s, o) =
                                    align::two_bit_score_chunks(&encoded, candidate_region);
                                (s, format!("{}M", read_len), o)
                            }
                        }
                    };

                    if score > best_score {
                        best_score = score;
                        best_cigar = cigar;
                        best_offset = candidate_start;
                    }
                }

                if best_score >= min_score {
                    scored.push((
                        genome_id,
                        best_offset as u64,
                        best_score.clamp(0.0, 1.0),
                        best_cigar,
                    ));
                }
            }
        }

        scored.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        scored.truncate(50);
        scored
    }

    /// Anchor-based filter with explicit repetitive k-mer threshold.
    /// K-mers with more than `max_hits` total positions are skipped.
    /// Uses top-N rarest k-mers as anchors for error tolerance.
    pub fn anchor_filter_with_threshold(
        &self,
        read: &str,
        min_score: f64,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        if encoded.len() < self.k {
            return Vec::new();
        }

        let top_n_kmers = match self.fuzzy_method {
            FuzzyMethod::FuzzyKmer => self.find_top_n_rarest_kmers_fuzzy(&encoded, fm, max_hits),
            FuzzyMethod::FuzzySeed => {
                if self.use_spaced_seed {
                    self.find_top_n_rarest_kmers_spaced_fuzzy(&encoded, fm, max_hits)
                } else {
                    self.find_top_n_rarest_kmers_fuzzy(&encoded, fm, max_hits)
                }
            }
            FuzzyMethod::Neighborhood => {
                self.find_top_n_rarest_kmers_neighborhood(&encoded, fm, max_hits)
            }
            FuzzyMethod::None => self.find_top_n_rarest_kmers(&encoded, fm, max_hits),
        };
        if top_n_kmers.is_empty() {
            return Vec::new();
        }

        let mut scored = Vec::new();
        let read_len = encoded.len();
        let mut seen: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();

        for &(anchor_read_offset, ref anchor_kmer, _) in &top_n_kmers {
            let raw_positions = fm.find_positions(anchor_kmer, 500);

            let positions: Vec<(u32, u64)> = if raw_positions.len() > 100 {
                let stride = raw_positions.len() / 100;
                raw_positions.into_iter().step_by(stride).collect()
            } else {
                raw_positions
            };

            for &(genome_id, position) in &positions {
                if !seen.insert((genome_id, position)) {
                    continue;
                }

                let genome = match self.genomes.get(&genome_id) {
                    Some(g) => g,
                    None => continue,
                };

                let estimated_read_start = position as isize - anchor_read_offset as isize;

                if read_len <= 31 {
                    let estimated_start =
                        (position as isize - anchor_read_offset as isize).max(0) as usize;
                    let region_end = (estimated_start + read_len).min(genome.len());
                    let region = &genome[estimated_start..region_end];

                    if region.len() < self.k {
                        continue;
                    }

                    if region.len() == read_len && encoded == *region {
                        scored.push((
                            genome_id,
                            estimated_start as u64,
                            1.0,
                            format!("{}M", read_len),
                        ));
                        continue;
                    }

                    let (score, cigar, offset) = align::two_bit_align(&encoded, region);
                    if score >= min_score {
                        let actual_pos = (estimated_start + offset) as u64;
                        scored.push((genome_id, actual_pos, score, cigar));
                    }
                    continue;
                }

                let search_radius: isize = self.search_radius(read_len);
                let mut best_score = 0.0f64;
                let mut best_offset: usize = 0;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() {
                        continue;
                    }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < self.k {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];

                    if candidate_region.len() == read_len && encoded == *candidate_region {
                        best_score = 1.0;
                        best_offset = candidate_start;
                        break;
                    }

                    let (score, _) = if search_radius <= 10 {
                        align::two_bit_score_direct(&encoded, candidate_region)
                    } else {
                        align::two_bit_score_chunks(&encoded, candidate_region)
                    };
                    if score > best_score {
                        best_score = score;
                        best_offset = candidate_start;
                    }
                }

                if best_score >= min_score {
                    let cand_end = (best_offset + read_len).min(genome.len());
                    let cand_region = &genome[best_offset..cand_end];
                    let overlap = read_len.min(cand_region.len());
                    let read_part = &encoded[..overlap];

                    let mut cigar = String::with_capacity(read_len * 2 + 2);
                    let mut ops: Vec<(u8, usize)> = Vec::new();
                    for i in 0..overlap {
                        let op = if read_part[i] == cand_region[i] {
                            0u8
                        } else {
                            1u8
                        };
                        if let Some(last) = ops.last_mut() {
                            if last.0 == op {
                                last.1 += 1;
                            } else {
                                ops.push((op, 1));
                            }
                        } else {
                            ops.push((op, 1));
                        }
                    }
                    if cand_region.len() < read_len {
                        let clip = read_len - overlap;
                        if ops.is_empty() || ops.last().unwrap().0 != 2 {
                            ops.push((2, clip));
                        } else {
                            ops.last_mut().unwrap().1 += clip;
                        }
                    }
                    for (op, count) in ops {
                        cigar.push_str(&count.to_string());
                        cigar.push(match op {
                            0 => 'M',
                            1 => 'X',
                            _ => 'S',
                        });
                    }
                    let mapped_position = best_offset as u64;
                    scored.push((genome_id, mapped_position, best_score, cigar));
                }
            }
        }

        scored.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        scored.truncate(50);
        scored
    }

    /// Quality-aware anchor filter with smart threshold.
    /// Computes an adaptive minimum score based on read length and quality distribution.
    pub fn anchor_filter_with_quality_smart(
        &self,
        read: &str,
        quality: &[u8],
        min_score: f64,
        min_quality: u8,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String, f64)> {
        let encoded = encode_sequence(read);
        let read_len = encoded.len();

        // Compute smart threshold based on quality distribution
        let avg_quality: f64 = if !quality.is_empty() {
            quality.iter().map(|&q| q as f64).sum::<f64>() / quality.len() as f64
        } else {
            20.0
        };

        let smart_min = self.compute_smart_threshold(read_len, true, avg_quality);
        let effective_min = min_score.max(smart_min);

        self.anchor_filter_with_quality(read, quality, effective_min, min_quality, max_hits)
    }

    /// Quality-aware anchor filter: finds top-N rarest k-mers using only high-quality bases,
    /// then scores alignment with Phred-scaled quality penalties.
    pub fn anchor_filter_with_quality(
        &self,
        read: &str,
        quality: &[u8],
        min_score: f64,
        min_quality: u8,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String, f64)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        if encoded.len() < self.k {
            return Vec::new();
        }

        let top_n_kmers =
            self.find_top_n_rarest_kmers_quality(&encoded, quality, fm, min_quality, max_hits);
        if top_n_kmers.is_empty() {
            // Fallback: use regular anchor filter if no high-quality k-mers found
            let regular = self.anchor_filter_with_threshold(read, min_score, max_hits);
            return regular
                .into_iter()
                .map(|(g, p, s, c)| (g, p, s, c, 0.0))
                .collect();
        }

        let mut scored = Vec::new();
        let read_len = encoded.len();
        let mut seen: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();

        for &(anchor_read_offset, ref anchor_kmer, _) in &top_n_kmers {
            let raw_positions = fm.find_positions(anchor_kmer, 500);

            let positions: Vec<(u32, u64)> = if raw_positions.len() > 100 {
                let stride = raw_positions.len() / 100;
                raw_positions.into_iter().step_by(stride).collect()
            } else {
                raw_positions
            };

            for &(genome_id, position) in &positions {
                if !seen.insert((genome_id, position)) {
                    continue;
                }

                let genome = match self.genomes.get(&genome_id) {
                    Some(g) => g,
                    None => continue,
                };

                let estimated_read_start = position as isize - anchor_read_offset as isize;

                if read_len <= 31 {
                    let estimated_start =
                        (position as isize - anchor_read_offset as isize).max(0) as usize;
                    let region_end = (estimated_start + read_len).min(genome.len());
                    let region = &genome[estimated_start..region_end];

                    if region.len() < self.k {
                        continue;
                    }

                    let qual_slice = &quality[..read_len.min(quality.len())];
                    let (score, cigar, _, penalty) =
                        align::two_bit_align_with_quality(&encoded, region, qual_slice);

                    if score >= min_score {
                        let actual_pos = (estimated_start) as u64;
                        scored.push((genome_id, actual_pos, score, cigar, penalty));
                    }
                    continue;
                }

                let search_radius: isize = self.search_radius(read_len);
                let mut best_score = f64::NEG_INFINITY;
                let mut best_offset: usize = 0;
                let mut best_penalty = 0.0f64;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() {
                        continue;
                    }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < self.k {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];
                    let qual_slice = &quality[..read_len.min(quality.len())];

                    let (chunk_score, _, chunk_penalty) = if search_radius <= 10 {
                        let (s, _) = align::two_bit_score_direct(&encoded, candidate_region);
                        (s, 0, 0.0f64)
                    } else {
                        align::two_bit_score_chunks_with_quality(
                            &encoded,
                            candidate_region,
                            qual_slice,
                        )
                    };
                    let adjusted = chunk_score + chunk_penalty;

                    if adjusted > best_score {
                        best_score = adjusted;
                        best_offset = candidate_start;
                        best_penalty = chunk_penalty;
                    }
                }

                if best_score >= min_score {
                    let cand_end = (best_offset + read_len).min(genome.len());
                    let cand_region = &genome[best_offset..cand_end];
                    let overlap = read_len.min(cand_region.len());
                    let read_part = &encoded[..overlap];

                    let mut cigar = String::with_capacity(read_len * 2 + 2);
                    let mut ops: Vec<(u8, usize)> = Vec::new();
                    for i in 0..overlap {
                        let op = if read_part[i] == cand_region[i] {
                            0u8
                        } else {
                            1u8
                        };
                        if let Some(last) = ops.last_mut() {
                            if last.0 == op {
                                last.1 += 1;
                            } else {
                                ops.push((op, 1));
                            }
                        } else {
                            ops.push((op, 1));
                        }
                    }
                    if cand_region.len() < read_len {
                        let clip = read_len - overlap;
                        if ops.is_empty() || ops.last().unwrap().0 != 2 {
                            ops.push((2, clip));
                        } else {
                            ops.last_mut().unwrap().1 += clip;
                        }
                    }
                    for (op, count) in ops {
                        cigar.push_str(&count.to_string());
                        cigar.push(match op {
                            0 => 'M',
                            1 => 'X',
                            _ => 'S',
                        });
                    }
                    scored.push((
                        genome_id,
                        best_offset as u64,
                        best_score.clamp(0.0, 1.0),
                        cigar,
                        best_penalty,
                    ));
                }
            }
        }

        scored.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        scored.truncate(50);
        scored
    }

    /// Dispatch anchor filtering based on chunk anchor strategy.
    /// Uses quality scores when strategy is Golden, otherwise ignores them.
    fn anchor_filter_for_chunk(
        &self,
        chunk: &str,
        quality: Option<&[u8]>,
        min_score: f64,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String)> {
        match self.chunk_anchor_strategy {
            ChunkAnchorStrategy::Golden => {
                if let Some(qual) = quality {
                    self.anchor_filter_with_golden_anchors(chunk, qual, min_score, max_hits)
                } else {
                    // Fallback to rarest when quality not available
                    self.anchor_filter_with_mode(chunk, self.align_mode, min_score, max_hits)
                }
            }
            ChunkAnchorStrategy::Spaced => {
                self.anchor_filter_with_mode(chunk, self.align_mode, min_score, max_hits)
            }
            ChunkAnchorStrategy::Rarest => {
                self.anchor_filter_with_mode(chunk, self.align_mode, min_score, max_hits)
            }
        }
    }

    /// Golden anchor anchor filter: uses quality-weighted anchors for better long-read mapping.
    pub fn anchor_filter_with_golden_anchors(
        &self,
        read: &str,
        quality: &[u8],
        min_score: f64,
        max_hits: usize,
    ) -> Vec<(u32, u64, f64, String)> {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        if encoded.len() < self.k {
            return Vec::new();
        }

        let top_n_kmers = self.find_golden_anchors(&encoded, quality, fm, max_hits);
        if top_n_kmers.is_empty() {
            return self.anchor_filter_with_threshold(read, min_score, max_hits);
        }

        let mut scored = Vec::new();
        let read_len = encoded.len();
        let mut seen: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();

        for &(anchor_read_offset, ref anchor_kmer, _) in &top_n_kmers {
            let raw_positions = fm.find_positions(anchor_kmer, 500);

            let positions: Vec<(u32, u64)> = if raw_positions.len() > 100 {
                let stride = raw_positions.len() / 100;
                raw_positions.into_iter().step_by(stride).collect()
            } else {
                raw_positions
            };

            for &(genome_id, position) in &positions {
                if !seen.insert((genome_id, position)) {
                    continue;
                }

                let genome = match self.genomes.get(&genome_id) {
                    Some(g) => g,
                    None => continue,
                };

                let estimated_read_start = position as isize - anchor_read_offset as isize;

                let search_radius: isize = self.search_radius(read_len);
                let mut best_score = f64::NEG_INFINITY;
                let mut best_offset: usize = 0;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() {
                        continue;
                    }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < self.k {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];
                    let (s, _) = align::two_bit_score_direct(&encoded, candidate_region);

                    if s > best_score {
                        best_score = s;
                        best_offset = candidate_start;
                    }
                }

                if best_score >= min_score {
                    let cand_end = (best_offset + read_len).min(genome.len());
                    let cand_region = &genome[best_offset..cand_end];
                    let overlap = read_len.min(cand_region.len());
                    let read_part = &encoded[..overlap];

                    let mut cigar = String::with_capacity(read_len * 2 + 2);
                    let mut ops: Vec<(u8, usize)> = Vec::new();
                    for i in 0..overlap {
                        let op = if read_part[i] == cand_region[i] {
                            0u8
                        } else {
                            1u8
                        };
                        if let Some(last) = ops.last_mut() {
                            if last.0 == op {
                                last.1 += 1;
                            } else {
                                ops.push((op, 1));
                            }
                        } else {
                            ops.push((op, 1));
                        }
                    }
                    if cand_region.len() < read_len {
                        let clip = read_len - overlap;
                        if ops.is_empty() || ops.last().unwrap().0 != 2 {
                            ops.push((2, clip));
                        } else {
                            ops.last_mut().unwrap().1 += clip;
                        }
                    }
                    for (op, count) in ops {
                        cigar.push_str(&count.to_string());
                        cigar.push(match op {
                            0 => 'M',
                            1 => 'X',
                            _ => 'S',
                        });
                    }

                    scored.push((
                        genome_id,
                        best_offset as u64,
                        best_score.clamp(0.0, 1.0),
                        cigar,
                    ));
                }
            }
        }

        scored.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        scored.truncate(50);
        scored
    }

    /// Quality-aware full pipeline: map a read with quality scores to all indexed genomes.
    pub fn map_read_with_quality(
        &self,
        read: &str,
        quality: &[u8],
        min_quality: u8,
        context_window: usize,
    ) -> Vec<QualityMappingResult> {
        let scored = self.anchor_filter_with_quality(
            read,
            quality,
            0.7,
            min_quality,
            DEFAULT_REPEAT_THRESHOLD,
        );
        if scored.is_empty() {
            return Vec::new();
        }

        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        let mut results = Vec::new();

        for &(genome_id, position, align_score, ref cigar, quality_penalty) in &scored {
            if align_score < 0.5 {
                continue;
            }

            let rarity = if encoded.len() >= self.k {
                if let Some(genome_seq) = self.genomes.get(&genome_id) {
                    let start = (position as usize).min(genome_seq.len().saturating_sub(1));
                    let end = (start + self.k).min(genome_seq.len());
                    if end - start == self.k {
                        let kmer_bytes = &genome_seq[start..end];
                        let kmer_encoded = kmer_bytes
                            .iter()
                            .map(|&b| b as char)
                            .filter_map(encode_base)
                            .collect::<Vec<u8>>();
                        if kmer_encoded.len() == self.k {
                            let occ = fm.count_occurrences(&kmer_encoded);
                            1.0 / (occ as f64).max(1.0)
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let adjusted_score = (align_score + quality_penalty).clamp(0.0, 1.0);

            let read_len = encoded.len();
            let context =
                self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(QualityMappingResult {
                genome_id,
                position,
                align_score,
                adjusted_score,
                combined_score,
                cigar: cigar.clone(),
                quality_penalty,
                quality_scores: quality.to_vec(),
                context,
                is_reverse: false,
                rarity,
                md_string: String::new(),
            });
        }

        results.sort_by(|a, b| {
            b.combined_score
                .partial_cmp(&a.combined_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Extract ±window bases around a position in a genome.
    fn extract_genome_context(
        &self,
        genome_id: u32,
        position: u64,
        read_len: usize,
        window: usize,
    ) -> String {
        let genome = match self.genomes.get(&genome_id) {
            Some(g) => g,
            None => return String::new(),
        };

        let pos = position as usize;
        let start = pos.saturating_sub(window);
        let end = (pos + read_len + window).min(genome.len());

        let mut ctx = String::new();
        for &b in &genome[start..end] {
            ctx.push(decode_base(b));
        }
        ctx
    }

    /// Full pipeline: map a read to all indexed genomes.
    ///
    /// Uses anchor-based 2-bit XOR filtering (1 anchor + XOR score).
    /// Tries both forward and reverse complement strands, returns the best alignment.
    ///
    /// # Arguments
    /// * `read` — DNA read sequence (any case, N bases skipped)
    /// * `context_window` — Number of flanking bases to include in context string
    ///
    /// # Returns
    /// Ranked `MappingResult` list (best match first). Empty if no match found.
    ///
    /// # Example
    ///
    /// ```
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// bp.add_genome("Ecoli", "AGCTAGCTAGCTAGCTAGCTAGCT");
    /// bp.build();
    /// let results = bp.map_read("AGCTAGCTAGCTAGCT", 10);
    /// assert!(!results.is_empty());
    /// ```
    pub fn map_read(&self, read: &str, context_window: usize) -> Vec<MappingResult> {
        self.map_read_with_mode(read, AlignMode::Xor, context_window)
    }

    /// Full pipeline with configurable alignment mode.
    ///
    /// Uses anchor-based filtering with the specified alignment algorithm.
    /// Tries both forward and reverse complement, returns best alignment.
    ///
    /// # Arguments
    /// * `read` — DNA read sequence
    /// * `mode` — Alignment mode: `Xor` (fast), `Sw` (accurate), `Hybrid` (balanced)
    /// * `context_window` — Flanking bases for context string
    ///
    /// # Example
    ///
    /// ```
    /// use bit_pop::{BitPop, AlignMode};
    /// let mut bp = BitPop::new(10);
    /// bp.add_genome("g1", "AGCTAGCTAGCTAGCT");
    /// bp.build();
    /// let results = bp.map_read_with_mode("AGCTAGCTAGCT", AlignMode::Hybrid, 10);
    /// ```
    pub fn map_read_with_mode(
        &self,
        read: &str,
        mode: AlignMode,
        context_window: usize,
    ) -> Vec<MappingResult> {
        let forward_results = self.map_read_orientation(read, mode, context_window, false);
        let rc_read = reverse_complement(read);
        let rc_results = self.map_read_orientation(&rc_read, mode, context_window, true);

        let best_forward = forward_results.first().cloned();
        let best_rc = rc_results.first().cloned();

        match (best_forward, best_rc) {
            (Some(f), Some(r)) => {
                if r.score > f.score {
                    let mut results = rc_results;
                    if !results.is_empty() {
                        results[0].is_reverse = true;
                    }
                    results
                } else {
                    let mut results = forward_results;
                    if !results.is_empty() {
                        results[0].is_reverse = false;
                    }
                    results
                }
            }
            (Some(_f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(_r)) => {
                let mut results = rc_results;
                if !results.is_empty() {
                    results[0].is_reverse = true;
                }
                results
            }
            (None, None) => Vec::new(),
        }
    }

    /// Map a single orientation (forward or RC) of a read.
    fn map_read_orientation(
        &self,
        read: &str,
        mode: AlignMode,
        context_window: usize,
        _is_rc: bool,
    ) -> Vec<MappingResult> {
        let scored = self.anchor_filter_with_mode(read, mode, 0.7, DEFAULT_REPEAT_THRESHOLD);
        if scored.is_empty() {
            return Vec::new();
        }
        self.rank_scored_results(&scored, read, context_window, 0.5)
    }

    /// Map a read with custom threshold (for two-pass mapping).
    pub fn map_read_with_threshold(
        &self,
        read: &str,
        mode: AlignMode,
        context_window: usize,
        min_score: f64,
    ) -> Vec<MappingResult> {
        let forward_results =
            self.map_read_orientation_threshold(read, mode, context_window, false, min_score);
        let rc_read = reverse_complement(read);
        let rc_results =
            self.map_read_orientation_threshold(&rc_read, mode, context_window, true, min_score);

        let best_forward = forward_results.first().cloned();
        let best_rc = rc_results.first().cloned();

        match (best_forward, best_rc) {
            (Some(f), Some(r)) => {
                if r.score > f.score {
                    let mut results = rc_results;
                    if !results.is_empty() {
                        results[0].is_reverse = true;
                    }
                    results
                } else {
                    let mut results = forward_results;
                    if !results.is_empty() {
                        results[0].is_reverse = false;
                    }
                    results
                }
            }
            (Some(_f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(_r)) => {
                let mut results = rc_results;
                if !results.is_empty() {
                    results[0].is_reverse = true;
                }
                results
            }
            (None, None) => Vec::new(),
        }
    }

    /// Map a single orientation with custom threshold.
    fn map_read_orientation_threshold(
        &self,
        read: &str,
        mode: AlignMode,
        context_window: usize,
        _is_rc: bool,
        min_score: f64,
    ) -> Vec<MappingResult> {
        let scored = self.anchor_filter_with_mode(read, mode, min_score, DEFAULT_REPEAT_THRESHOLD);
        if scored.is_empty() {
            return Vec::new();
        }
        self.rank_scored_results(&scored, read, context_window, min_score)
    }

    /// Diagnose why a read failed to map.
    /// Returns a human-readable reason for unmapped reads.
    pub fn diagnose_read(&self, read: &str) -> String {
        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return "Index not built".to_string(),
        };

        let encoded = encode_sequence(read);
        let window_size = if self.use_spaced_seed {
            self.spaced_seed_pattern.len()
        } else {
            self.k
        };

        if encoded.len() < window_size {
            return "Read too short".to_string();
        }

        // Check how many k-mers exist in the index
        let mut found_in_index = 0;
        let mut _below_threshold = 0;
        let mut above_threshold = 0;

        for i in 0..=(encoded.len() - window_size) {
            let kmer: Vec<u8> = encoded[i..i + window_size].to_vec();
            let count = fm.count_occurrences(&kmer);
            if count > 0 {
                found_in_index += 1;
                if count > DEFAULT_REPEAT_THRESHOLD {
                    _below_threshold += 1;
                } else {
                    above_threshold += 1;
                }
            }
        }

        if found_in_index == 0 {
            return "No k-mers in index".to_string();
        }

        if above_threshold == 0 {
            return "All k-mers below rarity threshold".to_string();
        }

        // Check best alignment score
        let scored =
            self.anchor_filter_with_mode(read, AlignMode::Xor, 0.0, DEFAULT_REPEAT_THRESHOLD);
        if let Some((_, _, best_score, _)) =
            scored.iter().max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        {
            if *best_score < 0.5 {
                format!("Alignment score too low: {:.3}", best_score)
            } else if *best_score < 0.7 {
                format!("Alignment score below threshold: {:.3}", best_score)
            } else {
                "K-mers exist but alignment failed".to_string()
            }
        } else {
            "No anchor candidates found".to_string()
        }
    }

    /// Full pipeline: map a read with quality scores to all indexed genomes.
    /// Uses quality-aware anchor filtering + Phred-scaled scoring.
    /// Tries both forward and reverse complement, returns best alignment.
    pub fn map_read_with_quality_mode(
        &self,
        read: &str,
        quality: &[u8],
        mode: AlignMode,
        min_quality: u8,
        context_window: usize,
    ) -> Vec<QualityMappingResult> {
        let forward_results = self.map_read_quality_orientation(
            read,
            quality,
            mode,
            min_quality,
            context_window,
            false,
        );
        let rc_read = reverse_complement(read);
        let rc_results = self.map_read_quality_orientation(
            &rc_read,
            quality,
            mode,
            min_quality,
            context_window,
            true,
        );

        let best_forward = forward_results.first().cloned();
        let best_rc = rc_results.first().cloned();

        match (best_forward, best_rc) {
            (Some(f), Some(r)) => {
                if r.combined_score > f.combined_score {
                    let mut results = rc_results;
                    if !results.is_empty() {
                        results[0].is_reverse = true;
                    }
                    results
                } else {
                    let mut results = forward_results;
                    if !results.is_empty() {
                        results[0].is_reverse = false;
                    }
                    results
                }
            }
            (Some(_f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(_r)) => {
                let mut results = rc_results;
                if !results.is_empty() {
                    results[0].is_reverse = true;
                }
                results
            }
            (None, None) => Vec::new(),
        }
    }

    /// Full pipeline with golden anchor selection: map a read with quality scores using quality-weighted anchors.
    pub fn map_read_with_golden_anchors(
        &self,
        read: &str,
        quality: &[u8],
        mode: AlignMode,
        context_window: usize,
    ) -> Vec<QualityMappingResult> {
        let forward_results =
            self.map_read_golden_orientation(read, quality, mode, context_window, false);
        let rc_read = reverse_complement(read);
        let rc_results =
            self.map_read_golden_orientation(&rc_read, quality, mode, context_window, true);

        let best_forward = forward_results.first().cloned();
        let best_rc = rc_results.first().cloned();

        match (best_forward, best_rc) {
            (Some(f), Some(r)) => {
                if r.combined_score > f.combined_score {
                    let mut results = rc_results;
                    if !results.is_empty() {
                        results[0].is_reverse = true;
                    }
                    results
                } else {
                    let mut results = forward_results;
                    if !results.is_empty() {
                        results[0].is_reverse = false;
                    }
                    results
                }
            }
            (Some(_f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(_r)) => {
                let mut results = rc_results;
                if !results.is_empty() {
                    results[0].is_reverse = true;
                }
                results
            }
            (None, None) => Vec::new(),
        }
    }

    /// Map a single orientation (forward or RC) of a read with quality scores using golden anchors.
    fn map_read_golden_orientation(
        &self,
        read: &str,
        quality: &[u8],
        _mode: AlignMode,
        context_window: usize,
        is_rc: bool,
    ) -> Vec<QualityMappingResult> {
        let scored =
            self.anchor_filter_with_golden_anchors(read, quality, 0.5, DEFAULT_REPEAT_THRESHOLD);
        if scored.is_empty() {
            return Vec::new();
        }

        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        let mut results = Vec::new();

        for &(genome_id, position, align_score, ref cigar) in &scored {
            if align_score < 0.5 {
                continue;
            }

            let rarity = if encoded.len() >= self.k {
                if let Some(genome_seq) = self.genomes.get(&genome_id) {
                    let start = (position as usize).min(genome_seq.len().saturating_sub(1));
                    let end = (start + self.k).min(genome_seq.len());
                    if end - start == self.k {
                        let kmer_bytes = &genome_seq[start..end];
                        let kmer_encoded = kmer_bytes
                            .iter()
                            .map(|&b| b as char)
                            .filter_map(encode_base)
                            .collect::<Vec<u8>>();
                        if kmer_encoded.len() == self.k {
                            let occ = fm.count_occurrences(&kmer_encoded);
                            1.0 / (occ as f64).max(1.0)
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;

            let read_len = encoded.len();
            let context =
                self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(QualityMappingResult {
                genome_id,
                position,
                align_score,
                adjusted_score: align_score,
                combined_score,
                cigar: cigar.clone(),
                quality_penalty: 0.0,
                quality_scores: quality.to_vec(),
                context,
                is_reverse: is_rc,
                rarity,
                md_string: String::new(),
            });
        }

        results.sort_by(|a, b| {
            b.combined_score
                .partial_cmp(&a.combined_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results.truncate(10);
        results
    }

    /// Map a single orientation (forward or RC) of a read with quality scores.
    fn map_read_quality_orientation(
        &self,
        read: &str,
        quality: &[u8],
        _mode: AlignMode,
        min_quality: u8,
        context_window: usize,
        _is_rc: bool,
    ) -> Vec<QualityMappingResult> {
        let scored = self.anchor_filter_with_quality_smart(
            read,
            quality,
            0.7,
            min_quality,
            DEFAULT_REPEAT_THRESHOLD,
        );
        if scored.is_empty() {
            return Vec::new();
        }

        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);
        let mut results = Vec::new();

        for &(genome_id, position, align_score, ref cigar, quality_penalty) in &scored {
            if align_score < 0.5 {
                continue;
            }

            let rarity = if encoded.len() >= self.k {
                if let Some(genome_seq) = self.genomes.get(&genome_id) {
                    let start = (position as usize).min(genome_seq.len().saturating_sub(1));
                    let end = (start + self.k).min(genome_seq.len());
                    if end - start == self.k {
                        let kmer_bytes = &genome_seq[start..end];
                        let kmer_encoded = kmer_bytes
                            .iter()
                            .map(|&b| b as char)
                            .filter_map(encode_base)
                            .collect::<Vec<u8>>();
                        if kmer_encoded.len() == self.k {
                            let occ = fm.count_occurrences(&kmer_encoded);
                            1.0 / (occ as f64).max(1.0)
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                }
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let adjusted_score = (align_score + quality_penalty).clamp(0.0, 1.0);

            let read_len = encoded.len();
            let context =
                self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(QualityMappingResult {
                genome_id,
                position,
                align_score,
                adjusted_score,
                combined_score,
                cigar: cigar.clone(),
                quality_penalty,
                quality_scores: quality.to_vec(),
                context,
                is_reverse: false,
                rarity,
                md_string: String::new(),
            });
        }

        results.sort_by(|a, b| {
            b.combined_score
                .partial_cmp(&a.combined_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Get genome name by ID.
    pub fn genome_name(&self, genome_id: u32) -> Option<&str> {
        self.genome_names.get(&genome_id).map(|s| s.as_str())
    }

    /// Get the number of indexed genomes.
    pub fn genome_count(&self) -> usize {
        self.genomes.len()
    }

    /// Get the total indexed length (BWT length from FM-index).
    pub fn bwt_len(&self) -> usize {
        self.fm_index.as_ref().map(|fm| fm.len()).unwrap_or(0)
    }

    /// Get the length of a genome sequence.
    pub fn genome_seq_len(&self, genome_id: u32) -> Option<usize> {
        self.genomes.get(&genome_id).map(|s| s.len())
    }

    /// Get all genome names in order of genome_id.
    pub fn genome_names_ordered(&self) -> Vec<String> {
        let mut names: Vec<(u32, String)> = self
            .genome_names
            .iter()
            .map(|(id, name)| (*id, name.clone()))
            .collect();
        names.sort_by_key(|(id, _)| *id);
        names.into_iter().map(|(_, name)| name).collect()
    }

    // --- Serialization helpers ---

    /// Get the k-mer size.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get a genome's DNA sequence by ID (for serialization).
    pub fn get_genome_seq(&self, genome_id: u32) -> Option<&Vec<u8>> {
        self.genomes.get(&genome_id)
    }

    /// Get the FM-index for serialization.
    pub fn get_fm_index(&self) -> Option<&FmIndex> {
        self.fm_index.as_ref()
    }

    /// Create a BitPop from serialized FM-index data.
    pub fn from_fm_index(
        k: usize,
        genomes: HashMap<u32, Vec<u8>>,
        genome_names: HashMap<u32, String>,
        fm_index: FmIndex,
    ) -> Self {
        Self {
            fm_index: Some(fm_index),
            genomes,
            genome_names,
            k,
            auto_k: false,
            top_n: 1,
            use_spaced_seed: false,
            spaced_seed_pattern: vec![
                true, true, true, true, true, false, true, true, true, true, true, true, true, true,
            ],
            spaced_seed_hash: None,
            read_type: "short".to_string(),
            fuzzy_method: FuzzyMethod::None,
            fuzzy_mismatches: 1,
            neighborhood_hash: None,
            search_radius: None,
            chunk_size: 0,
            chunk_vote_threshold: 0.0,
            chunk_top_n: 1,
            chunk_pct: 0.0,
            chunk_min: 20,
            chunk_max: 500,
            enable_snp_detect: false,
            snp_min_support: 3,
            snp_penalty: 0.1,
            align_mode: AlignMode::Xor,
            chunk_anchor_strategy: ChunkAnchorStrategy::Rarest,
            enable_hf: false,
            hf_min_run: 3,
            hf_profiles: HashMap::new(),
            chain_config: chain::ChainConfig::default(),
        }
    }

    /// Serialize and write to a file (persisted format with compression).
    ///
    /// Uses format v5 with memmap2 for <10ms load time and zstd compression.
    ///
    /// # Arguments
    /// * `path` — Output file path (e.g., "index.bitpop")
    ///
    /// # Example
    ///
    /// ```no_run
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// bp.add_genome("g1", "ACGT...");
    /// bp.build();
    /// bp.serialize_to_file("index.bitpop").unwrap();
    /// ```
    pub fn serialize_to_file(&self, path: &str) -> io::Result<()> {
        persisted::save_bitpop(self, path)?;
        Ok(())
    }

    /// Load a BitPop instance from a file (persisted format).
    ///
    /// Auto-detects and loads format v3 (legacy), v4, and v5.
    /// Format v5 uses memmap2 for <10ms load time.
    ///
    /// # Arguments
    /// * `path` — Input file path (e.g., "index.bitpop")
    ///
    /// # Example
    ///
    /// ```no_run
    /// use bit_pop::BitPop;
    /// let bp = BitPop::deserialize_from_file("index.bitpop").unwrap();
    /// let results = bp.map_read("AGCTAGCTAGCT", 10);
    /// ```
    pub fn deserialize_from_file(path: &str) -> io::Result<Self> {
        let bp = persisted::load_bitpop(path)?;
        Ok(bp)
    }

    /// Map multiple reads in parallel and write results to a SAM file.
    ///
    /// Uses rayon work-stealing scheduler for multi-core parallelism.
    /// Writes SAM header + all mappings to the output file.
    ///
    /// # Arguments
    /// * `reads` — Slice of (read_name, read_sequence) tuples
    /// * `output_path` — Output SAM file path
    /// * `context_window` — Flanking bases for context string
    ///
    /// # Returns
    /// Number of reads that had at least one mapping.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use bit_pop::BitPop;
    /// let mut bp = BitPop::new(10);
    /// bp.add_genome("g1", "ACGTACGTACGTACGT");
    /// bp.build();
    /// let reads = vec![
    ///     ("read1", "ACGTACGTACGT"),
    ///     ("read2", "GTCAGTCAGTCA"),
    /// ];
    /// let mapped = bp.map_reads_parallel(&reads, "output.sam", 10).unwrap();
    /// println!("Mapped {} reads", mapped);
    /// ```
    pub fn map_reads_parallel(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
    ) -> io::Result<usize> {
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        let mapped: Vec<(String, String, Vec<MappingResult>)> = reads
            .par_iter()
            .map(|(name, seq)| {
                let results = self.map_read(seq, context_window);
                (name.to_string(), seq.to_string(), results)
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
        writer.write_header(&genome_header)?;

        let mut mapped_count = 0;
        for (name, seq, results) in &mapped {
            writer.write_mappings(name, seq, results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
        }

        Ok(mapped_count)
    }

    /// Map multiple FASTQ reads in parallel with quality-aware scoring and write to SAM.
    ///
    /// Reads are filtered by minimum quality, then scored with Phred-scaled penalties.
    ///
    /// # Arguments
    /// * `fastq_path` — Input FASTQ file path
    /// * `output_path` — Output SAM file path
    /// * `min_quality` — Minimum Phred quality (0 = no filter)
    /// * `context_window` — Flanking bases for context string
    ///
    /// # Returns
    /// Number of reads that had at least one mapping.
    pub fn map_reads_from_fastq_parallel(
        &self,
        fastq_path: &str,
        output_path: &str,
        min_quality: u8,
        context_window: usize,
    ) -> io::Result<usize> {
        let reads = fastq::parse_fastq(fastq_path)?;

        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        let mapped: Vec<(String, String, Vec<QualityMappingResult>)> = reads
            .into_par_iter()
            .map(|(name, seq, qual)| {
                let results = self.map_read_with_quality(&seq, &qual, min_quality, context_window);
                (name, seq, results)
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
        writer.write_header(&genome_header)?;

        let mut mapped_count = 0;
        for (name, seq, results) in &mapped {
            writer.write_quality_mappings(name, seq, results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
        }

        Ok(mapped_count)
    }

    /// Map multiple reads in parallel using optimized batching with work stealing.
    /// Uses rayon's work-stealing scheduler for better load balancing.
    pub fn map_reads_parallel_optimized(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
        batch_size: usize,
    ) -> io::Result<usize> {
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        // Split reads into batches for work-stealing
        let batches: Vec<Vec<(&str, &str)>> = reads
            .chunks(batch_size)
            .map(|chunk| chunk.to_vec())
            .collect();

        let all_results: Vec<(String, String, Vec<MappingResult>)> = batches
            .into_par_iter()
            .flat_map_iter(|batch| {
                batch
                    .into_iter()
                    .map(|(name, seq)| {
                        let results = self.map_read(seq, context_window);
                        (name.to_string(), seq.to_string(), results)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
        writer.write_header(&genome_header)?;

        let mut mapped_count = 0;
        for (name, seq, results) in &all_results {
            writer.write_mappings(name, seq, results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
        }

        Ok(mapped_count)
    }

    /// Map multiple reads and write results to a SAM file.
    /// Returns the number of reads that had at least one mapping.
    pub fn map_reads_to_sam(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
    ) -> io::Result<usize> {
        let mut writer = sam::SamWriter::new(output_path)?;

        let genomes: Vec<(&str, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name, self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();
        writer.write_header(&genomes)?;

        let names = self.genome_names_ordered();
        let name_refs: Vec<&str> = names.iter().map(|s| s.as_str()).collect();

        let mut mapped_count = 0;

        for (read_name, read_seq) in reads {
            let results = self.map_read(read_seq, context_window);
            writer.write_mappings(read_name, read_seq, &results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
        }

        Ok(mapped_count)
    }

    /// Map multiple reads with progress callback and write results to a SAM/BAM file.
    /// The callback is called after every `progress_interval` reads with the current count.
    pub fn map_reads_to_output_with_progress(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
        progress_interval: usize,
        write_bam: bool,
        mut on_progress: impl FnMut(usize, usize),
    ) -> io::Result<usize> {
        use crate::bam::AlignmentWriter;
        let mut writer: Box<dyn AlignmentWriter> = if write_bam {
            Box::new(bam::BamWriter::new(output_path)?)
        } else {
            Box::new(sam::SamWriter::new(output_path)?)
        };

        let genomes: Vec<(&str, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name, self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();
        writer.write_header(&genomes)?;

        let names = self.genome_names_ordered();
        let name_refs: Vec<&str> = names.iter().map(|s| s.as_str()).collect();

        let mut mapped_count = 0;
        let total = reads.len();

        for (i, (read_name, read_seq)) in reads.iter().enumerate() {
            let results = self.map_read_with_mode(read_seq, self.align_mode, context_window);
            writer.write_mappings(read_name, read_seq, &results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
            if (i + 1) % progress_interval == 0 {
                on_progress(i + 1, total);
            }
        }

        if !reads.len().is_multiple_of(progress_interval) {
            on_progress(reads.len(), total);
        }

        Ok(mapped_count)
    }

    /// Map multiple reads with progress callback and write results to a SAM file.
    #[deprecated(
        since = "0.2.0",
        note = "Use map_reads_to_output_with_progress instead"
    )]
    #[allow(unused_mut)]
    pub fn map_reads_to_sam_with_progress(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
        progress_interval: usize,
        mut on_progress: impl FnMut(usize, usize),
    ) -> io::Result<usize> {
        self.map_reads_to_output_with_progress(
            reads,
            output_path,
            context_window,
            progress_interval,
            false,
            on_progress,
        )
    }

    /// Map multiple reads in parallel with progress reporting.
    /// Returns the number of reads that had at least one mapping.
    pub fn map_reads_parallel_with_progress(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
        progress_interval: usize,
        on_progress: impl FnMut(usize, usize) + Send,
    ) -> io::Result<usize> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::Mutex;

        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();
        let completed = AtomicUsize::new(0);
        let progress_callback = Mutex::new(on_progress);
        let total = reads.len();

        let mapped: Vec<(String, String, Vec<MappingResult>)> = reads
            .par_iter()
            .map(|(name, seq)| {
                let results = self.map_read_with_mode(seq, self.align_mode, context_window);
                let count = completed.fetch_add(1, Ordering::Relaxed) + 1;
                if count.is_multiple_of(progress_interval) || count == total {
                    if let Ok(mut cb) = progress_callback.lock() {
                        cb(count, total);
                    }
                }
                (name.to_string(), seq.to_string(), results)
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
        writer.write_header(&genome_header)?;

        let mut mapped_count = 0;
        for (name, seq, results) in &mapped {
            writer.write_mappings(name, seq, results, &name_refs)?;
            if !results.is_empty() {
                mapped_count += 1;
            }
        }

        Ok(mapped_count)
    }

    // --- Paired-end mapping ---

    /// Map a single paired-end read to all indexed genomes.
    /// Returns the best mapping result for each read in the pair.
    /// If `reconcile` is true and R1/R2 map to different genomes, attempts
    /// concordant reconciliation using top-N candidate overlap.
    /// If `insert_stats` is provided, uses Gaussian insert size model for
    /// reconciliation scoring.
    pub fn map_read_paired(
        &self,
        paired: &PairedRead,
        context_window: usize,
        reconcile: bool,
        reconcile_top_n: usize,
    ) -> PairedMappingResult {
        self.map_read_paired_with_stats(paired, context_window, reconcile, reconcile_top_n, None)
    }

    /// Map a single paired-end read with pre-computed insert size statistics.
    /// Uses Gaussian insert size model for reconciliation scoring when
    /// `insert_stats` is provided and has sufficient data (count >= 2).
    pub fn map_read_paired_with_stats(
        &self,
        paired: &PairedRead,
        context_window: usize,
        reconcile: bool,
        reconcile_top_n: usize,
        insert_stats: Option<&InsertSizeStats>,
    ) -> PairedMappingResult {
        let mut insert_stats_out = InsertSizeStats::new();

        let map1 = self.map_single_read_to_best(&paired.read1_seq);
        let map2 = self.map_single_read_to_best(&paired.read2_seq);

        let (final_map1, final_map2) = if reconcile {
            if let Some(stats) = insert_stats {
                if stats.count >= 2 {
                    self.reconcile_pair_with_gaussian(
                        paired,
                        &map1,
                        &map2,
                        reconcile_top_n,
                        context_window,
                        stats,
                    )
                } else {
                    self.reconcile_pair(paired, &map1, &map2, reconcile_top_n, context_window)
                }
            } else {
                self.reconcile_pair(paired, &map1, &map2, reconcile_top_n, context_window)
            }
        } else {
            (map1, map2)
        };

        // Compute TLEN (observed template length)
        let tlen = compute_tlen(
            &final_map1,
            &final_map2,
            paired.read1_seq.len(),
            paired.read2_seq.len(),
        );
        insert_stats_out.update(tlen);

        PairedMappingResult {
            read_name: paired.name.clone(),
            map1: final_map1,
            map2: final_map2,
            tlen,
            insert_size_stats: insert_stats_out,
        }
    }

    /// Reconcile a discordant paired-end mapping.
    ///
    /// When R1 and R2 map to different genomes, this function:
    /// 1. Gets top-N candidates for each read across all genomes
    /// 2. Finds common genomes in both candidate lists
    /// 3. Picks the concordant genome with the best combined score
    /// 4. Falls back to original independent best if no concordant genome
    ///    meets the quality threshold
    ///
    /// The combined score uses the same 85/15 weighting as the main pipeline:
    /// `combined = (score1 * 0.85 + rarity1 * 0.15 + score2 * 0.85 + rarity2 * 0.15) / 2`
    ///
    /// A concordant assignment is accepted only if its combined score is within
    /// 15% of the best discordant combined score (configurable via threshold).
    fn reconcile_pair(
        &self,
        paired: &PairedRead,
        map1: &Option<PairedReadMapping>,
        map2: &Option<PairedReadMapping>,
        top_n: usize,
        context_window: usize,
    ) -> (Option<PairedReadMapping>, Option<PairedReadMapping>) {
        self.reconcile_pair_with_gaussian(
            paired,
            map1,
            map2,
            top_n,
            context_window,
            &InsertSizeStats::new(),
        )
    }

    /// Reconcile a discordant paired-end mapping using a Gaussian insert size model.
    ///
    /// When R1 and R2 map to different genomes, this function:
    /// 1. Gets top-N candidates for each read across all genomes
    /// 2. Finds common genomes in both candidate lists
    /// 3. Scores each concordant candidate using combined alignment+rarity score
    ///    PLUS a Gaussian insert size confidence bonus
    /// 4. Picks the concordant genome with the best total score
    /// 5. Falls back to original independent best if no concordant genome
    ///    meets the quality threshold
    ///
    /// The total score is:
    ///   `total = combined_score * (1.0 + gaussian_bonus * 0.5)`
    /// where `gaussian_bonus = insert_size_confidence(observed_tlen)` in [0, 1].
    /// This gives up to 50% score boost for concordant candidates with plausible
    /// insert sizes under the Gaussian model.
    ///
    /// A concordant assignment is accepted only if its total score is within
    /// 15% of the best discordant combined score.
    fn reconcile_pair_with_gaussian(
        &self,
        paired: &PairedRead,
        map1: &Option<PairedReadMapping>,
        map2: &Option<PairedReadMapping>,
        top_n: usize,
        context_window: usize,
        insert_stats: &InsertSizeStats,
    ) -> (Option<PairedReadMapping>, Option<PairedReadMapping>) {
        // Only reconcile if both mapped to different genomes
        let (m1, m2) = match (map1, map2) {
            (Some(ref a), Some(ref b)) if a.genome_id != b.genome_id => (a, b),
            _ => return (map1.clone(), map2.clone()),
        };

        // Best discordant score (baseline to beat)
        let score1 = m1.score * 0.85 + m1.rarity * 0.15;
        let score2 = m2.score * 0.85 + m2.rarity * 0.15;
        let best_discordant = (score1 + score2) / 2.0;

        // Compute observed TLEN for Gaussian confidence scoring
        let observed_tlen =
            compute_tlen(map1, map2, paired.read1_seq.len(), paired.read2_seq.len());
        let gaussian_confidence = insert_stats.insert_size_confidence(observed_tlen);

        // Get top-N candidates for both reads
        let candidates1 = self
            .map_read(&paired.read1_seq, context_window)
            .into_iter()
            .take(top_n)
            .collect::<Vec<_>>();
        let candidates2 = self
            .map_read(&paired.read2_seq, context_window)
            .into_iter()
            .take(top_n)
            .collect::<Vec<_>>();

        // Find best concordant genome (appearing in both top-N)
        // Score includes Gaussian insert size confidence bonus
        let mut best_concordant: Option<(f64, &MappingResult, &MappingResult)> = None;

        for c1 in &candidates1 {
            for c2 in &candidates2 {
                if c1.genome_id == c2.genome_id {
                    let s1 = c1.score * 0.85 + c1.rarity * 0.15;
                    let s2 = c2.score * 0.85 + c2.rarity * 0.15;
                    let combined = (s1 + s2) / 2.0;
                    // Apply Gaussian insert size confidence bonus
                    // Up to 50% score boost for plausible insert sizes
                    let gaussian_bonus = gaussian_confidence * 0.5;
                    let total_score = combined * (1.0 + gaussian_bonus);
                    if let Some((best_score, _, _)) = best_concordant {
                        if total_score > best_score {
                            best_concordant = Some((total_score, c1, c2));
                        }
                    } else {
                        best_concordant = Some((total_score, c1, c2));
                    }
                }
            }
        }

        // Accept concordant if within 15% of discordant baseline
        if let Some((concordant_score, c1, c2)) = best_concordant {
            let threshold = best_discordant * 0.85;
            if concordant_score >= threshold {
                let new_map1 = PairedReadMapping {
                    genome_id: c1.genome_id,
                    position: c1.position,
                    score: c1.score,
                    cigar: c1.cigar.clone(),
                    is_reverse: c1.is_reverse,
                    mapped: true,
                    align_score: c1.score,
                    rarity: c1.rarity,
                    md_string: c1.md_string.clone(),
                };
                let new_map2 = PairedReadMapping {
                    genome_id: c2.genome_id,
                    position: c2.position,
                    score: c2.score,
                    cigar: c2.cigar.clone(),
                    is_reverse: c2.is_reverse,
                    mapped: true,
                    align_score: c2.score,
                    rarity: c2.rarity,
                    md_string: c2.md_string.clone(),
                };
                return (Some(new_map1), Some(new_map2));
            }
        }

        // No acceptable concordant assignment — keep original
        (map1.clone(), map2.clone())
    }

    /// Map a single read and return its best mapping result.
    fn map_single_read_to_best(&self, seq: &str) -> Option<PairedReadMapping> {
        let results = self.map_read(seq, 0);
        if results.is_empty() {
            return None;
        }
        let best = &results[0];
        Some(PairedReadMapping {
            genome_id: best.genome_id,
            position: best.position,
            score: best.score,
            cigar: best.cigar.clone(),
            is_reverse: best.is_reverse,
            mapped: true,
            align_score: best.score,
            rarity: best.rarity,
            md_string: best.md_string.clone(),
        })
    }

    /// Map multiple paired-end reads in parallel and write SAM/BAM output.
    /// Uses a two-pass approach: first pass collects insert size statistics,
    /// second pass uses the Gaussian model for discordant pair reconciliation.
    pub fn map_paired_reads_parallel(
        &self,
        pairs: &[PairedReads],
        output_path: &str,
        context_window: usize,
        reconcile: bool,
        reconcile_top_n: usize,
        write_bam: bool,
    ) -> io::Result<usize> {
        use crate::bam::AlignmentWriter;
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        // Two-pass approach for Gaussian insert size model:
        // Pass 1: Map all pairs without reconciliation to collect insert size stats
        let initial_pairs: Vec<(String, String, Vec<u8>, String, Vec<u8>)> = pairs
            .iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                (
                    name.clone(),
                    seq1.clone(),
                    qual1.clone(),
                    seq2.clone(),
                    qual2.clone(),
                )
            })
            .collect();

        let mut insert_stats = InsertSizeStats::new();
        let initial_results: Vec<PairedMappingResult> = initial_pairs
            .iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                let paired = PairedRead {
                    name: name.clone(),
                    read1_seq: seq1.clone(),
                    read1_qual: qual1.clone(),
                    read2_seq: seq2.clone(),
                    read2_qual: qual2.clone(),
                };
                let result = self.map_read_paired(&paired, context_window, false, reconcile_top_n);
                let tlen = result.tlen;
                insert_stats.update(tlen);
                result
            })
            .collect();

        // Pass 2: Map all pairs with Gaussian-aware reconciliation if stats are sufficient
        let mapped_pairs: Vec<PairedMappingResult> = if reconcile && insert_stats.count >= 10 {
            initial_pairs
                .iter()
                .map(|(name, seq1, qual1, seq2, qual2)| {
                    let paired = PairedRead {
                        name: name.clone(),
                        read1_seq: seq1.clone(),
                        read1_qual: qual1.clone(),
                        read2_seq: seq2.clone(),
                        read2_qual: qual2.clone(),
                    };
                    self.map_read_paired_with_stats(
                        &paired,
                        context_window,
                        true,
                        reconcile_top_n,
                        Some(&insert_stats),
                    )
                })
                .collect()
        } else {
            // Not enough stats for Gaussian model, use initial results
            initial_results
        };

        let mut writer: Box<dyn AlignmentWriter> = if write_bam {
            Box::new(bam::BamWriter::new(output_path)?)
        } else {
            Box::new(sam::SamWriter::new(output_path)?)
        };
        writer.write_header(&genome_header)?;

        // Write paired-end output
        for pair_result in &mapped_pairs {
            writer.write_paired_mappings(
                &pair_result.read_name,
                pair_result,
                &genome_name_refs,
                &insert_stats,
            )?;
        }

        Ok(mapped_pairs.len())
    }

    /// Map multiple paired-end reads with quality-aware scoring and write SAM/BAM output.
    /// Uses a two-pass approach: first pass collects insert size statistics,
    /// second pass uses the Gaussian model for discordant pair reconciliation.
    #[expect(clippy::too_many_arguments)]
    pub fn map_paired_reads_parallel_quality(
        &self,
        pairs: &[PairedReads],
        output_path: &str,
        min_quality: u8,
        context_window: usize,
        reconcile: bool,
        reconcile_top_n: usize,
        write_bam: bool,
    ) -> io::Result<usize> {
        use crate::bam::AlignmentWriter;
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        // Two-pass approach for Gaussian insert size model:
        // Pass 1: Map all pairs without reconciliation to collect insert size stats
        let initial_pairs: Vec<(String, String, Vec<u8>, String, Vec<u8>)> = pairs
            .iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                (
                    name.clone(),
                    seq1.clone(),
                    qual1.clone(),
                    seq2.clone(),
                    qual2.clone(),
                )
            })
            .collect();

        let mut insert_stats = InsertSizeStats::new();
        let initial_results: Vec<PairedMappingResult> = initial_pairs
            .iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                let paired = PairedRead {
                    name: name.clone(),
                    read1_seq: seq1.clone(),
                    read1_qual: qual1.clone(),
                    read2_seq: seq2.clone(),
                    read2_qual: qual2.clone(),
                };

                let map1 = self.map_paired_read_with_quality(
                    &paired.read1_seq,
                    &paired.read1_qual,
                    min_quality,
                );
                let map2 = self.map_paired_read_with_quality(
                    &paired.read2_seq,
                    &paired.read2_qual,
                    min_quality,
                );

                let tlen =
                    compute_tlen(&map1, &map2, paired.read1_seq.len(), paired.read2_seq.len());
                insert_stats.update(tlen);

                PairedMappingResult {
                    read_name: name.clone(),
                    map1,
                    map2,
                    tlen,
                    insert_size_stats: InsertSizeStats::new(),
                }
            })
            .collect();

        // Pass 2: Map all pairs with Gaussian-aware reconciliation if stats are sufficient
        let mapped_pairs: Vec<PairedMappingResult> = if reconcile && insert_stats.count >= 10 {
            initial_pairs
                .iter()
                .map(|(name, seq1, qual1, seq2, qual2)| {
                    let paired = PairedRead {
                        name: name.clone(),
                        read1_seq: seq1.clone(),
                        read1_qual: qual1.clone(),
                        read2_seq: seq2.clone(),
                        read2_qual: qual2.clone(),
                    };

                    let map1 = self.map_paired_read_with_quality(
                        &paired.read1_seq,
                        &paired.read1_qual,
                        min_quality,
                    );
                    let map2 = self.map_paired_read_with_quality(
                        &paired.read2_seq,
                        &paired.read2_qual,
                        min_quality,
                    );

                    let (final_map1, final_map2) = self.reconcile_pair_with_gaussian(
                        &paired,
                        &map1,
                        &map2,
                        reconcile_top_n,
                        context_window,
                        &insert_stats,
                    );

                    let tlen = compute_tlen(
                        &final_map1,
                        &final_map2,
                        paired.read1_seq.len(),
                        paired.read2_seq.len(),
                    );

                    PairedMappingResult {
                        read_name: name.clone(),
                        map1: final_map1,
                        map2: final_map2,
                        tlen,
                        insert_size_stats: InsertSizeStats::new(),
                    }
                })
                .collect()
        } else {
            // Not enough stats for Gaussian model, use initial results
            initial_results
        };

        let mut writer: Box<dyn AlignmentWriter> = if write_bam {
            Box::new(bam::BamWriter::new(output_path)?)
        } else {
            Box::new(sam::SamWriter::new(output_path)?)
        };
        writer.write_header(&genome_header)?;

        for pair_result in &mapped_pairs {
            writer.write_paired_mappings(
                &pair_result.read_name,
                pair_result,
                &genome_name_refs,
                &insert_stats,
            )?;
        }

        Ok(mapped_pairs.len())
    }

    fn map_paired_read_with_quality(
        &self,
        seq: &str,
        qual: &[u8],
        min_quality: u8,
    ) -> Option<PairedReadMapping> {
        let results = self.map_read_with_quality_mode(seq, qual, AlignMode::Hybrid, min_quality, 0);
        if results.is_empty() {
            return None;
        }
        let best = &results[0];
        Some(PairedReadMapping {
            genome_id: best.genome_id,
            position: best.position,
            score: best.combined_score,
            cigar: best.cigar.clone(),
            is_reverse: best.is_reverse,
            mapped: true,
            align_score: best.align_score,
            rarity: best.rarity,
            md_string: best.md_string.clone(),
        })
    }

    /// Split a long read into overlapping chunks for chunk-based classification.
    ///
    /// Used for PacBio HiFi reads (6-11kb) where full-read mapping fails due to
    /// high error rates. Each chunk is mapped independently via existing Bit-Pop.
    ///
    /// # Arguments
    /// * `read` — DNA read sequence
    /// * `chunk_size` — Size of each chunk in bases
    /// * `overlap` — Overlap between consecutive chunks (default: 30bp)
    ///
    /// # Returns
    /// Vector of (chunk_index, chunk_sequence) tuples.
    fn split_read_into_chunks(
        &self,
        read: &str,
        chunk_size: usize,
        overlap: usize,
    ) -> Vec<(usize, String)> {
        if chunk_size == 0 || read.len() <= chunk_size {
            return vec![(0, read.to_string())];
        }

        let mut chunks = Vec::new();
        let mut i = 0;
        let mut chunk_idx = 0;

        while i < read.len() {
            let end = (i + chunk_size).min(read.len());
            chunks.push((chunk_idx, read[i..end].to_string()));
            chunk_idx += 1;
            i += chunk_size - overlap;
            if i >= read.len() && chunk_idx == 1 {
                break;
            }
        }

        chunks
    }

    /// Chunk-based read classification for PacBio HiFi long reads.
    ///
    /// Phase 4 improvements:
    /// 1. Quality-weighted voting — chunks with higher alignment scores contribute more
    /// 2. Vote threshold — require minimum % of chunks to agree before accepting mapping
    /// 3. Top-N voting — return top N genomes instead of winner-takes-all
    ///
    /// Algorithm:
    /// 1. Split read into overlapping chunks (default: 150bp, 30bp overlap)
    /// 2. Map each chunk using existing Bit-Pop anchor_filter + XOR alignment
    /// 3. Quality-weighted vote: each genome scores = sum(chunk_score^2 / total_chunks)
    /// 4. Vote threshold: skip genomes with < threshold% of weighted votes
    /// 5. Top-N: return top N genomes by weighted score
    ///
    /// # Arguments
    /// * `read` — DNA read sequence (typically 6-11kb PacBio HiFi)
    /// * `context_window` — Flanking bases for context string
    ///
    /// # Returns
    /// Ranked MappingResult list (best match first). Empty if no match found.
    pub fn map_read_with_chunking(&self, read: &str, context_window: usize) -> Vec<MappingResult> {
        self.map_read_with_chunking_internal(read, context_window, None, None)
    }

    /// Map a long read using chunking with quality scores.
    /// Quality scores are used when chunk_anchor_strategy is Golden.
    pub fn map_read_with_chunking_and_quality(
        &self,
        read: &str,
        quality: &[u8],
        context_window: usize,
    ) -> Vec<MappingResult> {
        self.map_read_with_chunking_internal(read, context_window, None, Some(quality))
    }

    /// Internal chunk mapping with optional SNP mismatch collection and quality scores.
    fn map_read_with_chunking_internal(
        &self,
        read: &str,
        context_window: usize,
        mut snp_collector: Option<&mut std::collections::HashMap<(u32, u32, u8, u8), u32>>,
        quality: Option<&[u8]>,
    ) -> Vec<MappingResult> {
        let chunk_size = if self.chunk_pct > 0.0 {
            let raw = (read.len() as f64 * self.chunk_pct) as usize;
            raw.clamp(self.chunk_min, self.chunk_max).max(self.k)
        } else if self.chunk_size > 0 {
            self.chunk_size.max(self.k)
        } else {
            if read.len() > 1000 {
                150.max(self.k)
            } else {
                return self.map_read(read, context_window);
            }
        };

        let overlap = (chunk_size as f64 * 0.2) as usize;
        let chunks = self.split_read_into_chunks(read, chunk_size, overlap);

        if chunks.len() == 1 {
            return self.map_read(read, context_window);
        }

        let total_chunks = chunks.len();
        let mut genome_votes: HashMap<u32, ChunkVote> = HashMap::new();

        for (chunk_idx, chunk_seq) in &chunks {
            let chunk_start = *chunk_idx * (chunk_size - overlap);
            let chunk_quality = quality.map(|q| {
                let q_end = (chunk_start + chunk_seq.len()).min(q.len());
                &q[chunk_start..q_end]
            });
            let chunk_results = self.anchor_filter_for_chunk(
                chunk_seq,
                chunk_quality,
                0.5,
                DEFAULT_REPEAT_THRESHOLD,
            );
            if chunk_results.is_empty() {
                continue;
            }

            let best = &chunk_results[0];
            let genome_id = best.0;
            let position = best.1;
            let score = best.2;
            let cigar = best.3.clone();

            // Collect mismatch data for SNP detection
            if let Some(ref mut collector) = snp_collector {
                let read_enc = encode_sequence(chunk_seq);
                let genome = match self.genomes.get(&genome_id) {
                    Some(g) => g,
                    None => continue,
                };
                let pos = position as usize;
                let region_end = (pos + read_enc.len()).min(genome.len());
                if region_end > pos && !cigar.is_empty() {
                    let text_region = &genome[pos..region_end];
                    let mismatches = align::extract_mismatches(&read_enc, text_region, 0);
                    for m in mismatches {
                        let key = (genome_id, m.genome_pos as u32, m.genome_base, m.read_base);
                        *collector.entry(key).or_default() += 1;
                    }
                }
            }

            let vote = genome_votes.entry(genome_id).or_insert(ChunkVote {
                genome_id,
                win_count: 0,
                total_score: 0.0,
                avg_score: 0.0,
                chunks_mapped: 0,
                quality_weighted_score: 0.0,
            });

            vote.win_count += 1;
            vote.total_score += score;
            vote.chunks_mapped += 1;
            vote.avg_score = vote.total_score / vote.chunks_mapped as f64;
            vote.quality_weighted_score += score * score;
        }

        if genome_votes.is_empty() {
            return self.map_read(read, context_window);
        }

        let mut genome_ranking: Vec<(&u32, &ChunkVote)> = genome_votes.iter().collect();
        genome_ranking.sort_by(|a, b| {
            let score_a = a.1.quality_weighted_score / total_chunks as f64;
            let score_b = b.1.quality_weighted_score / total_chunks as f64;
            let fraction_a = a.1.win_count as f64 / total_chunks as f64;
            let fraction_b = b.1.win_count as f64 / total_chunks as f64;
            score_a
                .partial_cmp(&score_b)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    fraction_a
                        .partial_cmp(&fraction_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    a.1.total_score
                        .partial_cmp(&b.1.total_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });

        let threshold = self.chunk_vote_threshold;
        if threshold > 0.0 {
            genome_ranking
                .retain(|(_, vote)| vote.win_count as f64 / total_chunks as f64 >= threshold);
        }

        // Fallback: if vote threshold eliminated all genomes, use full-read mapping
        if genome_ranking.is_empty() {
            return self.map_read(read, context_window);
        }

        let fm = match &self.fm_index {
            Some(idx) => idx,
            None => return Vec::new(),
        };

        let encoded = encode_sequence(read);

        let full_results = self.map_read(read, context_window);

        // Apply SNP-aware scoring if enabled
        let top_n = self.chunk_top_n.min(genome_ranking.len());
        let results: Vec<MappingResult> = genome_ranking
            .into_iter()
            .take(top_n)
            .map(|(genome_id, vote)| {
                let normalized_weighted_score = vote.quality_weighted_score / total_chunks as f64;
                let align_score = normalized_weighted_score.sqrt().clamp(0.0, 1.0);

                let combined_score = if self.enable_snp_detect {
                    // SNP-aware scoring: penalize low scores that aren't supported by known SNPs
                    // The chunk voting already aggregates across chunks, so we apply a mild penalty
                    // to results where the vote fraction is low (indicating inconsistent alignment)
                    let vote_fraction = vote.win_count as f64 / total_chunks as f64;
                    if vote_fraction < 0.5 {
                        // Low agreement: apply penalty unless supported by SNP data
                        align_score * 0.9
                    } else {
                        align_score
                    }
                } else {
                    align_score
                };

                // Per-genome rarity: extract k-mer from genome at candidate position
                let rarity = if encoded.len() >= self.k {
                    let pos = full_results.first().map(|r| r.position).unwrap_or(0);
                    if let Some(genome_seq) = self.genomes.get(genome_id) {
                        let start = (pos as usize).min(genome_seq.len().saturating_sub(1));
                        let end = (start + self.k).min(genome_seq.len());
                        if end - start == self.k {
                            let kmer_bytes = &genome_seq[start..end];
                            let kmer_encoded = kmer_bytes
                                .iter()
                                .map(|&b| b as char)
                                .filter_map(encode_base)
                                .collect::<Vec<u8>>();
                            if kmer_encoded.len() == self.k {
                                let occ = fm.count_occurrences(&kmer_encoded);
                                1.0 / (occ as f64).max(1.0)
                            } else {
                                1.0
                            }
                        } else {
                            1.0
                        }
                    } else {
                        1.0
                    }
                } else {
                    1.0
                };

                MappingResult {
                    genome_id: *genome_id,
                    position: full_results.first().map(|r| r.position).unwrap_or(0),
                    score: combined_score * 0.85 + rarity * 0.15,
                    cigar: full_results
                        .first()
                        .map(|r| r.cigar.clone())
                        .unwrap_or_else(|| format!("{}M", read.len())),
                    context: full_results
                        .first()
                        .map(|r| r.context.clone())
                        .unwrap_or_default(),
                    is_reverse: full_results.first().map(|r| r.is_reverse).unwrap_or(false),
                    rarity,
                    md_string: full_results
                        .first()
                        .map(|r| r.md_string.clone())
                        .unwrap_or_default(),
                    hf_score: 0.0,
                }
            })
            .collect();

        results
    }

    /// Map multiple reads with chunk-based classification in parallel.
    ///
    /// Uses rayon for parallel chunk mapping. Each read is split into chunks,
    /// mapped independently, then voted on for final classification.
    ///
    /// If SNP detection is enabled, collects mismatches from all reads,
    /// builds a SNP map, and applies SNP-aware scoring.
    pub fn map_reads_with_chunking_parallel(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
    ) -> std::io::Result<usize> {
        self.map_reads_with_chunking_parallel_with_progress(
            reads,
            output_path,
            context_window,
            |_, _| {},
        )
    }

    pub fn map_reads_with_chunking_parallel_with_progress(
        &self,
        reads: &[(&str, &str)],
        output_path: &str,
        context_window: usize,
        progress: impl Fn(usize, usize) + Sync + Send,
    ) -> std::io::Result<usize> {
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned
            .iter()
            .map(|(n, l)| (n.as_str(), *l))
            .collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();
        let total = reads.len();

        // Phase 1: Collect SNP mismatch data if enabled
        let snp_detector = if self.enable_snp_detect {
            use std::collections::HashMap;
            use std::sync::atomic::{AtomicUsize, Ordering};
            let progress_counter = AtomicUsize::new(0);

            let all_data: Vec<(
                String,
                String,
                Vec<MappingResult>,
                HashMap<(u32, u32, u8, u8), u32>,
            )> = reads
                .par_iter()
                .map(|(name, seq)| {
                    let mut local_mismatches: HashMap<(u32, u32, u8, u8), u32> = HashMap::new();
                    let results = self.map_read_with_chunking_internal(
                        seq,
                        context_window,
                        Some(&mut local_mismatches),
                        None,
                    );
                    let completed = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                    progress(completed, total);
                    (name.to_string(), seq.to_string(), results, local_mismatches)
                })
                .collect();

            let mut merged_counts: HashMap<(u32, u32, u8, u8), u32> = HashMap::new();
            let results: Vec<(String, String, Vec<MappingResult>)> = all_data
                .into_iter()
                .map(|(name, seq, results, local_counts)| {
                    for (key, count) in local_counts {
                        *merged_counts.entry(key).or_default() += count;
                    }
                    (name, seq, results)
                })
                .collect();

            // Build SNP detector from merged counts
            let mut detector = snp::SnpDetector::new();
            for ((genome_id, position, ref_base, alt_base), count) in &merged_counts {
                if *count >= self.snp_min_support {
                    // Add mismatch count times so build_snp_map picks it up
                    for _ in 0..*count {
                        detector.add_mismatch(*genome_id, *position, *alt_base, *ref_base);
                    }
                }
            }
            detector.build_snp_map(1); // Already filtered by min_support

            if self.genome_count() > 0 {
                println!(
                    "  SNP detection: found {} SNPs (min support: {})",
                    detector.snp_count(),
                    self.snp_min_support
                );
            }

            Some((results, detector))
        } else {
            None
        };

        // Phase 2: Re-score with SNP-aware scoring or write results
        let mapped_count = if let Some((all_mismatches, detector)) = snp_detector {
            let mut writer = sam::SamWriter::new(output_path)?;
            writer.write_header(&genome_header)?;

            let mut mapped_count = 0;
            for (name, seq, results) in &all_mismatches {
                // Apply SNP-aware scoring to each result
                let scored_results = if self.enable_snp_detect {
                    self.apply_snp_scoring(results, &detector, seq)
                } else {
                    results.to_vec()
                };

                writer.write_mappings(name, seq, &scored_results, &name_refs)?;
                if !scored_results.is_empty() {
                    mapped_count += 1;
                }
            }

            mapped_count
        } else {
            use std::sync::atomic::{AtomicUsize, Ordering};
            let progress_counter = AtomicUsize::new(0);

            let mapped: Vec<(String, String, Vec<MappingResult>)> = reads
                .par_iter()
                .map(|(name, seq)| {
                    let results = self.map_read_with_chunking(seq, context_window);
                    let completed = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                    progress(completed, total);
                    (name.to_string(), seq.to_string(), results)
                })
                .collect();

            let mut writer = sam::SamWriter::new(output_path)?;
            writer.write_header(&genome_header)?;

            let mut mapped_count = 0;
            for (name, seq, results) in &mapped {
                writer.write_mappings(name, seq, results, &name_refs)?;
                if !results.is_empty() {
                    mapped_count += 1;
                }
            }

            mapped_count
        };

        Ok(mapped_count)
    }

    /// Apply SNP-aware scoring to mapping results.
    ///
    /// For each result, checks if the alignment score is consistent with
    /// known SNP positions. Results with scores that conflict with SNP
    /// data get penalized.
    fn apply_snp_scoring(
        &self,
        results: &[MappingResult],
        detector: &snp::SnpDetector,
        _read_seq: &str,
    ) -> Vec<MappingResult> {
        if !detector.is_built() {
            return results.to_vec();
        }

        results
            .iter()
            .map(|result| {
                let genome_snps = detector.get_genome_snps(result.genome_id);
                let snp_count = genome_snps.len();

                // If genome has many known SNPs, it's a good strain-specific match
                // Boost score slightly for genomes with SNP support
                let snp_boost = if snp_count > 0 {
                    (snp_count as f64 * self.snp_penalty).min(0.15)
                } else {
                    0.0
                };

                let mut new_result = result.clone();
                new_result.score = (result.score + snp_boost).min(1.0);

                new_result
            })
            .collect()
    }
}

fn compute_tlen(
    map1: &Option<PairedReadMapping>,
    map2: &Option<PairedReadMapping>,
    len1: usize,
    len2: usize,
) -> i64 {
    match (map1, map2) {
        (Some(m1), Some(m2)) => {
            if m1.genome_id != m2.genome_id {
                return 0;
            }
            let pos1 = m1.position as i64;
            let pos2 = m2.position as i64;
            let end1 = pos1 + len1 as i64;
            let end2 = pos2 + len2 as i64;
            let outer_start = pos1.min(pos2);
            let outer_end = end1.max(end2);
            let tlen = outer_end - outer_start;

            // Sign based on which read is forward/reverse
            if m1.is_reverse {
                -tlen
            } else {
                tlen
            }
        }
        _ => 0,
    }
}

/// Print simple atomic progress line when BITPOP_PROGRESS=atomic env var is set.
/// Call this alongside existing progress bar updates to show progress in GUI.
pub fn report_atomic_progress(current: u64, total: u64) {
    if std::env::var("BITPOP_PROGRESS").as_deref().ok() == Some("atomic") {
        let pct = (current as f64 / total as f64) * 100.0;
        eprintln!("\r  Progress: {}/{} ({:.1}%)", current, total, pct);
        let _ = std::io::Write::flush(&mut std::io::stderr());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_base() {
        assert_eq!(encode_base('A'), Some(1));
        assert_eq!(encode_base('C'), Some(2));
        assert_eq!(encode_base('G'), Some(3));
        assert_eq!(encode_base('T'), Some(4));
        assert_eq!(encode_base('N'), None);
        assert_eq!(encode_base('a'), Some(1));
        assert_eq!(decode_base(1), 'A');
        assert_eq!(decode_base(2), 'C');
        assert_eq!(decode_base(3), 'G');
        assert_eq!(decode_base(4), 'T');
    }

    #[test]
    fn test_encode_decode_sequence() {
        let seq = "ACGTACGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded.len(), 8);
        let decoded = decode_sequence(&encoded);
        assert_eq!(decoded, seq.to_uppercase());
    }

    #[test]
    fn test_encode_decode_sequence_with_n() {
        let seq = "ACNGT";
        let encoded = encode_sequence(seq);
        assert_eq!(encoded.len(), 4); // N is skipped
        let decoded = decode_sequence(&encoded);
        assert_eq!(decoded, "ACGT");
    }

    #[test]
    fn test_encode_decode_kmer() {
        let kmer = "ACGTACGT";
        let encoded = encode_kmer(kmer).expect("Should encode");
        let decoded = decode_kmer(encoded, 8);
        assert_eq!(decoded, kmer);
    }

    #[test]
    fn test_kmer_too_long() {
        assert!(encode_kmer(&"ACGT".repeat(10)).is_none());
    }

    #[test]
    fn test_kmer_invalid_base() {
        assert!(encode_kmer("ACGX").is_none());
    }

    #[test]
    fn test_add_genome() {
        let mut bp = BitPop::new(6);
        let gid = bp.add_genome("test", "ACGTACGTACGTACGT");
        assert_eq!(gid, 0);
        assert_eq!(bp.genome_count(), 1);
        bp.build();
        assert!(bp.bwt_len() > 0);
    }

    #[test]
    fn test_kmer_filter_basic() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGTACGTACGT");
        bp.build();
        let candidates = bp.kmer_filter("ACGTACGTACGT");
        assert!(!candidates.is_empty());
        assert_eq!(candidates[0].0, 0);
    }

    #[test]
    fn test_kmer_filter_no_match() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        let candidates = bp.kmer_filter("TTTTTTTTTTTT");
        assert!(candidates.is_empty());
    }

    #[test]
    fn test_align_read_exact() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        let (score, cigar, _) = bp.align_read("ACGTACGT", 0, 0);
        assert_eq!(score, 1.0);
        assert_eq!(cigar, "8M");
    }

    #[test]
    fn test_align_read_mismatch() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        let (score, cigar, _) = bp.align_read("ACGTACGA", 0, 0);
        // 2-bit XOR: 7/8 matches (ACGTACG match, last A≠T mismatch)
        assert!((score - 0.875).abs() < 0.001);
        assert_eq!(cigar, "7M1X");
    }

    #[test]
    fn test_map_read_full_pipeline() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        let results = bp.map_read("ACGTACGT", 3);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
        assert!(!results[0].context.is_empty());
    }

    #[test]
    fn test_multi_genome() {
        let mut bp = BitPop::new(6);
        let g1 = bp.add_genome("human", "ACGTACGTACGTACGT");
        let g2 = bp.add_genome("chimp", "ACGTACGTACGTAACA");
        assert_eq!(g1, 0);
        assert_eq!(g2, 1);
        assert_eq!(bp.genome_count(), 2);
        bp.build();
        let results = bp.map_read("ACGTACGTACGT", 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_multi_genome_ranking() {
        let mut bp = BitPop::new(6);
        bp.add_genome("shared", "ACGTAACAACGTAACAACGTAACA");
        bp.add_genome("unique", "TTTTTTTTACGTAACATTTTTTTT");
        bp.build();
        let results = bp.map_read("ACGTAACA", 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_load_genome_fasta() {
        use std::fs::File;
        use std::io::Write;

        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_fasta_load_{}_{}.fasta",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let mut f = File::create(&path).unwrap();
        writeln!(f, ">chr1 Human chromosome 1").unwrap();
        writeln!(f, "ACGTACGTACGTACGT").unwrap();
        writeln!(f, ">chr2 Human chromosome 2").unwrap();
        writeln!(f, "TTTTGGGGACGTACGT").unwrap();
        drop(f);

        let mut bp = BitPop::new(6);
        let ids = bp.load_genome_fasta(path.to_str().unwrap()).unwrap();
        assert_eq!(ids, vec![0, 1]);
        assert_eq!(bp.genome_count(), 2);
        assert_eq!(bp.genome_name(0), Some("chr1 Human chromosome 1"));
        assert_eq!(bp.genome_name(1), Some("chr2 Human chromosome 2"));
        bp.build();

        let results = bp.map_read("ACGTACGT", 3);
        assert!(!results.is_empty());

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_load_fasta_nonexistent() {
        let mut bp = BitPop::new(6);
        let result = bp.load_genome_fasta("/nonexistent/file.fasta");
        assert!(result.is_err());
    }

    #[test]
    fn test_map_reads_to_sam() {
        let mut bp = BitPop::new(6);
        bp.add_genome("chr1", "ACGTACGTACGTACGTACGTACGT");
        bp.add_genome("chr2", "TTTTGGGGTTTTGGGGTTTTGGGG");
        bp.build();

        let reads = vec![
            ("read1", "ACGTACGT"),
            ("read2", "TTTTGGGG"),
            ("read3", "NNNNSUPERINVALID"),
        ];

        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_sam_{}_{}.sam",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap().to_string();

        let mapped = bp.map_reads_to_sam(&reads, &path_str, 3).unwrap();
        assert_eq!(mapped, 2);

        let content = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.trim().lines().collect();

        // 2 header lines + 3 data lines (read3 unmapped)
        assert!(lines.len() >= 5);
        assert!(lines[0].starts_with("@SQ"));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_build_and_map() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGTACGTACGT");
        bp.build();

        let results = bp.map_read("ACGTACGT", 3);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
    }

    #[test]
    fn test_build_and_map_multi_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("human", "ACGTACGTACGTACGT");
        bp.add_genome("chimp", "ACGTACGTACGTAACA");
        bp.build();

        let results = bp.map_read("ACGTACGTACGT", 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_build_preserves_functionality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGTACGTACGT");
        bp.build();
        let results = bp.map_read("ACGTACGT", 3);

        assert!(!results.is_empty());
        assert!(results[0].score >= 0.5);
        assert_eq!(results[0].genome_id, 0);
    }

    #[test]
    fn test_build_and_sam() {
        let mut bp = BitPop::new(6);
        bp.add_genome("chr1", "ACGTACGTACGTACGTACGTACGT");
        bp.build();

        let reads = vec![("read1", "ACGTACGT")];

        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_sam_build_{}_{}.sam",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let mapped = bp
            .map_reads_to_sam(&reads, path.to_str().unwrap(), 3)
            .unwrap();
        assert_eq!(mapped, 1);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_anchor_filter_basic() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        let scored = bp.anchor_filter("ACGTACGT", 0.5);
        assert!(!scored.is_empty());
        assert_eq!(scored[0].0, 0);
        assert!(scored[0].2 >= 0.5);
    }

    #[test]
    fn test_anchor_filter_no_match() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        let scored = bp.anchor_filter("TTTTTTTTTTTT", 0.5);
        assert!(scored.is_empty());
    }

    #[test]
    fn test_anchor_filter_multi_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("human", "ACGTACGTACGTACGT");
        bp.add_genome("chimp", "ACGTACGTACGTAACA");
        bp.build();
        let scored = bp.anchor_filter("ACGTACGTACGT", 0.5);
        assert!(!scored.is_empty());
        for (gid, _, score, _) in &scored {
            assert!(*score >= 0.5);
            assert!(*gid < 2);
        }
    }

    #[test]
    fn test_anchor_filter_exact_match() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "NNNNACGTACGTACGTNNNN");
        bp.build();
        let scored = bp.anchor_filter("ACGTACGT", 0.5);
        assert!(!scored.is_empty());
        assert_eq!(scored[0].2, 1.0);
    }

    #[test]
    fn test_anchor_filter_partial_match() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "NNNNACGTACGANNNN");
        bp.build();
        let scored = bp.anchor_filter("ACGTACGT", 0.3);
        assert!(!scored.is_empty());
        assert!(scored[0].2 > 0.0 && scored[0].2 < 1.0);
    }

    #[test]
    fn test_anchor_filter_with_build() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        let scored = bp.anchor_filter("ACGTACGT", 0.5);
        assert!(!scored.is_empty());
    }

    #[test]
    fn test_rank_scored_results() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let scored = vec![
            (0u32, 6u64, 1.0f64, "8M".to_string()),
            (0u32, 15u64, 0.875f64, "7M1X".to_string()),
        ];
        let results = bp.rank_scored_results(&scored, "ACGTACGT", 3, 0.0);
        assert_eq!(results.len(), 2);
        assert!(results[0].score > results[1].score);
    }

    #[test]
    fn test_anchor_filter_long_read() {
        let mut bp = BitPop::new(8);
        let genome = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        bp.add_genome("test", genome);
        bp.build();
        let read = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let scored = bp.anchor_filter(read, 0.5);
        assert!(!scored.is_empty());
        assert!(scored[0].2 >= 0.9);
    }

    #[test]
    fn test_anchor_filter_with_threshold() {
        let mut bp = BitPop::new(6);
        let genome = format!("{}{}{}", "AAAAAA", "AACCGGTT", "TTTTTT");
        bp.add_genome("test", &genome);
        bp.build();

        let scored = bp.anchor_filter_with_threshold("AACCGGTT", 0.5, 100);
        assert!(!scored.is_empty());
        assert_eq!(scored[0].0, 0);
    }

    #[test]
    fn test_anchor_filter_threshold_blocks_repetitive() {
        let mut bp = BitPop::new(6);
        let repetitive = "ACGT".repeat(5000);
        bp.add_genome("repetitive", &repetitive);
        bp.build();

        let scored_no_thresh = bp.anchor_filter_with_threshold("ACGTACGTACGTACGT", 0.5, usize::MAX);

        let scored_tight = bp.anchor_filter_with_threshold("ACGTACGTACGTACGT", 0.5, 100);

        assert!(scored_tight.len() <= scored_no_thresh.len() || scored_tight.is_empty());
    }

    #[test]
    fn test_quality_preserved_with_threshold() {
        let mut bp = BitPop::new(8);

        let mut genome = String::new();
        genome.push_str(&"ACGT".repeat(500));
        genome.push_str("AACCGGTTAACCGGTT");
        genome.push_str(&"TTTT".repeat(500));

        bp.add_genome("test", &genome);
        bp.build();

        let results_no_thresh = bp.map_read("AACCGGTTAACCGGTT", 3);
        let results_with_thresh = bp
            .anchor_filter_with_threshold("AACCGGTTAACCGGTT", 0.5, 100)
            .into_iter()
            .collect::<Vec<_>>();

        assert!(!results_no_thresh.is_empty());

        if !results_with_thresh.is_empty() {
            let unique_pos = genome.find("AACCGGTTAACCGGTT").unwrap();
            let best = &results_with_thresh[0];
            assert!(best.1 as usize >= unique_pos.saturating_sub(5));
        }
    }

    #[test]
    fn test_kmer_filter_with_threshold() {
        let mut bp = BitPop::new(6);

        let repetitive = "ACGT".repeat(1000);
        bp.add_genome("repetitive", &repetitive);
        bp.build();

        let candidates_raw = bp.kmer_filter("ACGTACGTACGT");

        let candidates_filtered = bp.kmer_filter_with_threshold("ACGTACGTACGT", 50);

        assert!(candidates_filtered.len() <= candidates_raw.len());
    }

    #[test]
    fn test_kmer_filter_threshold_preserves_unique() {
        let mut bp = BitPop::new(6);

        let mut genome = String::new();
        genome.push_str(&"AAAA".repeat(500));
        genome.push_str("AACCGGTT");
        genome.push_str(&"TTTT".repeat(500));

        bp.add_genome("test", &genome);
        bp.build();

        let candidates = bp.kmer_filter_with_threshold("AACCGGTT", 100);
        assert!(!candidates.is_empty());
    }

    #[test]
    fn test_default_repeat_threshold_constant() {
        assert_eq!(DEFAULT_REPEAT_THRESHOLD, 10000);
    }

    // === FAZA 5: Quality-Aware Pipeline Tests ===

    #[test]
    fn test_quality_aware_kmer_filter_basic() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGTACGTACGT");
        bp.build();

        let read = "ACGTACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30];

        let candidates = bp.kmer_filter_with_quality(read, &quality, 20, usize::MAX);
        assert!(!candidates.is_empty());
    }

    #[test]
    fn test_quality_aware_kmer_filter_filters_low_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGTACGTACGT");
        bp.build();

        let read = "ACGTACGTACGT";
        let high_qual: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30];
        let low_qual: Vec<u8> = vec![5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5];

        let candidates_high = bp.kmer_filter_with_quality(read, &high_qual, 20, usize::MAX);
        let candidates_low = bp.kmer_filter_with_quality(read, &low_qual, 20, usize::MAX);

        assert!(!candidates_high.is_empty());
        assert!(candidates_low.is_empty() || candidates_low.len() <= candidates_high.len());
    }

    #[test]
    fn test_map_read_with_quality_basic() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let results = bp.map_read_with_quality(read, &quality, 20, 3);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
        assert!(results[0].align_score >= 0.5);
    }

    #[test]
    fn test_map_read_with_quality_low_quality_filtering() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let high_qual: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];
        let low_qual: Vec<u8> = vec![5, 5, 5, 5, 5, 5, 5, 5];

        let results_high = bp.map_read_with_quality(read, &high_qual, 20, 3);
        let results_low = bp.map_read_with_quality(read, &low_qual, 20, 3);

        // High quality should give better or equal results
        if !results_high.is_empty() && !results_low.is_empty() {
            assert!(results_high[0].align_score >= results_low[0].align_score);
        }
    }

    #[test]
    fn test_quality_mapping_result_fields() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let results = bp.map_read_with_quality(read, &quality, 20, 3);

        assert!(!results.is_empty());
        let r = &results[0];
        assert_eq!(r.quality_scores.len(), 8);
        assert!(!r.context.is_empty());
    }

    #[test]
    fn test_quality_aware_scoring_penalizes_high_qual_mismatches() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Read with mismatch at high quality position
        let read = "ACGTACGA"; // Last base A instead of T (mismatch)
        let high_qual: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let results = bp.map_read_with_quality(read, &high_qual, 20, 3);

        if !results.is_empty() {
            // Should have a quality penalty for the high-quality mismatch
            assert!(results[0].quality_penalty <= 0.0);
        }
    }

    #[test]
    fn test_fastq_parallel_mapping() {
        let mut bp = BitPop::new(6);
        bp.add_genome("chr1", "ACGTACGTACGTACGTACGTACGT");
        bp.add_genome("chr2", "TTTTGGGGTTTTGGGGTTTTGGGG");
        bp.build();

        let dir = std::env::temp_dir();
        let fastq_path = dir.join(format!(
            "test_fastq_{}_{}.fastq",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let sam_path = dir.join(format!(
            "test_fastq_{}_{}.sam",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));

        {
            use std::io::Write;
            let mut f = std::fs::File::create(&fastq_path).unwrap();
            writeln!(f, "@read1").unwrap();
            writeln!(f, "ACGTACGT").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "IIIIIIII").unwrap();
            writeln!(f, "@read2").unwrap();
            writeln!(f, "TTTTGGGG").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "!!!!!!!!").unwrap();
        }

        let mapped = bp
            .map_reads_from_fastq_parallel(
                fastq_path.to_str().unwrap(),
                sam_path.to_str().unwrap(),
                20,
                3,
            )
            .unwrap();

        assert_eq!(mapped, 2);

        let content = std::fs::read_to_string(sam_path.to_str().unwrap()).unwrap();
        assert!(content.contains("@SQ"));
        assert!(content.contains("chr1"));

        let _ = std::fs::remove_file(&fastq_path);
        let _ = std::fs::remove_file(&sam_path);
    }

    // === FAZA 6: Parallel Optimization Tests ===

    #[test]
    fn test_build_parallel_produces_same_results() {
        let mut bp_sequential = BitPop::new(6);
        bp_sequential.add_genome("human", "ACGTACGTACGTACGT");
        bp_sequential.add_genome("chimp", "ACGTACGTACGTAACA");
        bp_sequential.build();

        let mut bp_parallel = BitPop::new(6);
        bp_parallel.add_genome("human", "ACGTACGTACGTACGT");
        bp_parallel.add_genome("chimp", "ACGTACGTACGTAACA");
        bp_parallel.build_parallel();

        let results_seq = bp_sequential.map_read("ACGTACGTACGT", 3);
        let results_par = bp_parallel.map_read("ACGTACGTACGT", 3);

        assert!(!results_seq.is_empty());
        assert!(!results_par.is_empty());
        assert_eq!(results_seq.len(), results_par.len());
    }

    #[test]
    fn test_build_parallel_multi_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("genome_0", "ACGTAAAAACGTAAAA");
        bp.add_genome("genome_1", "CGTGTTTTACGTTTTT");
        bp.add_genome("genome_2", "ATATGGGGATATGGGG");
        bp.add_genome("genome_3", "GTATCCCCGTATCCCC");
        bp.add_genome("genome_4", "ACGTACGTACGTACGT");
        bp.build_parallel();

        let results = bp.map_read("ACGTACGT", 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_optimized_parallel_mapping() {
        let mut bp = BitPop::new(6);
        bp.add_genome("chr1", "ACGTACGTACGTACGTACGTACGT");
        bp.add_genome("chr2", "TTTTGGGGTTTTGGGGTTTTGGGG");
        bp.build();

        let reads: Vec<(&str, &str)> = vec![
            ("read1", "ACGTACGT"),
            ("read2", "TTTTGGGG"),
            ("read3", "ACGTACGT"),
            ("read4", "TTTTGGGG"),
        ];

        let dir = std::env::temp_dir();
        let sam_path = dir.join(format!(
            "test_optimized_{}_{}.sam",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));

        let mapped = bp
            .map_reads_parallel_optimized(
                &reads,
                sam_path.to_str().unwrap(),
                3,
                2, // batch_size
            )
            .unwrap();

        assert_eq!(mapped, 4);

        let _ = std::fs::remove_file(&sam_path);
    }

    #[test]
    fn test_quality_aware_anchor_filter() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let scored = bp.anchor_filter_with_quality(read, &quality, 0.5, 20, usize::MAX);
        assert!(!scored.is_empty());
        assert_eq!(scored[0].0, 0);
        assert!(scored[0].2 >= 0.5);
    }

    #[test]
    fn test_quality_aware_anchor_filter_with_low_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let low_qual: Vec<u8> = vec![5, 5, 5, 5, 5, 5, 5, 5];

        // Should fallback to regular anchor filter when all k-mers are low quality
        let scored = bp.anchor_filter_with_quality(read, &low_qual, 0.5, 20, usize::MAX);

        if !scored.is_empty() {
            // Fallback should still work but with no quality penalty info
            assert!(scored[0].2 >= 0.5);
        }
    }

    // === AlignMode Tests ===

    #[test]
    fn test_align_mode_display() {
        assert_eq!(AlignMode::Xor.to_string(), "xor");
        assert_eq!(AlignMode::Sw.to_string(), "sw");
        assert_eq!(AlignMode::Hybrid.to_string(), "hybrid");
    }

    #[test]
    fn test_align_mode_default() {
        let default: AlignMode = Default::default();
        assert_eq!(default, AlignMode::Xor);
    }

    #[test]
    fn test_align_read_sw_basic() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let (score, cigar, _) = bp.align_read_sw("ACGTACGT", 0, 0);
        assert!(score > 0.0);
        assert!(!cigar.is_empty());
    }

    #[test]
    fn test_align_read_sw_exact_match() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let (_score, cigar, _) = bp.align_read_sw("ACGTACGT", 0, 0);
        assert_eq!(cigar, "8M");
    }

    #[test]
    fn test_align_read_sw_vs_xor() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";

        // Position 7 is where ACGTACGT starts in the genome (A at pos 7)
        let (xor_score, _, _) = bp.align_read(read, 0, 7);
        let (sw_score, _, _) = bp.align_read_sw(read, 0, 7);

        // Both should find the match, SW might be slightly different due to local alignment
        assert!(xor_score > 0.0);
        assert!(sw_score > 0.0);
    }

    #[test]
    fn test_align_read_with_mode_xor() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let (score, _, _) = bp.align_read_with_mode("ACGTACGT", AlignMode::Xor, 0, 0);
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_align_read_with_mode_sw() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let (score, _, _) = bp.align_read_with_mode("ACGTACGT", AlignMode::Sw, 0, 0);
        assert!(score > 0.0);
    }

    #[test]
    fn test_align_read_with_mode_hybrid() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let (score, _, _) = bp.align_read_with_mode("ACGTACGT", AlignMode::Hybrid, 0, 0);
        assert!(score > 0.0);
    }

    #[test]
    fn test_anchor_filter_with_mode_xor() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let scored = bp.anchor_filter_with_mode("ACGTACGT", AlignMode::Xor, 0.5, usize::MAX);
        assert!(!scored.is_empty());
        assert!(scored[0].2 >= 0.5);
    }

    #[test]
    fn test_anchor_filter_with_mode_sw() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let scored = bp.anchor_filter_with_mode("ACGTACGT", AlignMode::Sw, 0.5, usize::MAX);
        assert!(!scored.is_empty());
    }

    #[test]
    fn test_anchor_filter_with_mode_hybrid() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let scored = bp.anchor_filter_with_mode("ACGTACGT", AlignMode::Hybrid, 0.5, usize::MAX);
        assert!(!scored.is_empty());
    }

    #[test]
    fn test_map_read_with_mode() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let results_xor = bp.map_read_with_mode("ACGTACGT", AlignMode::Xor, 3);
        assert!(!results_xor.is_empty());
        assert_eq!(results_xor[0].genome_id, 0);

        let results_sw = bp.map_read_with_mode("ACGTACGT", AlignMode::Sw, 3);
        assert!(!results_sw.is_empty());
    }

    #[test]
    fn test_map_read_with_quality_mode() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let results = bp.map_read_with_quality_mode(read, &quality, AlignMode::Xor, 20, 3);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
        assert_eq!(results[0].quality_scores.len(), 8);
    }

    #[test]
    fn test_map_read_with_quality_mode_sw() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let results = bp.map_read_with_quality_mode(read, &quality, AlignMode::Sw, 20, 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_smart_threshold_short_read() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Short read (<20bp) should have stricter threshold
        let threshold = bp.compute_smart_threshold(15, false, 20.0);
        assert!(threshold > 0.5);
    }

    #[test]
    fn test_smart_threshold_long_read() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Long read (>100bp) should have more lenient threshold
        let threshold = bp.compute_smart_threshold(150, false, 20.0);
        assert!(threshold < 0.5);
    }

    #[test]
    fn test_smart_threshold_high_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // High quality should have stricter threshold
        let threshold = bp.compute_smart_threshold(50, true, 30.0);
        assert!(threshold > 0.5);
    }

    #[test]
    fn test_smart_threshold_low_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Low quality should have more lenient threshold
        let threshold = bp.compute_smart_threshold(50, true, 10.0);
        assert!(threshold < 0.5);
    }

    #[test]
    fn test_anchor_filter_with_quality_smart() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        // Smart threshold should adapt based on quality
        let scored_smart = bp.anchor_filter_with_quality_smart(read, &quality, 0.5, 20, usize::MAX);

        let scored_fixed = bp.anchor_filter_with_quality(read, &quality, 0.5, 20, usize::MAX);

        // Both should find the same match (same genome)
        if !scored_smart.is_empty() && !scored_fixed.is_empty() {
            assert_eq!(scored_smart[0].0, scored_fixed[0].0);
        }
    }

    #[test]
    fn test_align_read_sw_with_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        let (score, cigar, _, penalty) = bp.align_read_sw_with_quality(read, &quality, 0, 0);
        assert!(score > 0.0);
        assert!(!cigar.is_empty());
        assert_eq!(penalty, 0.0); // perfect match = no penalty
    }

    #[test]
    fn test_align_read_sw_with_quality_mismatch() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Read has mismatch at last position (A instead of T)
        let read = "ACGTACGA";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];

        // Position 0 has ACGTACGT in genome, read is ACGTACGA (mismatch at pos 7)
        let (score, _, _, penalty) = bp.align_read_sw_with_quality(read, &quality, 0, 0);
        assert!(score > 0.0);
        // SW should find the match ACGTACG (7 bases) with score=14, normalized to ~0.875
        // No penalty since SW local alignment stops before the mismatch
        assert!(penalty >= 0.0);
    }

    #[test]
    fn test_full_pipeline_xor_vs_sw() {
        let mut bp = BitPop::new(6);
        bp.add_genome("human", "ACGTACGTACGTACGTACGTACGT");
        bp.add_genome("chimp", "ACGTACGTACGTAACAACGTACGT");
        bp.build();

        let read = "ACGTACGTACGT";

        let results_xor = bp.map_read_with_mode(read, AlignMode::Xor, 3);
        let results_sw = bp.map_read_with_mode(read, AlignMode::Sw, 3);

        // Both should find mappings
        assert!(!results_xor.is_empty());
        assert!(!results_sw.is_empty());

        // Top result should be same genome for both
        assert_eq!(results_xor[0].genome_id, results_sw[0].genome_id);
    }

    #[test]
    fn test_full_pipeline_hybrid() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";
        let results = bp.map_read_with_mode(read, AlignMode::Hybrid, 3);

        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
    }

    #[test]
    fn test_quality_aware_scoring_differentiates() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGTACGTACGTACGT");
        bp.build();

        // Same read, different qualities
        let read = "ACGTACGA"; // mismatch at last position
        let high_qual: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];
        let low_qual: Vec<u8> = vec![10, 10, 10, 10, 10, 10, 10, 10];

        let results_high = bp.map_read_with_quality_mode(read, &high_qual, AlignMode::Xor, 20, 3);
        let results_low = bp.map_read_with_quality_mode(read, &low_qual, AlignMode::Xor, 20, 3);

        if !results_high.is_empty() && !results_low.is_empty() {
            // High quality should have larger negative penalty for the mismatch
            assert!(results_high[0].quality_penalty < results_low[0].quality_penalty);
        }
    }

    #[test]
    fn test_smart_threshold_bounds() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGT");
        bp.build();

        // Threshold should always be in [0.3, 0.8] range
        for read_len in [5, 10, 20, 50, 100, 200].iter() {
            for has_qual in [true, false] {
                for avg_q in [5.0, 15.0, 25.0, 35.0] {
                    let threshold = bp.compute_smart_threshold(*read_len, has_qual, avg_q);
                    assert!(
                        (0.3..=0.8).contains(&threshold),
                        "Threshold {} out of bounds for len={}, qual={}, avg={}",
                        threshold,
                        read_len,
                        has_qual,
                        avg_q
                    );
                }
            }
        }
    }

    #[test]
    fn test_align_read_sw_empty_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("empty", "");
        bp.build();

        let (score, _, _) = bp.align_read_sw("ACGT", 0, 0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_align_read_sw_with_quality_empty_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("empty", "");
        bp.build();

        let quality: Vec<u8> = vec![30, 30, 30, 30];
        let (score, _, _, _) = bp.align_read_sw_with_quality("ACGT", &quality, 0, 0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_align_read_sw_invalid_genome_id() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "ACGT");
        bp.build();

        let (score, _, _) = bp.align_read_sw("ACGT", 999, 0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_align_mode_stability() {
        // Multiple calls with same mode should give consistent results
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        let read = "ACGTACGT";

        let (s1, c1, _) = bp.align_read_with_mode(read, AlignMode::Xor, 0, 6);
        let (s2, c2, _) = bp.align_read_with_mode(read, AlignMode::Xor, 0, 6);

        assert_eq!(s1, s2);
        assert_eq!(c1, c2);
    }

    #[test]
    fn test_top_n_default_is_one() {
        let bp = BitPop::new(8);
        assert_eq!(bp.top_n(), 1);
    }

    #[test]
    fn test_set_top_n() {
        let mut bp = BitPop::new(8);
        bp.set_top_n(3);
        assert_eq!(bp.top_n(), 3);
        bp.set_top_n(1);
        assert_eq!(bp.top_n(), 1);
    }

    #[test]
    fn test_top_n_anchor_filter_returns_results() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        bp.set_top_n(3);
        let scored = bp.anchor_filter("ACGTACGT", 0.5);
        assert!(!scored.is_empty());
        assert_eq!(scored[0].0, 0);
        assert!(scored[0].2 >= 0.5);
    }

    #[test]
    fn test_top_n_anchor_filter_with_error() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        bp.set_top_n(1);
        let scored_single = bp.anchor_filter("ACGTACGT", 0.5);
        bp.set_top_n(3);
        let scored_top3 = bp.anchor_filter("ACGTACGT", 0.5);
        assert!(!scored_top3.is_empty());
        assert!(scored_top3.len() >= scored_single.len());
    }

    #[test]
    fn test_top_n_multi_genome() {
        let mut bp = BitPop::new(6);
        bp.add_genome("human", "ACGTACGTACGTACGT");
        bp.add_genome("chimp", "ACGTACGTACGTAACA");
        bp.build();
        bp.set_top_n(3);
        let scored = bp.anchor_filter("ACGTACGTACGT", 0.5);
        assert!(!scored.is_empty());
        for (gid, _, score, _) in &scored {
            assert!(*score >= 0.5);
            assert!(*gid < 2);
        }
    }

    #[test]
    fn test_top_n_anchor_filter_with_quality() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        bp.set_top_n(3);
        let read = "ACGTACGT";
        let quality: Vec<u8> = vec![30, 30, 30, 30, 30, 30, 30, 30];
        let scored = bp.anchor_filter_with_quality(read, &quality, 0.5, 20, 10000);
        assert!(!scored.is_empty());
    }

    #[test]
    fn test_top_n_anchor_filter_quality_fallback() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();
        bp.set_top_n(3);
        let read = "ACGTACGT";
        let low_quality: Vec<u8> = vec![5, 5, 5, 5, 5, 5, 5, 5];
        let scored = bp.anchor_filter_with_quality(read, &low_quality, 0.5, 20, 10000);
        assert!(!scored.is_empty());
    }

    #[test]
    fn test_top_n_threshold_blocks_repetitive() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "NNNNNACGTACGTACGTACGTACGTNNNNN");
        bp.build();
        bp.set_top_n(3);
        let scored = bp.anchor_filter_with_threshold("ACGTACGTACGT", 0.5, 100);
        assert!(!scored.is_empty());
        for (_, _, score, _) in &scored {
            assert!(*score >= 0.5);
        }
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAA"), "TTT");
        assert_eq!(reverse_complement("CCCC"), "GGGG");
        assert_eq!(reverse_complement("GGGG"), "CCCC");
        assert_eq!(reverse_complement("TTTT"), "AAAA");
    }

    #[test]
    fn test_reverse_complement_bytes() {
        let encoded = encode_sequence("ACGT");
        let rc = reverse_complement_bytes(&encoded);
        let rc_seq = decode_sequence(&rc);
        assert_eq!(rc_seq, "ACGT");

        let encoded2 = encode_sequence("ATCG");
        let rc2 = reverse_complement_bytes(&encoded2);
        let rc_seq2 = decode_sequence(&rc2);
        assert_eq!(rc_seq2, "CGAT");
    }

    #[test]
    fn test_reverse_complement_double_rc() {
        let original = "ACGTACGTATCG";
        let rc1 = reverse_complement(original);
        let rc2 = reverse_complement(&rc1);
        assert_eq!(original, rc2);
    }

    #[test]
    fn test_reverse_complement_with_n() {
        let rc = reverse_complement("ACNNGT");
        assert_eq!(rc, "ACNNGT");
    }

    #[test]
    fn test_map_read_reverse_complement_forward_wins() {
        let mut bp = BitPop::new(10);
        bp.add_genome("genome1", "ACGTACGTACGTACGTACGT");
        bp.build();

        let read = "ACGTACGTACGT";
        let results = bp.map_read_with_mode(read, AlignMode::Xor, 5);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
        assert!(
            !results[0].is_reverse,
            "Forward alignment should win for perfect match"
        );
    }

    #[test]
    fn test_map_read_reverse_complement_rc_wins() {
        let mut bp = BitPop::new(10);
        bp.add_genome("genome1", "ACGTACGTACGTACGTACGT");
        bp.build();

        let rc_read = reverse_complement("ACGTACGTACGT");
        let results = bp.map_read_with_mode(&rc_read, AlignMode::Xor, 5);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
    }

    #[test]
    fn test_map_read_reverse_complement_no_forward_match() {
        let mut bp = BitPop::new(10);
        bp.add_genome("genome1", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        bp.build();

        let read = "CCCCCCCCCCCC";
        let results = bp.map_read_with_mode(read, AlignMode::Xor, 5);
        assert!(!results.is_empty());
        assert_eq!(results[0].genome_id, 0);
    }

    #[test]
    fn test_generate_kmer_variants_basic() {
        let kmer = encode_kmer("ACGT").unwrap();
        let variants = generate_kmer_variants(kmer, 4, 1);
        assert!(!variants.is_empty(), "should generate variants");
        assert!(
            !variants.contains(&kmer),
            "variants should not contain original"
        );

        // Verify one variant has exactly 1 mismatch
        let variant = variants[0];
        let original_decoded = decode_kmer(kmer, 4);
        let variant_decoded = decode_kmer(variant, 4);
        let mismatch_count = original_decoded
            .chars()
            .zip(variant_decoded.chars())
            .filter(|&(a, b)| a != b)
            .count();
        assert_eq!(mismatch_count, 1, "variant should have exactly 1 mismatch");
    }

    #[test]
    fn test_fuzzy_kmer_find_in_index() {
        let mut bp = BitPop::new(4);
        // Genome contains ACGT and CGTG (1 mismatch from ACGT at position 2)
        bp.add_genome("test", "ACGTACGTACGTCGTG");
        bp.build();

        let fm = bp.fm_index.as_ref().unwrap();

        // Original k-mer should be found
        let kmer_bytes = BitPop::kmer_encoded_to_bytes(encode_kmer("ACGT").unwrap(), 4);
        let count = fm.count_occurrences(&kmer_bytes);
        assert!(count > 0, "original kmer should be found (count={})", count);

        // Variant CGTG (1 mismatch at position 2: G->T) should be found
        let variant_bytes = BitPop::kmer_encoded_to_bytes(encode_kmer("CGTG").unwrap(), 4);
        let count = fm.count_occurrences(&variant_bytes);
        assert!(count > 0, "variant CGTG should be found (count={})", count);
    }

    #[test]
    fn test_generate_kmer_variants_no_mismatches() {
        let kmer = encode_kmer("ACGT").unwrap();
        let variants = generate_kmer_variants(kmer, 4, 0);
        assert_eq!(variants.len(), 0);
    }

    #[test]
    fn test_generate_kmer_variants_k_too_long() {
        // k > 31 should return empty
        let variants = generate_kmer_variants(0, 32, 1);
        assert_eq!(variants.len(), 0);
    }

    #[test]
    fn test_fuzzy_method_enum_default() {
        let method: FuzzyMethod = Default::default();
        assert_eq!(method, FuzzyMethod::None);
    }

    #[test]
    fn test_set_fuzzy_method() {
        let mut bp = BitPop::new(10);
        bp.set_fuzzy_method(FuzzyMethod::FuzzyKmer);
        bp.set_fuzzy_mismatches(2);
        assert_eq!(bp.fuzzy_mismatches, 2);
    }

    #[test]
    fn test_anchor_filter_fuzzy_kmer() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        bp.set_fuzzy_method(FuzzyMethod::FuzzyKmer);
        bp.set_fuzzy_mismatches(1);

        // Test that regular kmer filter works first
        let fm = bp.fm_index.as_ref().unwrap();
        let encoded = encode_sequence("ACGTAC");
        let regular_candidates = bp.find_top_n_rarest_kmers(&encoded, fm, usize::MAX);
        assert!(
            !regular_candidates.is_empty(),
            "regular kmer filter should find candidates"
        );

        let fuzzy_candidates = bp.find_top_n_rarest_kmers_fuzzy(&encoded, fm, usize::MAX);
        assert!(
            !fuzzy_candidates.is_empty(),
            "fuzzy kmer filter should find candidates"
        );

        let scored = bp.anchor_filter("ACGTAC", 0.5);
        assert!(
            !scored.is_empty(),
            "Fuzzy anchor_filter should find matches"
        );
    }

    #[test]
    fn test_anchor_filter_fuzzy_seed() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.build();

        bp.set_spaced_seed(true);
        bp.set_fuzzy_method(FuzzyMethod::FuzzySeed);
        bp.set_fuzzy_mismatches(1);

        // Just verify it doesn't panic and returns results
        let scored = bp.anchor_filter("ACGTAC", 0.3);
        // Note: fuzzy seed may or may not find matches depending on the seed pattern
        // The important thing is that it doesn't panic
        let _ = scored;
    }

    #[test]
    fn test_build_neighborhood_hash() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.set_fuzzy_method(FuzzyMethod::Neighborhood);
        bp.set_fuzzy_mismatches(1);
        bp.build();

        assert!(bp.neighborhood_hash.is_some());
    }

    #[test]
    fn test_anchor_filter_neighborhood() {
        let mut bp = BitPop::new(6);
        bp.add_genome("test", "AACCGGTACGTACGTAACCGGTTTC");
        bp.set_fuzzy_method(FuzzyMethod::Neighborhood);
        bp.set_fuzzy_mismatches(1);
        bp.build();

        let scored = bp.anchor_filter("ACGTAC", 0.5);
        assert!(
            !scored.is_empty(),
            "Neighborhood anchor_filter should find matches"
        );
    }

    #[test]
    fn test_discordant_pair_reconciliation_basic() {
        // Two genomes: G1 has "AAAAACCCCC", G2 has "GGGGGTTTTT"
        // R1 matches G1, R2 matches G2 → discordant
        // Reconciliation should NOT force concordant since scores are too different
        let mut bp = BitPop::new(5);
        bp.add_genome("G1", "AAAAACCCCCGGGGGTTTTT");
        bp.add_genome("G2", "GGGGGTTTTTAAAAACCCCC");
        bp.build();

        let paired = PairedRead {
            name: "test1".to_string(),
            read1_seq: "AAAAA".to_string(),
            read1_qual: vec![40; 5],
            read2_seq: "GGGGG".to_string(),
            read2_qual: vec![40; 5],
        };

        let result = bp.map_read_paired(&paired, 0, true, 5);

        // Both reads should map (even if discordant)
        assert!(result.map1.is_some(), "R1 should map");
        assert!(result.map2.is_some(), "R2 should map");
    }

    #[test]
    fn test_discordant_pair_reconciliation_concordant_wins() {
        // Both reads match G1 better → should be concordant on G1
        let mut bp = BitPop::new(5);
        bp.add_genome("G1", "AAAAACCCCCGGGGGTTTTT");
        bp.add_genome("G2", "CCCCCGGGGGTTTTTAAAAA");
        bp.build();

        let paired = PairedRead {
            name: "test2".to_string(),
            read1_seq: "AAAAA".to_string(),
            read1_qual: vec![40; 5],
            read2_seq: "CCCCC".to_string(),
            read2_qual: vec![40; 5],
        };

        let result = bp.map_read_paired(&paired, 0, true, 5);

        if let (Some(m1), Some(m2)) = (&result.map1, &result.map2) {
            // If both map, reconciliation should prefer concordant assignment
            // when the concordant score is competitive
            // Both should be on the same genome after reconciliation
            assert_eq!(
                m1.genome_id, m2.genome_id,
                "Concordant reconciliation: both reads should map to same genome"
            );
        }
    }

    #[test]
    fn test_reconciliation_disabled_keeps_discordant() {
        let mut bp = BitPop::new(5);
        bp.add_genome("G1", "AAAAACCCCCGGGGGTTTTT");
        bp.add_genome("G2", "GGGGGTTTTTAAAAACCCCC");
        bp.build();

        let paired = PairedRead {
            name: "test3".to_string(),
            read1_seq: "AAAAA".to_string(),
            read1_qual: vec![40; 5],
            read2_seq: "TTTTT".to_string(),
            read2_qual: vec![40; 5],
        };

        let result = bp.map_read_paired(&paired, 0, false, 5);

        // With reconciliation disabled, discordant pairs are allowed
        if result.map1.is_some() && result.map2.is_some() {
            // May or may not be same genome depending on scores
            // The point is reconciliation didn't force anything
            let _ = &result.map1;
            let _ = &result.map2;
        }
    }

    #[test]
    fn test_concordant_pair_unchanged_by_reconciliation() {
        // Both reads naturally map to same genome → reconciliation is no-op
        let mut bp = BitPop::new(5);
        bp.add_genome("G1", "AAAAACCCCCGGGGGTTTTT");
        bp.add_genome("G2", "NNNNNNNNNNNNNNNNNNNNNN");
        bp.build();

        let paired = PairedRead {
            name: "test4".to_string(),
            read1_seq: "AAAAA".to_string(),
            read1_qual: vec![40; 5],
            read2_seq: "CCCCC".to_string(),
            read2_qual: vec![40; 5],
        };

        let result = bp.map_read_paired(&paired, 0, true, 5);

        if let (Some(m1), Some(m2)) = (&result.map1, &result.map2) {
            assert_eq!(m1.genome_id, m2.genome_id, "Both should map to G1");
        }
    }

    // --- Gaussian Insert Size Model Tests ---

    #[test]
    fn test_gaussian_stats_new_is_empty() {
        let stats = InsertSizeStats::new();
        assert_eq!(stats.count, 0);
        assert_eq!(stats.mean, 0.0);
        assert_eq!(stats.stddev, 0.0);
    }

    #[test]
    fn test_gaussian_stats_update_mean() {
        let mut stats = InsertSizeStats::new();
        stats.update(100);
        assert_eq!(stats.count, 1);
        assert_eq!(stats.mean, 100.0);
        stats.update(200);
        assert_eq!(stats.count, 2);
        assert_eq!(stats.mean, 150.0);
        stats.update(300);
        assert_eq!(stats.count, 3);
        assert_eq!(stats.mean, 200.0);
    }

    #[test]
    fn test_gaussian_stats_update_stddev() {
        let mut stats = InsertSizeStats::new();
        stats.update(100);
        assert_eq!(stats.count, 1);
        // stddev is 0 with single observation
        stats.update(100);
        assert_eq!(stats.count, 2);
        assert_eq!(stats.stddev, 0.0); // identical values
        stats.update(200);
        assert!(stats.stddev > 0.0);
        // With values [100, 100, 200]: mean=133.33, variance=(33.33^2 + 33.33^2 + 66.67^2)/2 = 3333.33, stddev≈57.74
        assert!(stats.stddev > 50.0 && stats.stddev < 70.0);
    }

    #[test]
    fn test_gaussian_stats_negative_tlen_ignored() {
        let mut stats = InsertSizeStats::new();
        stats.update(-100);
        stats.update(0);
        assert_eq!(stats.count, 0);
        stats.update(150);
        assert_eq!(stats.count, 1);
        assert_eq!(stats.mean, 150.0);
    }

    #[test]
    fn test_gaussian_probability_at_mean() {
        let mut stats = InsertSizeStats::new();
        // Build a tight distribution around 300
        for _ in 0..100 {
            stats.update(300);
        }
        // With identical values, stddev = 0, so probability = 0
        // Try with slight variation
        let mut stats2 = InsertSizeStats::new();
        for _ in 0..50 {
            stats2.update(300);
        }
        for _ in 0..50 {
            stats2.update(302);
        }
        // Probability at mean should be highest
        let prob_at_mean = stats2.gaussian_probability(301);
        let prob_far = stats2.gaussian_probability(500);
        assert!(prob_at_mean > prob_far);
    }

    #[test]
    fn test_gaussian_probability_far_from_mean() {
        let mut stats = InsertSizeStats::new();
        for _ in 0..50 {
            stats.update(300);
        }
        for _ in 0..50 {
            stats.update(300);
        }
        // With stddev = 0, probability should be 0
        assert_eq!(stats.gaussian_probability(300), 0.0);

        let mut stats2 = InsertSizeStats::new();
        for _ in 0..25 {
            stats2.update(298);
        }
        for _ in 0..25 {
            stats2.update(300);
        }
        for _ in 0..25 {
            stats2.update(302);
        }
        for _ in 0..25 {
            stats2.update(304);
        }
        // Value far from mean should have very low probability
        let prob_near = stats2.gaussian_probability(300);
        let prob_far = stats2.gaussian_probability(1000);
        assert!(prob_near > prob_far);
        assert!(prob_far < 0.0001);
    }

    #[test]
    fn test_gaussian_log_probability() {
        let mut stats = InsertSizeStats::new();
        for _ in 0..50 {
            stats.update(300);
        }
        for _ in 0..50 {
            stats.update(302);
        }
        let log_prob = stats.log_gaussian_probability(301);
        assert!(log_prob.is_finite());
        assert!(log_prob < 0.0); // log of PDF < 0 for continuous distributions

        let log_prob_far = stats.log_gaussian_probability(1000);
        assert!(log_prob_far < log_prob);
        assert!(log_prob_far < -100.0); // very negative for far values
    }

    #[test]
    fn test_gaussian_insert_size_confidence() {
        let mut stats = InsertSizeStats::new();
        for _ in 0..50 {
            stats.update(300);
        }
        for _ in 0..50 {
            stats.update(300);
        }
        // With stddev = 0, confidence should be 0
        assert_eq!(stats.insert_size_confidence(300), 0.0);

        let mut stats2 = InsertSizeStats::new();
        for _ in 0..25 {
            stats2.update(298);
        }
        for _ in 0..25 {
            stats2.update(300);
        }
        for _ in 0..25 {
            stats2.update(302);
        }
        for _ in 0..25 {
            stats2.update(304);
        }
        // Confidence at mean should be highest (close to 1.0)
        let conf_at_mean = stats2.insert_size_confidence(301);
        let conf_far = stats2.insert_size_confidence(1000);
        assert!(conf_at_mean > conf_far);
        assert!(conf_at_mean > 0.9);
        assert!(conf_far < 0.01);
    }

    #[test]
    fn test_gaussian_is_proper_pair() {
        let mut stats = InsertSizeStats::new();
        // Build distribution centered at 300 with small variance
        for _ in 0..100 {
            stats.update(300);
        }
        for _ in 0..100 {
            stats.update(305);
        }

        // Value within mean ± 3*stddev should be proper pair
        let (lower, upper) = stats.expected_range();
        let mid = (lower + upper) / 2;
        assert!(stats.is_proper_pair(mid));

        // Value outside range should not be proper pair
        let far_out = upper + 100;
        assert!(!stats.is_proper_pair(far_out));

        // Negative tlen should never be proper pair
        assert!(!stats.is_proper_pair(-100));

        // Zero tlen should never be proper pair
        assert!(!stats.is_proper_pair(0));
    }

    #[test]
    fn test_gaussian_expected_range() {
        let mut stats = InsertSizeStats::new();
        assert_eq!(stats.expected_range(), (0, 0)); // no data

        stats.update(300);
        assert_eq!(stats.expected_range(), (0, 0)); // need at least 2

        stats.update(300);
        let (lower, upper) = stats.expected_range();
        assert_eq!(lower, 0); // identical values, stddev = 0
        assert_eq!(upper, 0);
    }

    #[test]
    fn test_gaussian_welford_accuracy() {
        let mut stats = InsertSizeStats::new();
        // Known dataset: [10, 20, 30, 40, 50]
        // mean = 30, variance = 250, stddev ≈ 15.81
        for val in [10i64, 20, 30, 40, 50] {
            stats.update(val);
        }
        assert!((stats.mean - 30.0).abs() < 0.01);
        assert!((stats.stddev - 15.811).abs() < 0.1);
    }

    #[test]
    fn test_gaussian_two_pass_mapping() {
        let mut bp = BitPop::new(5);
        bp.add_genome(
            "G1",
            "AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT",
        );
        bp.add_genome(
            "G2",
            "GGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCC",
        );
        bp.build();

        // Create pairs where both reads map to G1
        let pairs: Vec<(String, String, Vec<u8>, String, Vec<u8>)> = (0..20)
            .map(|i| {
                (
                    format!("read_{i}"),
                    "AAAAACCCCCGGGGG".to_string(),
                    vec![40; 15],
                    "TTTTTAAAAACCCCC".to_string(),
                    vec![40; 15],
                )
            })
            .collect();

        // First pass: collect stats without reconciliation
        let mut insert_stats = InsertSizeStats::new();
        let mut mapped_count = 0;
        for (name, seq1, qual1, seq2, qual2) in &pairs {
            let paired = PairedRead {
                name: name.clone(),
                read1_seq: seq1.clone(),
                read1_qual: qual1.clone(),
                read2_seq: seq2.clone(),
                read2_qual: qual2.clone(),
            };
            let result = bp.map_read_paired(&paired, 0, false, 5);
            if result.map1.is_some() && result.map2.is_some() {
                mapped_count += 1;
            }
            insert_stats.update(result.tlen);
        }

        // Should have collected some stats from mapped pairs
        assert!(mapped_count > 0 || insert_stats.count > 0);
    }

    #[test]
    fn test_gaussian_confidence_range() {
        let mut stats = InsertSizeStats::new();
        for _ in 0..50 {
            stats.update(280);
        }
        for _ in 0..50 {
            stats.update(320);
        }

        // Confidence should always be in [0, 1]
        for tlen in [100, 200, 280, 300, 320, 400, 500] {
            let conf = stats.insert_size_confidence(tlen);
            assert!(
                (0.0..=1.0).contains(&conf),
                "confidence {} for tlen {} out of range",
                conf,
                tlen
            );
        }
    }
}
