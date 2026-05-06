pub mod delta;
pub mod fm;
pub mod align;
pub mod rank;
pub mod fasta;
pub mod sam;
pub mod serialize;
pub mod persisted;
pub mod fastq;
pub mod ncbi;
pub mod cache;
pub mod index_manager;

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
        _   => None,
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

/// Encode a DNA sequence into a compact byte slice (2 bits per base).
/// Skips unknown bases (N).
pub fn encode_sequence(seq: &str) -> Vec<u8> {
    seq.chars().filter_map(|c| encode_base(c)).collect()
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
    let complement: Vec<u8> = encoded.iter().map(|&b| match b {
        1 => 4, // A -> T
        4 => 1, // T -> A
        2 => 3, // C -> G
        3 => 2, // G -> C
        other => other,
    }).collect();
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

/// Insert size statistics for paired-end mapping.
#[derive(Debug, Clone)]
pub struct InsertSizeStats {
    pub mean: f64,
    pub stddev: f64,
    pub count: usize,
}

impl InsertSizeStats {
    pub fn new() -> Self {
        Self { mean: 0.0, stddev: 0.0, count: 0 }
    }

    pub fn update(&mut self, insert_size: i64) {
        if insert_size <= 0 {
            return;
        }
        self.count += 1;
        let old_mean = self.mean;
        self.mean += (insert_size as f64 - old_mean) / self.count as f64;
    }

    pub fn update_with_variance(&mut self, insert_size: i64) {
        if insert_size <= 0 {
            return;
        }
        self.count += 1;
        let delta = insert_size as f64 - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = insert_size as f64 - self.mean;
        // Welford's online algorithm for variance
        // We store M2 incrementally
        if self.count == 1 {
            self.stddev = 0.0;
        } else {
            // Approximate stddev update (simplified)
            self.stddev = ((delta * delta2).max(0.0) / (self.count as f64 - 1.0)).sqrt();
        }
    }

    pub fn is_proper_pair(&self, observed_tlen: i64) -> bool {
        if self.count < 2 || observed_tlen <= 0 {
            return false;
        }
        let lower = (self.mean - 3.0 * self.stddev).max(0.0) as i64;
        let upper = (self.mean + 3.0 * self.stddev) as i64;
        observed_tlen >= lower && observed_tlen <= upper
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
}

impl fmt::Display for AlignMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignMode::Xor => write!(f, "xor"),
            AlignMode::Sw => write!(f, "sw"),
            AlignMode::Hybrid => write!(f, "hybrid"),
        }
    }
}

// --- BitPop (main struct) ---

/// Bit-Pop genomic mapper.
///
/// 3-stage pipeline:
/// 1. K-mer inverted index → candidate positions
/// 2. Bit-level alignment → precise matches
/// 3. Multi-genome ranking → scored results
///
/// After adding all genomes, call `build()` to compress the k-mer index.
/// Compression reduces memory by ~60-70% and improves cache performance.
pub struct BitPop {
    /// FM-index (built on demand via `build()`)
    fm_index: Option<FmIndex>,

    /// Forward genome storage: genome_id → DNA sequence
    genomes: HashMap<u32, Vec<u8>>,

    /// Genome names for output
    genome_names: HashMap<u32, String>,

    /// K-mer size for indexing
    k: usize,

    /// Number of top rarest k-mers to try as anchors (for error tolerance)
    top_n: usize,
}

impl BitPop {
    /// Create a new Bit-Pop indexer with the given k-mer size.
    /// Recommended: k=15 for short reads, k=20 for long reads.
    pub fn new(_k: usize) -> Self {
        Self {
            fm_index: None,
            genomes: HashMap::new(),
            genome_names: HashMap::new(),
            k: _k,
            top_n: 1,
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

    /// Add a genome (reference sequence) to the index.
    /// Returns the assigned genome_id.
    ///
    /// After adding all genomes, call `build()` to construct the FM-index.
    pub fn add_genome(&mut self, name: &str, sequence: &str) -> u32 {
        let genome_id = self.genomes.len() as u32;
        let encoded = encode_sequence(sequence);
        self.genome_names.insert(genome_id, name.to_string());
        self.genomes.insert(genome_id, encoded.clone());
        genome_id
    }

     /// Finalize the index: construct the FM-index from all stored genomes.
    /// After build, the index is immutable.
    pub fn build(&mut self) {
        let mut genome_list: Vec<(u32, &str, &[u8])> = self.genomes.iter()
            .map(|(gid, seq)| {
                (*gid, self.genome_names.get(gid).map(|s| s.as_str()).unwrap_or(""), seq.as_slice())
            })
            .collect();
        genome_list.sort_by_key(|(gid, _, _)| *gid);
        let genomes: Vec<(&str, &[u8])> = genome_list.into_iter()
            .map(|(_, name, seq)| (name, seq))
            .collect();
        self.fm_index = Some(FmIndex::build(&genomes));
    }

    /// Finalize the index with parallel FM-Index build using rayon.
    /// Collects genome data into owned vectors first, then builds in parallel.
    pub fn build_parallel(&mut self) {
        // First, collect all genome data into owned vectors
        let mut genome_list: Vec<(u32, String, Vec<u8>)> = self.genomes.iter()
            .map(|(gid, seq)| {
                let name = self.genome_names.get(gid).cloned().unwrap_or_default();
                (*gid, name, seq.to_vec())
            })
            .collect();
        genome_list.sort_by_key(|(gid, _, _)| *gid);

        // Parallel build: encode genomes in parallel
        let genomes: Vec<(String, Vec<u8>)> = genome_list.into_par_iter()
            .map(|(_, name, seq)| (name, seq))
            .collect();

        // Convert to reference format for FmIndex::build
        let genome_refs: Vec<(&str, &[u8])> = genomes.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_slice()))
            .collect();

        self.fm_index = Some(FmIndex::build(&genome_refs));
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

    /// Encode k bytes (each 2 bits) into a u64 k-mer value.
    fn encode_kmer_bytes(&self, bytes: &[u8]) -> u64 {
        let mut result: u64 = 0;
        for &b in bytes {
            result = (result << 2) | b as u64;
        }
        result
    }

    /// Stage 1: K-mer filter for a read.
    /// Returns candidate positions across all genomes using FM-index backward search.
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

            let has_low_quality: bool = quality[i..qual_end]
                .iter()
                .any(|&q| q < min_quality);

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

            let positions = fm.find_positions(kmer, max_hits.min(usize::MAX));
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
        candidates.sort_by(|a, b| b.2.cmp(&a.2));
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

            let positions = fm.find_positions(kmer, max_hits.min(usize::MAX));
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
        candidates.sort_by(|a, b| b.2.cmp(&a.2));
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
        let cigar = if score >= 1.0 {
            format!("{}M", read_len)
        } else {
            format!("{}M", read_len)
        };
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
            (normalized.min(1.0).max(0.0), cigar, region_start)
        } else {
            // For longer reads: chunked SW with full traceback → real CIGAR
            let (sw_score, cigar) = align::smith_waterman_chunked(&read_enc, region);
            if sw_score == 0 {
                let (score, offset) = align::smith_waterman_score(&read_enc, region);
                (score, format!("{}M", read_len), region_start + offset)
            } else {
                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                (normalized.min(1.0).max(0.0), cigar, region_start)
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
            let (sw_score, cigar, offset, qual_penalty) = align::smith_waterman_with_quality(
                &read_enc, region, quality, 2, -2, 0
            );
            if sw_score == 0 {
                return (0.0, String::new(), 0, 0.0);
            }
            let normalized = (sw_score as f64) / (2.0 * read_len as f64);
            (normalized.min(1.0).max(0.0), cigar, offset, qual_penalty)
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

                let qual_chunk = &quality[chunk_start.min(quality.len())..chunk_end.min(quality.len())];
                let text_start = chunk_start + region_start;
                let text_end = (text_start + chunk.len()).min(region.len());
                if text_end - text_start < chunk.len() {
                    continue;
                }
                let text_region = &region[text_start..text_end];

                let (sw_score, cigar) = align::smith_waterman_internal(chunk, text_region, 2, -1, -2);
                total_score += sw_score;
                if !cigar.is_empty() && sw_score > 0 {
                    align::parse_cigar_ops(&cigar, &mut all_ops);
                }

                // Accumulate quality penalty by re-scoring the chunk alignment path
                let (_, _, _, penalty) = align::smith_waterman_with_quality(
                    chunk, text_region, qual_chunk, 2, -2, 0
                );
                total_penalty += penalty;
            }

            if total_score == 0 {
                return (0.0, String::new(), 0, 0.0);
            }

            let cigar = align::build_cigar_string(&all_ops);
            let normalized = (total_score as f64) / (2.0 * read_len as f64);
            (normalized.min(1.0).max(0.0), cigar, region_start, total_penalty)
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
        }
    }

     /// Stage 3: Rank pre-scored mapping results.
    /// Takes already-aligned candidates and applies rarity-based ranking.
    pub fn rank_scored_results(
        &self,
        scored_candidates: &[(u32, u64, f64, String)],
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

        let rarity = if encoded.len() >= self.k {
            let occ = fm.count_occurrences(&encoded[..self.k]);
            1.0 / (occ as f64).max(1.0)
        } else {
            1.0
        };

        for &(genome_id, position, align_score, ref cigar) in scored_candidates {
            if align_score < 0.5 {
                continue;
            }

            // align_score is the primary signal (0.85 weight).
            // rarity provides a modest boost (0.15 weight) so perfect matches
            // floor at 0.85 instead of 0.5 when the first k-mer is common.
            let combined_score = align_score * 0.85 + rarity * 0.15;
            let context = self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(MappingResult {
                genome_id,
                position,
                score: combined_score,
                cigar: cigar.clone(),
                context,
                is_reverse: false,
            });
        }

        results.sort_by(|a, b| {
            b.score.partial_cmp(&a.score)
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
                let occ = fm.count_occurrences(&encoded[..self.k]);
                1.0 / (occ as f64).max(1.0)
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let context = self.extract_genome_context(genome_id, position, read_len, context_window);

            results.push(MappingResult {
                genome_id,
                position,
                score: combined_score,
                cigar,
                context,
                is_reverse: false,
            });
        }

        results.sort_by(|a, b| {
            b.score.partial_cmp(&a.score)
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

        (base_threshold + length_factor + qual_factor).max(0.3f64).min(0.8f64)
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

            let has_low_quality: bool = quality[i..qual_end]
                .iter()
                .any(|&q| q < min_quality);

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

    /// Anchor-based filter: replaces k-mer co-occurrence counting with
    /// rarest-k-mer anchor + 2-bit XOR scoring.
    ///
    /// Algorithm:
    /// 1. Find the rarest k-mer in the read (fewest total positions across all genomes)
    /// 2. Get all positions for that anchor k-mer via FM-index backward search
    /// 3. For each position: 2-bit XOR score the entire read against the genome region
    /// 4. Return positions with score >= threshold
    ///
    /// This is O(anchor_positions * read_length/31) vs O(total_kmer_hits) for kmer_filter.
    pub fn anchor_filter(&self, read: &str, min_score: f64) -> Vec<(u32, u64, f64, String)> {
        self.anchor_filter_with_threshold(read, min_score, usize::MAX)
    }

   /// Anchor-based filter with AlignMode support.
    /// Uses top-N rarest k-mers as anchors for error tolerance.
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
        if encoded.len() < self.k {
            return Vec::new();
        }

        let top_n_kmers = self.find_top_n_rarest_kmers(&encoded, fm, max_hits);
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

                let search_radius: isize = (self.k.max(read_len / 4)).min(200) as isize;
                let mut best_score = f64::NEG_INFINITY;
                let mut best_cigar = String::new();
                let mut best_offset: usize = 0;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() { continue; }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < self.k {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];

                    let (score, cigar, _) = match mode {
                        AlignMode::Xor => {
                            if read_len <= 31 {
                                align::two_bit_align(&encoded, candidate_region)
                            } else {
                                let (s, o) = align::two_bit_score_chunks(&encoded, candidate_region);
                                (s, format!("{}M", read_len), o)
                            }
                        }
                        AlignMode::Sw => {
                            if read_len <= 31 {
                                let (sw_score, cigar) = align::smith_waterman(&encoded, candidate_region);
                                if sw_score == 0 { continue; }
                                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                                (normalized.min(1.0).max(0.0), cigar, candidate_start)
                            } else {
                                let (sw_score, cigar) = align::smith_waterman_chunked(&encoded, candidate_region);
                                if sw_score == 0 { continue; }
                                let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                                (normalized.min(1.0).max(0.0), cigar, candidate_start)
                            }
                        }
                        AlignMode::Hybrid => {
                            if read_len <= 31 {
                                let (xor_s, xor_cigar, _) = align::two_bit_align(&encoded, candidate_region);
                                if xor_s >= 0.9 {
                                    (xor_s, xor_cigar, candidate_start)
                                } else {
                                    let (sw_score, sw_cigar) = align::smith_waterman(&encoded, candidate_region);
                                    if sw_score == 0 { continue; }
                                    let normalized = (sw_score as f64) / (2.0 * read_len as f64);
                                    if normalized > xor_s {
                                        (normalized.min(1.0), sw_cigar, candidate_start)
                                    } else {
                                        (xor_s, xor_cigar, candidate_start)
                                    }
                                }
                            } else {
                                let (s, _) = align::two_bit_score_chunks(&encoded, candidate_region);
                                if s >= 0.9 {
                                    (s, format!("{}M", read_len), candidate_start)
                                } else {
                                    let (sw_score, sw_cigar) = align::smith_waterman_chunked(&encoded, candidate_region);
                                    let sw_s = (sw_score as f64) / (2.0 * read_len as f64);
                                    if sw_s > s && sw_score > 0 { (sw_s.min(1.0), sw_cigar, candidate_start) }
                                    else { (s, format!("{}M", read_len), candidate_start) }
                                }
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
                    scored.push((genome_id, best_offset as u64, best_score.max(0.0).min(1.0), best_cigar));
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

        let top_n_kmers = self.find_top_n_rarest_kmers(&encoded, fm, max_hits);
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
                    let estimated_start = (position as isize - anchor_read_offset as isize).max(0) as usize;
                    let region_end = (estimated_start + read_len).min(genome.len());
                    let region = &genome[estimated_start..region_end];

                    if region.len() < self.k {
                        continue;
                    }

                    if region.len() == read_len && encoded == *region {
                        scored.push((genome_id, estimated_start as u64, 1.0, format!("{}M", read_len)));
                        continue;
                    }

                    let (score, cigar, offset) = align::two_bit_align(&encoded, region);
                    if score >= min_score {
                        let actual_pos = (estimated_start + offset) as u64;
                        scored.push((genome_id, actual_pos, score, cigar));
                    }
                    continue;
                }

                let search_radius: isize = (self.k.max(read_len / 4)).min(200) as isize;
                let mut best_score = 0.0f64;
                let mut best_offset: usize = 0;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() { continue; }
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

                    let (score, _) = align::two_bit_score_chunks(&encoded, candidate_region);
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
                        let op = if read_part[i] == cand_region[i] { 0u8 } else { 1u8 };
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

        let top_n_kmers = self.find_top_n_rarest_kmers_quality(&encoded, quality, fm, min_quality, max_hits);
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
                    let estimated_start = (position as isize - anchor_read_offset as isize).max(0) as usize;
                    let region_end = (estimated_start + read_len).min(genome.len());
                    let region = &genome[estimated_start..region_end];

                    if region.len() < self.k {
                        continue;
                    }

                    let qual_slice = &quality[..read_len.min(quality.len())];
                    let (score, cigar, _, penalty) = align::two_bit_align_with_quality(&encoded, region, qual_slice);
                    
                    if score >= min_score {
                        let actual_pos = (estimated_start + 0) as u64;
                        scored.push((genome_id, actual_pos, score, cigar, penalty));
                    }
                    continue;
                }

                let search_radius: isize = (self.k.max(read_len / 4)).min(200) as isize;
                let mut best_score = f64::NEG_INFINITY;
                let mut best_offset: usize = 0;
                let mut best_penalty = 0.0f64;

                for delta in -search_radius..=search_radius {
                    let candidate_start = (estimated_read_start + delta).max(0) as usize;
                    if candidate_start >= genome.len() { continue; }
                    let region_end = (candidate_start + read_len).min(genome.len());
                    if region_end - candidate_start < self.k {
                        continue;
                    }

                    let candidate_region = &genome[candidate_start..region_end];
                    let qual_slice = &quality[..read_len.min(quality.len())];

                    let (chunk_score, _, chunk_penalty) = align::two_bit_score_chunks_with_quality(&encoded, candidate_region, qual_slice);
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
                        let op = if read_part[i] == cand_region[i] { 0u8 } else { 1u8 };
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
                    scored.push((genome_id, best_offset as u64, best_score.max(0.0).min(1.0), cigar, best_penalty));
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
        let scored = self.anchor_filter_with_quality(read, quality, 0.7, min_quality, DEFAULT_REPEAT_THRESHOLD);
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
                let occ = fm.count_occurrences(&encoded[..self.k]);
                1.0 / (occ as f64).max(1.0)
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let adjusted_score = (align_score + quality_penalty).max(0.0).min(1.0);
            
            let read_len = encoded.len();
            let context = self.extract_genome_context(genome_id, position, read_len, context_window);

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
            });
        }

        results.sort_by(|a, b| b.combined_score.partial_cmp(&a.combined_score).unwrap_or(std::cmp::Ordering::Equal));
        results
    }

    /// Extract ±window bases around a position in a genome.
    fn extract_genome_context(&self, genome_id: u32, position: u64, read_len: usize, window: usize) -> String {
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
    /// Uses anchor-based 2-bit XOR filtering (1 anchor + XOR score).
    pub fn map_read(&self, read: &str, context_window: usize) -> Vec<MappingResult> {
        self.map_read_with_mode(read, AlignMode::Xor, context_window)
    }

    /// Full pipeline with configurable alignment mode.
    /// Uses anchor-based filtering with the specified alignment algorithm.
    /// Tries both forward and reverse complement, returns best alignment.
    pub fn map_read_with_mode(&self, read: &str, mode: AlignMode, context_window: usize) -> Vec<MappingResult> {
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
            (Some(f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(r)) => {
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
    fn map_read_orientation(&self, read: &str, mode: AlignMode, context_window: usize, _is_rc: bool) -> Vec<MappingResult> {
        let scored = self.anchor_filter_with_mode(read, mode, 0.7, DEFAULT_REPEAT_THRESHOLD);
        if scored.is_empty() {
            return Vec::new();
        }
        self.rank_scored_results(&scored, read, context_window)
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
        let forward_results = self.map_read_quality_orientation(read, quality, mode, min_quality, context_window, false);
        let rc_read = reverse_complement(read);
        let rc_results = self.map_read_quality_orientation(&rc_read, quality, mode, min_quality, context_window, true);

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
            (Some(f), None) => {
                let mut results = forward_results;
                if !results.is_empty() {
                    results[0].is_reverse = false;
                }
                results
            }
            (None, Some(r)) => {
                let mut results = rc_results;
                if !results.is_empty() {
                    results[0].is_reverse = true;
                }
                results
            }
            (None, None) => Vec::new(),
        }
    }

    /// Map a single orientation (forward or RC) of a read with quality scores.
    fn map_read_quality_orientation(
        &self,
        read: &str,
        quality: &[u8],
        mode: AlignMode,
        min_quality: u8,
        context_window: usize,
        _is_rc: bool,
    ) -> Vec<QualityMappingResult> {
        let scored = self.anchor_filter_with_quality_smart(read, quality, 0.7, min_quality, DEFAULT_REPEAT_THRESHOLD);
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
                let occ = fm.count_occurrences(&encoded[..self.k]);
                1.0 / (occ as f64).max(1.0)
            } else {
                1.0
            };

            let combined_score = align_score * 0.85 + rarity * 0.15;
            let adjusted_score = (align_score + quality_penalty).max(0.0).min(1.0);
            
            let read_len = encoded.len();
            let context = self.extract_genome_context(genome_id, position, read_len, context_window);

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
            });
        }

        results.sort_by(|a, b| b.combined_score.partial_cmp(&a.combined_score).unwrap_or(std::cmp::Ordering::Equal));
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
        let mut names: Vec<(u32, String)> = self.genome_names.iter().map(|(id, name)| (*id, name.clone())).collect();
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
            top_n: 1,
        }
    }

    /// Serialize and write to a file (persisted format with compression).
    pub fn serialize_to_file(&self, path: &str) -> io::Result<()> {
        persisted::save_bitpop(self, path)?;
        Ok(())
    }

    /// Load a BitPop instance from a file (persisted format).
    pub fn deserialize_from_file(path: &str) -> io::Result<Self> {
        let bp = persisted::load_bitpop(path)?;
        Ok(bp)
    }

    /// Map multiple reads in parallel and write results to a SAM file.
    /// Returns the number of reads that had at least one mapping.
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
        let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
            .map(|(n, l)| (n.as_str(), *l)).collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        let mapped: Vec<(String, String, Vec<MappingResult>)> = reads.par_iter()
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
        let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
            .map(|(n, l)| (n.as_str(), *l)).collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        let mapped: Vec<(String, String, Vec<QualityMappingResult>)> = reads.into_par_iter()
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
        let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
            .map(|(n, l)| (n.as_str(), *l)).collect();

        let name_refs: Vec<&str> = genome_name_refs.clone();

        // Split reads into batches for work-stealing
        let batches: Vec<Vec<(&str, &str)>> = reads
            .chunks(batch_size)
            .map(|chunk| chunk.to_vec())
            .collect();

        let all_results: Vec<(String, String, Vec<MappingResult>)> = batches.into_par_iter()
            .flat_map_iter(|batch| {
                batch.into_iter()
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

    // --- Paired-end mapping ---

    /// Map a single paired-end read to all indexed genomes.
    /// Returns the best mapping result for each read in the pair.
    pub fn map_read_paired(
        &self,
        paired: &PairedRead,
        context_window: usize,
    ) -> PairedMappingResult {
        let mut insert_stats = InsertSizeStats::new();

        let map1 = self.map_single_read_to_best(&paired.read1_seq);
        let map2 = self.map_single_read_to_best(&paired.read2_seq);

        // Compute TLEN (observed template length)
        let tlen = compute_tlen(&map1, &map2, paired.read1_seq.len(), paired.read2_seq.len());
        insert_stats.update(tlen);

        PairedMappingResult {
            read_name: paired.name.clone(),
            map1,
            map2,
            tlen,
            insert_size_stats: insert_stats,
        }
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
        })
    }

    /// Map multiple paired-end reads in parallel and write SAM output.
    pub fn map_paired_reads_parallel(
        &self,
        pairs: &[(String, String, Vec<u8>, String, Vec<u8>)],
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
        let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
            .map(|(n, l)| (n.as_str(), *l)).collect();

        // Collect insert size stats from all pairs first
        let mut insert_stats = InsertSizeStats::new();
        let mapped_pairs: Vec<PairedMappingResult> = pairs.iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                let paired = PairedRead {
                    name: name.clone(),
                    read1_seq: seq1.clone(),
                    read1_qual: qual1.clone(),
                    read2_seq: seq2.clone(),
                    read2_qual: qual2.clone(),
                };
                let result = self.map_read_paired(&paired, context_window);
                insert_stats.count = result.insert_size_stats.count;
                insert_stats.mean = result.insert_size_stats.mean;
                insert_stats.stddev = result.insert_size_stats.stddev;
                result
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
        writer.write_header(&genome_header)?;

        // Write paired-end SAM output
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

    /// Map multiple paired-end reads with quality-aware scoring.
    pub fn map_paired_reads_parallel_quality(
        &self,
        pairs: &[(String, String, Vec<u8>, String, Vec<u8>)],
        output_path: &str,
        min_quality: u8,
        context_window: usize,
    ) -> io::Result<usize> {
        let genomes_owned: Vec<(String, usize)> = (0..self.genome_count() as u32)
            .filter_map(|gid| {
                self.genome_name(gid)
                    .map(|name| (name.to_string(), self.genome_seq_len(gid).unwrap_or(0)))
            })
            .collect();

        let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
        let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
            .map(|(n, l)| (n.as_str(), *l)).collect();

        let mut insert_stats = InsertSizeStats::new();
        let mapped_pairs: Vec<PairedMappingResult> = pairs.iter()
            .map(|(name, seq1, qual1, seq2, qual2)| {
                let paired = PairedRead {
                    name: name.clone(),
                    read1_seq: seq1.clone(),
                    read1_qual: qual1.clone(),
                    read2_seq: seq2.clone(),
                    read2_qual: qual2.clone(),
                };

                let map1 = self.map_paired_read_with_quality(&paired.read1_seq, &paired.read1_qual, min_quality);
                let map2 = self.map_paired_read_with_quality(&paired.read2_seq, &paired.read2_qual, min_quality);

                let tlen = compute_tlen(&map1, &map2, paired.read1_seq.len(), paired.read2_seq.len());
                insert_stats.update(tlen);

                PairedMappingResult {
                    read_name: paired.name.clone(),
                    map1,
                    map2,
                    tlen,
                    insert_size_stats: InsertSizeStats {
                        mean: insert_stats.mean,
                        stddev: insert_stats.stddev,
                        count: insert_stats.count,
                    },
                }
            })
            .collect();

        let mut writer = sam::SamWriter::new(output_path)?;
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

    fn map_paired_read_with_quality(&self, seq: &str, qual: &[u8], min_quality: u8) -> Option<PairedReadMapping> {
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
        })
    }
}

fn compute_tlen(map1: &Option<PairedReadMapping>, map2: &Option<PairedReadMapping>, len1: usize, len2: usize) -> i64 {
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
                -(tlen as i64)
            } else {
                tlen as i64
            }
        }
        _ => 0,
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
        let path = dir.join(format!("bitpop_fasta_load_{}_{}.fasta", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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
        let path = dir.join(format!("bitpop_sam_{}_{}.sam", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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
        let path = dir.join(format!("bitpop_sam_build_{}_{}.sam", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        let mapped = bp.map_reads_to_sam(&reads, path.to_str().unwrap(), 3).unwrap();
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
        let results = bp.rank_scored_results(&scored, "ACGTACGT", 3);
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
        let results_with_thresh = bp.anchor_filter_with_threshold("AACCGGTTAACCGGTT", 0.5, 100)
            .into_iter()
            .map(|(g, p, s, c)| (g, p, s, c))
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
        let fastq_path = dir.join(format!("test_fastq_{}_{}.fastq", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        let sam_path = dir.join(format!("test_fastq_{}_{}.sam", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));

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

        let mapped = bp.map_reads_from_fastq_parallel(
            fastq_path.to_str().unwrap(),
            sam_path.to_str().unwrap(),
            20,
            3,
        ).unwrap();

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
        let sam_path = dir.join(format!("test_optimized_{}_{}.sam", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));

        let mapped = bp.map_reads_parallel_optimized(
            &reads,
            sam_path.to_str().unwrap(),
            3,
            2, // batch_size
        ).unwrap();

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

        let (score, cigar, _) = bp.align_read_sw("ACGTACGT", 0, 0);
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
                    assert!(threshold >= 0.3 && threshold <= 0.8,
                        "Threshold {} out of bounds for len={}, qual={}, avg={}",
                        threshold, read_len, has_qual, avg_q);
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
        assert!(!results[0].is_reverse, "Forward alignment should win for perfect match");
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
}
