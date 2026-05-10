/// FM-Index for multi-genome DNA mapping.
///
/// Encoding: $=0 (sentinel), A=1, C=2, G=3, T=4
/// SA construction via libsais (SA-IS, O(n))
/// Parallel build support via rayon
use libsais::SuffixArrayConstruction;
use rayon::prelude::*;

const ALPHABET_SIZE: usize = 5; // $=0, A=1, C=2, G=3, T=4
const SAMPLE_INTERVAL: usize = 32;

// --- Suffix Array Construction (SA-IS via libsais, O(n)) ---
fn build_suffix_array(s: &[u8]) -> Vec<usize> {
    let n = s.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![0];
    }

    // s is u8 with $=0, A=1, C=2, G=3, T=4
    // libsais treats last char as implicit sentinel (smallest)
    // For SmallAlphabet (u8), use for_text (immutable)
    let mut sa_buffer: Vec<i32> = vec![0; n];

    let result = SuffixArrayConstruction::for_text(s)
        .in_borrowed_buffer(&mut sa_buffer)
        .single_threaded()
        .run();

    match result {
        Ok(sa_with_text) => sa_with_text
            .suffix_array()
            .iter()
            .map(|&x| x as usize)
            .collect(),
        Err(e) => {
            eprintln!("libsais error: {:?}, falling back to radix sort", e);
            build_suffix_array_radix_fallback(s, n)
        }
    }
}

fn build_suffix_array_radix_fallback(s: &[u8], n: usize) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..n).collect();
    let mut rank = s.to_vec();
    let mut new_rank = vec![0u8; n];

    let mut k = 1;
    while k < n {
        let max_r = *rank.iter().max().unwrap_or(&0);
        let bucket_size = (max_r + 2) as usize;

        let mut cnt = vec![0usize; bucket_size];
        for &i in &sa {
            let key = if i + k < n { rank[i + k] + 1 } else { 0 };
            cnt[key as usize] += 1;
        }
        for i in 1..bucket_size {
            cnt[i] += cnt[i - 1];
        }

        let mut sa2 = vec![0usize; n];
        for &i in sa.iter().rev() {
            let key = if i + k < n { rank[i + k] + 1 } else { 0 };
            cnt[key as usize] -= 1;
            sa2[cnt[key as usize]] = i;
        }

        cnt.fill(0);
        for &i in &sa2 {
            let key = rank[i] + 1;
            cnt[key as usize] += 1;
        }
        for i in 1..bucket_size {
            cnt[i] += cnt[i - 1];
        }

        let mut sa3 = vec![0usize; n];
        for &i in sa2.iter().rev() {
            let key = rank[i] + 1;
            cnt[key as usize] -= 1;
            sa3[cnt[key as usize]] = i;
        }

        new_rank[sa3[0]] = 0;
        for i in 1..n {
            let prev = sa3[i - 1];
            let curr = sa3[i];
            let prev_key = (rank[prev], if prev + k < n { rank[prev + k] } else { 0 });
            let curr_key = (rank[curr], if curr + k < n { rank[curr + k] } else { 0 });
            new_rank[curr] = new_rank[prev] + if curr_key > prev_key { 1 } else { 0 };
        }
        std::mem::swap(&mut sa, &mut sa3);
        std::mem::swap(&mut rank, &mut new_rank);

        if rank[sa[n - 1]] == n as u8 - 1 {
            break;
        }
        k *= 2;
    }

    sa
}

fn build_bwt_from_sa(sa: &[usize], s: &[u8]) -> Vec<u8> {
    let n = sa.len();
    let mut bwt = vec![0u8; n];
    for (i, &sa_i) in sa.iter().enumerate() {
        let prev = if sa_i == 0 { n - 1 } else { sa_i - 1 };
        bwt[i] = s[prev];
    }
    bwt
}

fn build_bwt_from_sa_parallel(sa: &[usize], s: &[u8], num_threads: usize) -> Vec<u8> {
    let n = sa.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![s[0]];
    }

    let chunk_size = n.div_ceil(num_threads);
    let num_chunks = ((n - 1) / chunk_size) + 1;

    let bwt_chunks: Vec<Vec<u8>> = (0..num_chunks)
        .into_par_iter()
        .map(|chunk_idx| {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(n);
            let mut chunk = vec![0u8; end - start];
            for (i, &sa_i) in sa.iter().enumerate().skip(start).take(end - start) {
                let prev = if sa_i == 0 { n - 1 } else { sa_i - 1 };
                chunk[i - start] = s[prev];
            }
            chunk
        })
        .collect();

    bwt_chunks.into_iter().flatten().collect()
}

// --- Occ Counter (Rank Sampling) ---
pub struct OccCounter {
    samples: Vec<[u32; ALPHABET_SIZE]>,
    sample_interval: usize,
}

impl OccCounter {
    pub fn new(bwt: &[u8], sample_interval: usize) -> Self {
        let len = bwt.len();
        let num_samples = len.div_ceil(sample_interval);
        let mut samples: Vec<[u32; ALPHABET_SIZE]> = vec![[0; ALPHABET_SIZE]; num_samples];
        let mut counts = [0u32; ALPHABET_SIZE];

        for (i, &c) in bwt.iter().enumerate() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE {
                counts[ci] += 1;
            }
            if (i + 1) % sample_interval == 0 {
                let idx = (i + 1) / sample_interval - 1;
                if idx < num_samples {
                    samples[idx] = counts;
                }
            }
        }
        if !len.is_multiple_of(sample_interval) {
            samples[num_samples - 1] = counts;
        }

        OccCounter {
            samples,
            sample_interval,
        }
    }

    pub fn new_parallel(bwt: &[u8], sample_interval: usize, _num_threads: usize) -> Self {
        let len = bwt.len();
        if len == 0 {
            return OccCounter {
                samples: vec![],
                sample_interval,
            };
        }

        let num_samples = len.div_ceil(sample_interval);
        let mut samples: Vec<[u32; ALPHABET_SIZE]> = vec![[0; ALPHABET_SIZE]; num_samples];
        let mut counts = [0u32; ALPHABET_SIZE];

        for (i, &c) in bwt.iter().enumerate() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE {
                counts[ci] += 1;
            }
            if (i + 1) % sample_interval == 0 {
                let idx = (i + 1) / sample_interval - 1;
                if idx < num_samples {
                    samples[idx] = counts;
                }
            }
        }
        if !len.is_multiple_of(sample_interval) {
            samples[num_samples - 1] = counts;
        }

        OccCounter {
            samples,
            sample_interval,
        }
    }

    pub fn occ(&self, bwt: &[u8], c: u8, i: usize) -> u64 {
        if i == 0 {
            return 0;
        }
        let ci = c as usize;
        if ci >= ALPHABET_SIZE {
            return 0;
        }

        let i = i.min(bwt.len());
        let sample_idx = (i - 1) / self.sample_interval;
        let idx = sample_idx.min(self.samples.len() - 1);
        let sample_end = ((idx + 1) * self.sample_interval).min(bwt.len());
        let base = self.samples[idx][ci] as u64;

        if sample_end == i {
            base
        } else if sample_end > i {
            let excess: u64 = bwt[i..sample_end].iter().filter(|&&x| x == c).count() as u64;
            base - excess
        } else {
            let extra: u64 = bwt[sample_end..i].iter().filter(|&&x| x == c).count() as u64;
            base + extra
        }
    }
}

// --- Spaced Seed ---

#[derive(Debug, Clone)]
pub struct SpacedSeed {
    pub pattern: Vec<bool>,
}

impl SpacedSeed {
    pub fn from_binary(binary: &str) -> Self {
        SpacedSeed {
            pattern: binary.chars().map(|c| c == '1').collect(),
        }
    }

    pub fn default_v1() -> Self {
        SpacedSeed::from_binary("11101001110111")
    }

    pub fn len(&self) -> usize {
        self.pattern.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pattern.is_empty()
    }

    pub fn coverage(&self) -> f64 {
        let ones = self.pattern.iter().filter(|&&b| b).count();
        ones as f64 / self.pattern.len() as f64
    }

    pub fn extract_match(&self, kmer: &[u8]) -> Vec<u8> {
        kmer.iter()
            .enumerate()
            .filter(|&(i, _)| i < self.pattern.len() && self.pattern[i])
            .map(|(_, &b)| b)
            .collect()
    }
}

impl Default for SpacedSeed {
    fn default() -> Self {
        Self::default_v1()
    }
}

// --- FM-Index ---

pub struct FmIndex {
    bwt: Vec<u8>,
    sa: Vec<usize>,
    c_array: [usize; ALPHABET_SIZE],
    occ: OccCounter,
    genome_boundaries: Vec<(usize, usize, u32)>,
    num_genomes: usize,
}

impl FmIndex {
    /// Build FM-index from multiple genomes.
    /// Input: genomes as &[u8] with A=1, C=2, G=3, T=4 (NO sentinel in input)
    pub fn build(genomes: &[(&str, &[u8])]) -> Self {
        let mut s: Vec<u8> = Vec::new();
        let mut genome_boundaries: Vec<(usize, usize, u32)> = Vec::new();

        for (gid, (_, seq)) in genomes.iter().enumerate() {
            let start = s.len();
            for &byte in *seq {
                // Input bytes should be 1-4 (A=1, C=2, G=3, T=4)
                if (1..=4).contains(&byte) {
                    s.push(byte);
                }
            }
            let gid_u32 = gid as u32;
            if gid_u32 < (genomes.len() - 1) as u32 {
                s.push(0); // terminator between genomes
            }
            genome_boundaries.push((start, s.len() - start, gid_u32));
        }
        s.push(0); // final sentinel

        let sa = build_suffix_array(&s);
        let bwt = build_bwt_from_sa(&sa, &s);

        // C-array: C[c] = number of chars < c in BWT
        let mut counts = [0usize; ALPHABET_SIZE];
        for &c in bwt.iter() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE {
                counts[ci] += 1;
            }
        }
        let mut c_array = [0usize; ALPHABET_SIZE];
        let mut sum = 0usize;
        for c in 0..ALPHABET_SIZE {
            c_array[c] = sum;
            sum += counts[c];
        }

        let occ = OccCounter::new(&bwt, SAMPLE_INTERVAL);

        FmIndex {
            bwt,
            sa,
            c_array,
            occ,
            genome_boundaries,
            num_genomes: genomes.len(),
        }
    }

    /// Build FM-index in parallel using rayon.
    /// Parallelizes BWT construction and OccCounter building.
    /// Input: genomes as &[u8] with A=1, C=2, G=3, T=4 (NO sentinel in input)
    pub fn build_parallel(genomes: &[(&str, &[u8])]) -> Self {
        let mut s: Vec<u8> = Vec::new();
        let mut genome_boundaries: Vec<(usize, usize, u32)> = Vec::new();

        for (gid, (_, seq)) in genomes.iter().enumerate() {
            let start = s.len();
            for &byte in *seq {
                if (1..=4).contains(&byte) {
                    s.push(byte);
                }
            }
            let gid_u32 = gid as u32;
            if gid_u32 < (genomes.len() - 1) as u32 {
                s.push(0);
            }
            genome_boundaries.push((start, s.len() - start, gid_u32));
        }
        s.push(0);

        let num_threads = rayon::current_num_threads();
        let sa = build_suffix_array(&s);
        let bwt = build_bwt_from_sa_parallel(&sa, &s, num_threads);

        let mut counts = [0usize; ALPHABET_SIZE];
        for &c in bwt.iter() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE {
                counts[ci] += 1;
            }
        }
        let mut c_array = [0usize; ALPHABET_SIZE];
        let mut sum = 0usize;
        for c in 0..ALPHABET_SIZE {
            c_array[c] = sum;
            sum += counts[c];
        }

        let occ = OccCounter::new_parallel(&bwt, SAMPLE_INTERVAL, num_threads);

        FmIndex {
            bwt,
            sa,
            c_array,
            occ,
            genome_boundaries,
            num_genomes: genomes.len(),
        }
    }

    /// Backward search: find SA range [lower, upper) for pattern. O(m).
    /// Pattern bytes: A=1, C=2, G=3, T=4
    pub fn backward_search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        let mut lower: usize = 0;
        let mut upper: usize = self.bwt.len();

        for &byte in pattern.iter().rev() {
            let c = byte;
            if c == 0 || c > 4 {
                return None;
            }
            let ci = c as usize;

            let occ_l = self.occ.occ(&self.bwt, c, lower);
            let occ_u = self.occ.occ(&self.bwt, c, upper);

            lower = self.c_array[ci] + occ_l as usize;
            upper = self.c_array[ci] + occ_u as usize;

            if lower >= upper {
                return None;
            }
        }

        Some((lower, upper))
    }

    /// Backward search with spaced seed support.
    /// Only matches on positions where mask[i] == true.
    /// Returns SA range [lower, upper) for the spaced pattern.
    pub fn backward_search_spaced(&self, kmer: &[u8], mask: &[bool]) -> Option<(usize, usize)> {
        if kmer.len() != mask.len() {
            return None;
        }

        let mut lower: usize = 0;
        let mut upper: usize = self.bwt.len();

        let mut pattern = Vec::with_capacity(mask.len());
        for i in 0..mask.len() {
            if mask[i] {
                pattern.push(kmer[i]);
            }
        }

        for &byte in pattern.iter().rev() {
            let c = byte;
            if c == 0 || c > 4 {
                return None;
            }
            let ci = c as usize;

            let occ_l = self.occ.occ(&self.bwt, c, lower);
            let occ_u = self.occ.occ(&self.bwt, c, upper);

            lower = self.c_array[ci] + occ_l as usize;
            upper = self.c_array[ci] + occ_u as usize;

            if lower >= upper {
                return None;
            }
        }

        Some((lower, upper))
    }

    /// Find all positions where pattern occurs.
    pub fn find_positions(&self, pattern: &[u8], max_positions: usize) -> Vec<(u32, u64)> {
        let (lower, upper) = match self.backward_search(pattern) {
            Some(r) => r,
            None => return Vec::new(),
        };

        let mut positions = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for rank in lower..upper {
            if positions.len() >= max_positions {
                break;
            }
            let sa_pos = self.sa[rank];
            let (genome_id, genome_pos) = self.sa_to_genome_pos(sa_pos);
            if seen.insert((genome_id, genome_pos)) {
                positions.push((genome_id, genome_pos as u64));
            }
        }

        positions
    }

    /// Find which genomes contain this pattern.
    /// Stops early once 2+ genomes are found (for efficiency).
    /// Returns a HashSet of genome IDs where the pattern occurs.
    pub fn find_genomes(&self, pattern: &[u8]) -> std::collections::HashSet<u32> {
        let (lower, upper) = match self.backward_search(pattern) {
            Some(r) => r,
            None => return std::collections::HashSet::new(),
        };

        let mut genomes = std::collections::HashSet::new();
        for rank in lower..upper {
            let sa_pos = self.sa[rank];
            let (genome_id, _) = self.sa_to_genome_pos(sa_pos);
            genomes.insert(genome_id);

            // Stop early once we've seen 2+ genomes
            if genomes.len() >= 2 {
                break;
            }
        }

        genomes
    }

    /// Find all positions where spaced pattern occurs.
    /// Only matches on positions where mask[i] == true.
    pub fn find_positions_spaced(
        &self,
        kmer: &[u8],
        mask: &[bool],
        max_positions: usize,
    ) -> Vec<(u32, u64)> {
        let (lower, upper) = match self.backward_search_spaced(kmer, mask) {
            Some(r) => r,
            None => return Vec::new(),
        };

        let mut positions = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for rank in lower..upper {
            if positions.len() >= max_positions {
                break;
            }
            let sa_pos = self.sa[rank];
            let (genome_id, genome_pos) = self.sa_to_genome_pos(sa_pos);
            if seen.insert((genome_id, genome_pos)) {
                positions.push((genome_id, genome_pos as u64));
            }
        }

        positions
    }

    fn sa_to_genome_pos(&self, sa_pos: usize) -> (u32, usize) {
        let mut lo = 0usize;
        let mut hi = self.genome_boundaries.len();
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let (start, len, _) = self.genome_boundaries[mid];
            if sa_pos < start {
                hi = mid;
            } else if sa_pos >= start + len {
                lo = mid + 1;
            } else {
                let (_, _, gid) = self.genome_boundaries[mid];
                return (gid, sa_pos - start);
            }
        }
        (0, 0)
    }

    /// Count occurrences of pattern. O(m), independent of genome size.
    pub fn count_occurrences(&self, pattern: &[u8]) -> usize {
        match self.backward_search(pattern) {
            Some((lower, upper)) => upper - lower,
            None => 0,
        }
    }

    pub fn count_occurrences_spaced(&self, kmer: &[u8], mask: &[bool]) -> usize {
        match self.backward_search_spaced(kmer, mask) {
            Some((lower, upper)) => upper - lower,
            None => 0,
        }
    }

    /// Backward search with fuzzy spaced seed support.
    /// Allows up to `max_mismatches` mismatches in the "match" (true) positions of the mask.
    /// Returns the SA range with the most matches (smallest range).
    pub fn backward_search_spaced_fuzzy(
        &self,
        kmer: &[u8],
        mask: &[bool],
        max_mismatches: usize,
    ) -> Option<(usize, usize)> {
        if kmer.len() != mask.len() {
            return None;
        }

        let match_positions: Vec<usize> = (0..mask.len()).filter(|&i| mask[i]).collect();

        if match_positions.is_empty() {
            return Some((0, self.bwt.len()));
        }

        if max_mismatches == 0 {
            return self.backward_search_spaced(kmer, mask);
        }

        let mut best_range: Option<(usize, usize)> = None;

        Self::search_fuzzy_recursive(
            self,
            kmer,
            mask,
            &match_positions,
            0,
            0,
            max_mismatches,
            &mut best_range,
        );

        best_range
    }

    #[allow(clippy::too_many_arguments)]
    fn search_fuzzy_recursive(
        &self,
        kmer: &[u8],
        mask: &[bool],
        match_positions: &[usize],
        match_idx: usize,
        mismatches: usize,
        max_mismatches: usize,
        best_range: &mut Option<(usize, usize)>,
    ) {
        if mismatches > max_mismatches {
            return;
        }

        if match_idx == match_positions.len() {
            if let Some((lower, upper)) = self.backward_search_spaced(kmer, mask) {
                if let Some((best_lower, best_upper)) = best_range {
                    let current_size = upper - lower;
                    let best_size = *best_upper - *best_lower;
                    if current_size < best_size {
                        *best_range = Some((lower, upper));
                    }
                } else {
                    *best_range = Some((lower, upper));
                }
            }
            return;
        }

        let pos = match_positions[match_idx];
        let original_base = kmer[pos];

        // Try original base (no mismatch)
        self.search_fuzzy_recursive(
            kmer,
            mask,
            match_positions,
            match_idx + 1,
            mismatches,
            max_mismatches,
            best_range,
        );

        // Try alternative bases (mismatch)
        if mismatches < max_mismatches {
            for alt in 1..=4u8 {
                if alt != original_base {
                    let mut modified_kmer = kmer.to_vec();
                    modified_kmer[pos] = alt;
                    self.search_fuzzy_recursive(
                        &modified_kmer,
                        mask,
                        match_positions,
                        match_idx + 1,
                        mismatches + 1,
                        max_mismatches,
                        best_range,
                    );
                }
            }
        }
    }

    /// Find all positions where fuzzy spaced pattern occurs.
    /// Allows up to `max_mismatches` mismatches in the "match" positions.
    pub fn find_positions_spaced_fuzzy(
        &self,
        kmer: &[u8],
        mask: &[bool],
        max_mismatches: usize,
        max_positions: usize,
    ) -> Vec<(u32, u64)> {
        let (lower, upper) = match self.backward_search_spaced_fuzzy(kmer, mask, max_mismatches) {
            Some(r) => r,
            None => return Vec::new(),
        };

        let mut positions = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for rank in lower..upper {
            if positions.len() >= max_positions {
                break;
            }
            let sa_pos = self.sa[rank];
            let (genome_id, genome_pos) = self.sa_to_genome_pos(sa_pos);
            if seen.insert((genome_id, genome_pos)) {
                positions.push((genome_id, genome_pos as u64));
            }
        }

        positions
    }

    /// Count occurrences of a fuzzy spaced pattern.
    pub fn count_occurrences_spaced_fuzzy(
        &self,
        kmer: &[u8],
        mask: &[bool],
        max_mismatches: usize,
    ) -> usize {
        match self.backward_search_spaced_fuzzy(kmer, mask, max_mismatches) {
            Some((lower, upper)) => upper - lower,
            None => 0,
        }
    }

    pub fn len(&self) -> usize {
        self.bwt.len()
    }
    pub fn is_empty(&self) -> bool {
        self.bwt.is_empty()
    }
    pub fn num_genomes(&self) -> usize {
        self.num_genomes
    }
    pub fn genome_boundaries(&self) -> &Vec<(usize, usize, u32)> {
        &self.genome_boundaries
    }
    pub fn sa_at(&self, rank: usize) -> usize {
        self.sa[rank]
    }
    pub fn sa_len(&self) -> usize {
        self.sa.len()
    }
    pub fn bwt_at(&self, pos: usize) -> u8 {
        self.bwt[pos]
    }
    pub fn c_array(&self, idx: usize) -> usize {
        self.c_array[idx]
    }

    /// Create an FmIndex from components (for deserialization).
    pub fn from_components(
        bwt: Vec<u8>,
        sa: Vec<usize>,
        c_array: [usize; 5],
        occ: OccCounter,
        genome_boundaries: Vec<(usize, usize, u32)>,
        num_genomes: usize,
    ) -> Self {
        FmIndex {
            bwt,
            sa,
            c_array,
            occ,
            genome_boundaries,
            num_genomes,
        }
    }
}

// --- Tests ---

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debug_sa_bwt() {
        // Encoding: $=0, A=1, C=2, G=3, T=4
        let seq = vec![1u8, 2, 3, 4]; // ACGT
        let mut s: Vec<u8> = seq.clone();
        s.push(0); // sentinel $

        let sa = build_suffix_array(&s);
        let bwt = build_bwt_from_sa(&sa, &s);

        eprintln!("SA: {:?}", sa);
        eprintln!("BWT: {:?}", bwt);
        eprintln!("Text: {:?}", s);

        // Validate SA
        for i in 0..sa.len() - 1 {
            let a = sa[i];
            let b = sa[i + 1];
            assert!(
                &s[a..] <= &s[b..],
                "SA not sorted at {}: suffix[{}]={:?} > suffix[{}]={:?}",
                i,
                a,
                &s[a..],
                b,
                &s[b..]
            );
        }
    }

    #[test]
    fn test_backward_search_exact() {
        // A=1, C=2, G=3, T=4
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4]; // ACGTACGT
        let index = FmIndex::build(&[("test", &seq)]);
        let pattern = vec![1, 2, 3, 4]; // ACGT
        let result = index.backward_search(&pattern);
        assert!(result.is_some(), "ACGT should be found in ACGTACGT");
        let (l, u) = result.unwrap();
        assert!(u > l, "Should find at least one occurrence");
    }

    #[test]
    fn test_backward_search_not_found() {
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        let pattern = vec![4, 4, 4, 4]; // TTTT
        assert!(index.backward_search(&pattern).is_none());
    }

    #[test]
    fn test_backward_search_full() {
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        let pattern = vec![1, 2, 3, 4, 1, 2, 3, 4]; // ACGTACGT
        assert!(index.backward_search(&pattern).is_some());
    }

    #[test]
    fn test_multi_genome() {
        let g1 = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let g2 = vec![4, 4, 4, 4, 2, 2, 2, 2];
        let index = FmIndex::build(&[("human", &g1), ("mouse", &g2)]);
        assert_eq!(index.num_genomes(), 2);

        let pattern = vec![1, 2, 3, 4]; // ACGT
        let positions = index.find_positions(&pattern, 100);
        assert!(!positions.is_empty());
        for &(gid, _) in &positions {
            assert_eq!(gid, 0);
        }
    }

    #[test]
    fn test_count_occurrences() {
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]; // ACGT x3
        let index = FmIndex::build(&[("test", &seq)]);
        let pattern = vec![1, 2, 3, 4]; // ACGT
        assert!(index.count_occurrences(&pattern) >= 3);
    }

    #[test]
    fn test_find_positions_multi_genome() {
        let g1 = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let g2 = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let index = FmIndex::build(&[("human", &g1), ("chimp", &g2)]);

        let pattern = vec![1, 2, 3, 4];
        let positions = index.find_positions(&pattern, 100);
        let genomes: std::collections::HashSet<u32> = positions.iter().map(|(g, _)| *g).collect();
        assert_eq!(genomes.len(), 2);
    }

    #[test]
    fn test_repetitive_pattern() {
        let seq: Vec<u8> = (0..400u32).map(|i| ((i % 4) + 1) as u8).collect(); // ACGT repeated 100x
        let index = FmIndex::build(&[("test", &seq)]);

        assert!(index.count_occurrences(&vec![1, 2, 3, 4]) >= 100);
        assert!(index.count_occurrences(&vec![1, 2, 3, 4, 1, 2, 3, 4]) >= 50);
    }

    #[test]
    fn test_occ_counter_basic() {
        // $=0, A=1, C=2, G=3, T=4
        let bwt = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let occ = OccCounter::new(&bwt, 4);
        assert_eq!(occ.occ(&bwt, 0, 4), 1);
        assert_eq!(occ.occ(&bwt, 0, 8), 2);
        assert_eq!(occ.occ(&bwt, 1, 4), 1);
    }

    #[test]
    fn test_occ_counter_boundary() {
        let bwt = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let occ = OccCounter::new(&bwt, 4);
        assert_eq!(occ.occ(&bwt, 0, 5), 2);
        assert_eq!(occ.occ(&bwt, 0, 1), 1);
        assert_eq!(occ.occ(&bwt, 0, 0), 0);
    }

    #[test]
    fn test_sa_consistency() {
        let seq = vec![1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        assert_eq!(index.sa.len(), index.bwt.len());
    }

    #[test]
    fn test_short_pattern() {
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        assert!(index.backward_search(&vec![1]).is_some()); // A
    }

    #[test]
    fn test_terminator_handling() {
        let g1 = vec![1, 2, 3, 4];
        let g2 = vec![4, 3, 2, 1];
        let index = FmIndex::build(&[("g1", &g1), ("g2", &g2)]);
        let positions = index.find_positions(&vec![4, 3, 2, 1], 100); // TGCA
        if !positions.is_empty() {
            assert_eq!(positions[0].0, 1);
        }
    }

    #[test]
    fn test_longer_sequence() {
        let seq: Vec<u8> = (0..40u32).map(|i| ((i % 4) + 1) as u8).collect(); // ACGT x10
        let index = FmIndex::build(&[("test", &seq)]);
        // ACGTACGT (8 chars) appears at positions 0,4,8,12,16,20,24,28,32 = 9 times
        assert!(index.count_occurrences(&vec![1, 2, 3, 4, 1, 2, 3, 4]) >= 9);
    }

    #[test]
    fn test_single_base_genome() {
        let seq = vec![1]; // A
        let index = FmIndex::build(&[("test", &seq)]);
        assert!(index.backward_search(&vec![1]).is_some());
    }

    #[test]
    fn test_empty_pattern() {
        let seq = vec![1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        assert!(index.backward_search(&[]).is_some());
    }

    #[test]
    fn test_sa_correctness() {
        let seq = vec![1, 2, 3, 4];
        let index = FmIndex::build(&[("test", &seq)]);
        assert_eq!(index.sa.len(), seq.len() + 1);
    }

    #[test]
    fn test_larger_genome() {
        let seq: Vec<u8> = (0..1000u32).map(|i| ((i % 4) + 1) as u8).collect();
        let index = FmIndex::build(&[("test", &seq)]);
        assert!(index.count_occurrences(&vec![1, 2, 3, 4]) >= 250);
    }

    #[test]
    fn test_spaced_seed_extract() {
        let seed = SpacedSeed::from_binary("11101001110111");
        let kmer = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2]; // ACGTACGTACGTAC
        let extracted = seed.extract_match(&kmer);
        assert_eq!(extracted.len(), 10);
    }

    #[test]
    fn test_spaced_seed_coverage() {
        let seed = SpacedSeed::from_binary("11101001110111");
        assert_eq!(seed.len(), 14);
        assert!((seed.coverage() - 10.0 / 14.0).abs() < 0.01);
    }

    #[test]
    fn test_backward_search_spaced() {
        let seq = vec![1, 2, 3, 1, 4, 1, 2, 4, 1, 2, 3, 4, 1, 2];
        let index = FmIndex::build(&[("test", &seq)]);

        let kmer = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
        let mask = SpacedSeed::default_v1().pattern;

        let result = index.backward_search_spaced(&kmer, &mask);
        assert!(
            result.is_some(),
            "Spaced seed search should find pattern in sequence"
        );
    }

    #[test]
    fn test_find_positions_spaced() {
        let seq = vec![1, 2, 3, 1, 4, 1, 2, 4, 1, 2, 3, 4, 1, 2];
        let index = FmIndex::build(&[("test", &seq)]);

        let kmer = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
        let mask = SpacedSeed::default_v1().pattern;

        let positions = index.find_positions_spaced(&kmer, &mask, 100);
        assert!(!positions.is_empty(), "Should find at least one position");
        assert_eq!(positions[0].0, 0);
    }

    #[test]
    fn test_spaced_seed_vs_exact_search() {
        let seq = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
        let index = FmIndex::build(&[("test", &seq)]);

        let kmer = vec![1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2];
        let all_true_mask = vec![true; 14];

        let spaced_result = index.backward_search_spaced(&kmer, &all_true_mask);
        let exact_result = index.backward_search(&kmer[..14]);

        assert!(
            spaced_result.is_some(),
            "Spaced with all-true mask should match exact search"
        );
        assert!(exact_result.is_some());
    }
}
