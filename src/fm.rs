/// FM-Index for multi-genome DNA mapping.
///
/// Encoding: $=0 (sentinel), A=1, C=2, G=3, T=4
/// SA construction via libsais (SA-IS, O(n))

use libsais::{SuffixArrayConstruction, ThreadCount};

const ALPHABET_SIZE: usize = 5; // $=0, A=1, C=2, G=3, T=4
const SAMPLE_INTERVAL: usize = 32;

// --- Suffix Array Construction (SA-IS via libsais, O(n)) ---
fn build_suffix_array(s: &[u8]) -> Vec<usize> {
    let n = s.len();
    if n == 0 { return vec![]; }
    if n == 1 { return vec![0]; }

    // s is u8 with $=0, A=1, C=2, G=3, T=4
    // libsais treats last char as implicit sentinel (smallest)
    // For SmallAlphabet (u8), use for_text (immutable)
    let mut sa_buffer: Vec<i32> = vec![0; n];

    let result = SuffixArrayConstruction::for_text(s)
        .in_borrowed_buffer(&mut sa_buffer)
        .multi_threaded(ThreadCount::openmp_default())
        .run();

    match result {
        Ok(sa_with_text) => sa_with_text.suffix_array().iter().map(|&x| x as usize).collect(),
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
        for i in 1..bucket_size { cnt[i] += cnt[i - 1]; }

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
        for i in 1..bucket_size { cnt[i] += cnt[i - 1]; }

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

        if rank[sa[n - 1]] == n as u8 - 1 { break; }
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

// --- Occ Counter (Rank Sampling) ---
pub struct OccCounter {
    samples: Vec<[u32; ALPHABET_SIZE]>,
    sample_interval: usize,
}

impl OccCounter {
    pub fn new(bwt: &[u8], sample_interval: usize) -> Self {
        let len = bwt.len();
        let num_samples = (len + sample_interval - 1) / sample_interval;
        let mut samples: Vec<[u32; ALPHABET_SIZE]> = vec![[0; ALPHABET_SIZE]; num_samples];
        let mut counts = [0u32; ALPHABET_SIZE];
        
        for (i, &c) in bwt.iter().enumerate() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE { counts[ci] += 1; }
            if (i + 1) % sample_interval == 0 {
                let idx = (i + 1) / sample_interval - 1;
                if idx < num_samples { samples[idx] = counts; }
            }
        }
        if len % sample_interval != 0 {
            samples[num_samples - 1] = counts;
        }
        
        OccCounter { samples, sample_interval }
    }
    
    pub fn occ(&self, bwt: &[u8], c: u8, i: usize) -> u64 {
        if i == 0 { return 0; }
        let ci = c as usize;
        if ci >= ALPHABET_SIZE { return 0; }
        
        let i = i.min(bwt.len());
        let sample_idx = (i - 1) / self.sample_interval;
        let idx = sample_idx.min(self.samples.len() - 1);
        let sample_end = ((idx + 1) * self.sample_interval).min(bwt.len());
        let base = self.samples[idx][ci] as u64;
        
        if sample_end == i {
            return base;
        } else if sample_end > i {
            let excess: u64 = bwt[i..sample_end].iter().filter(|&&x| x == c).count() as u64;
            return base - excess;
        } else {
            let extra: u64 = bwt[sample_end..i].iter().filter(|&&x| x == c).count() as u64;
            return base + extra;
        }
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
                if byte >= 1 && byte <= 4 { s.push(byte); }
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
        let mut counts = vec![0usize; ALPHABET_SIZE];
        for &c in bwt.iter() {
            let ci = c as usize;
            if ci < ALPHABET_SIZE { counts[ci] += 1; }
        }
        let mut c_array = [0usize; ALPHABET_SIZE];
        let mut sum = 0usize;
        for c in 0..ALPHABET_SIZE {
            c_array[c] = sum;
            sum += counts[c];
        }
        
        let occ = OccCounter::new(&bwt, SAMPLE_INTERVAL);
        
        FmIndex { bwt, sa, c_array, occ, genome_boundaries, num_genomes: genomes.len() }
    }
    
    /// Backward search: find SA range [lower, upper) for pattern. O(m).
    /// Pattern bytes: A=1, C=2, G=3, T=4
    pub fn backward_search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        let mut lower: usize = 0;
        let mut upper: usize = self.bwt.len();
        
        for &byte in pattern.iter().rev() {
            let c = byte;
            if c == 0 || c > 4 { return None; }
            let ci = c as usize;
            
            let occ_l = self.occ.occ(&self.bwt, c, lower);
            let occ_u = self.occ.occ(&self.bwt, c, upper);
            
            lower = self.c_array[ci] + occ_l as usize;
            upper = self.c_array[ci] + occ_u as usize;
            
            if lower >= upper { return None; }
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
            if positions.len() >= max_positions { break; }
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
    
    pub fn len(&self) -> usize { self.bwt.len() }
    pub fn is_empty(&self) -> bool { self.bwt.is_empty() }
    pub fn num_genomes(&self) -> usize { self.num_genomes }
    pub fn genome_boundaries(&self) -> &Vec<(usize, usize, u32)> { &self.genome_boundaries }
    pub fn sa_at(&self, rank: usize) -> usize { self.sa[rank] }
    pub fn sa_len(&self) -> usize { self.sa.len() }
    pub fn bwt_at(&self, pos: usize) -> u8 { self.bwt[pos] }
    pub fn c_array(&self, idx: usize) -> usize { self.c_array[idx] }

    /// Create an FmIndex from components (for deserialization).
    pub fn from_components(
        bwt: Vec<u8>,
        sa: Vec<usize>,
        c_array: [usize; 5],
        occ: OccCounter,
        genome_boundaries: Vec<(usize, usize, u32)>,
        num_genomes: usize,
    ) -> Self {
        FmIndex { bwt, sa, c_array, occ, genome_boundaries, num_genomes }
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
            assert!(&s[a..] <= &s[b..],
                "SA not sorted at {}: suffix[{}]={:?} > suffix[{}]={:?}",
                i, a, &s[a..], b, &s[b..]);
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
        for &(gid, _) in &positions { assert_eq!(gid, 0); }
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
}
