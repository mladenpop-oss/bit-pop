/// Multi-genome ranking for mapping results.
/// Extends vibe-build's Rarity × Proximity to genomic context.
use crate::MappingResult;

/// Rank mapping candidates by combined genomic score.
///
/// Score = alignment_score × relevance
///   relevance = rarity_weight × proximity_weight × coverage_weight
///   rarity = 1 / genomes_containing_kmer (rarer = more specific)
///   proximity = kmer_count / read_length (more k-mers matched = tighter)
///   coverage = matched_length / read_length
pub fn rank_mappings(
    candidates: &[(u32, u64, usize)],
    align_scores: &[(u32, u64, f64)],
    genome_kmer_counts: &std::collections::HashMap<u64, usize>,
    read_length: usize,
) -> Vec<MappingResult> {
    if candidates.is_empty() {
        return Vec::new();
    }

    let mut results = Vec::new();

    for &(genome_id, position, kmer_count) in candidates {
        // Find alignment score for this candidate
        let align_score = align_scores
            .iter()
            .find(|&&(gid, pos, _)| gid == genome_id && pos == position)
            .map(|&(_, _, s)| s)
            .unwrap_or(0.0);

        if align_score < 0.3 {
            continue;
        }

        // Rarity: how many genomes share the top k-mer?
        // Higher genome_kmer_count = more common = lower rarity score
        let kmer_occ = genome_kmer_counts.values().sum::<usize>();
        let rarity = 1.0 / (kmer_occ.max(1) as f64);

        // Proximity: more k-mers matched relative to read length = better
        let proximity = (kmer_count as f64 / read_length as f64).min(1.0);

        // Combined relevance
        let relevance = (0.4 * rarity) + (0.6 * proximity);
        let final_score = align_score * (0.5 + 0.5 * relevance);

        results.push(MappingResult {
            genome_id,
            position,
            score: final_score,
            cigar: format!("{}M", read_length),
            context: String::new(),
            is_reverse: false,
        });
    }

    results.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results
}

/// Filter results: keep only those above a score threshold.
pub fn filter_by_score(results: Vec<MappingResult>, threshold: f64) -> Vec<MappingResult> {
    results
        .into_iter()
        .filter(|r| r.score >= threshold)
        .collect()
}

/// Deduplicate results: keep highest score per (genome_id, position window).
pub fn deduplicate(results: Vec<MappingResult>, window: u64) -> Vec<MappingResult> {
    let mut sorted = results;
    sorted.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut kept = Vec::new();
    let mut used: Vec<(u32, u64)> = Vec::new();

    for r in sorted {
        let dominated = used
            .iter()
            .any(|&(gid, pos)| gid == r.genome_id && r.position.abs_diff(pos) < window);
        if !dominated {
            used.push((r.genome_id, r.position));
            kept.push(r);
        }
    }

    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_rank_mappings_basic() {
        let candidates = vec![(0, 100, 5), (0, 200, 3), (1, 300, 4)];
        let scores = vec![(0, 100, 0.95), (0, 200, 0.80), (1, 300, 0.90)];
        let mut kmer_counts = HashMap::new();
        kmer_counts.insert(0x1A, 2);
        let results = rank_mappings(&candidates, &scores, &kmer_counts, 10);
        assert!(!results.is_empty());
        assert!(results[0].score >= results[1].score);
    }

    #[test]
    fn test_rank_mappings_empty() {
        let results = rank_mappings(&[], &[], &HashMap::new(), 10);
        assert!(results.is_empty());
    }

    #[test]
    fn test_rank_mappings_low_score_filtered() {
        let candidates = vec![(0, 100, 1)];
        let scores = vec![(0, 100, 0.1)];
        let results = rank_mappings(&candidates, &scores, &HashMap::new(), 10);
        assert!(results.is_empty());
    }

    #[test]
    fn test_filter_by_score() {
        let results = vec![
            MappingResult {
                genome_id: 0,
                position: 100,
                score: 0.9,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
            MappingResult {
                genome_id: 0,
                position: 200,
                score: 0.5,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
            MappingResult {
                genome_id: 0,
                position: 300,
                score: 0.3,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
        ];
        let filtered = filter_by_score(results, 0.6);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].score, 0.9);
    }

    #[test]
    fn test_deduplicate() {
        let results = vec![
            MappingResult {
                genome_id: 0,
                position: 100,
                score: 0.9,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
            MappingResult {
                genome_id: 0,
                position: 102,
                score: 0.85,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
            MappingResult {
                genome_id: 0,
                position: 500,
                score: 0.8,
                cigar: "10M".into(),
                context: String::new(),
                is_reverse: false,
            },
        ];
        let deduped = deduplicate(results, 10);
        assert_eq!(deduped.len(), 2);
        assert_eq!(deduped[0].position, 100);
        assert_eq!(deduped[1].position, 500);
    }
}
