/// SNP detection and mapping for strain resolution.
///
/// Detects single nucleotide polymorphisms (SNPs) by analyzing mismatch patterns
/// across multiple reads. SNPs are positions where >min_support reads have
/// the same mismatch against the reference genome.
///
/// # Example
///
/// ```
/// use bit_pop::snp::SnpDetector;
///
/// let mut detector = SnpDetector::new();
/// detector.add_mismatch(0, 100, 1, 3); // genome 0, position 100, A->G
/// detector.add_mismatch(0, 100, 1, 3);
/// detector.add_mismatch(0, 100, 1, 3);
/// detector.build_snp_map(3);
/// assert!(detector.is_known_snp(0, 100).is_some());
/// ```

use std::collections::HashMap;

/// A single nucleotide polymorphism detected from mismatch patterns.
#[derive(Debug, Clone)]
pub struct Snp {
    /// Genome ID where this SNP occurs.
    pub genome_id: u32,
    /// Position in the genome (0-indexed).
    pub position: u32,
    /// Reference base (from genome).
    pub ref_base: u8,
    /// Alternative base (from read).
    pub alt_base: u8,
    /// Number of reads supporting this SNP.
    pub support_count: u32,
}

/// SNP detector that analyzes mismatch patterns to identify strain-specific SNPs.
///
/// Collects mismatches from multiple reads, groups them by (genome_id, position),
/// and identifies SNPs where >min_support reads have the same base change.
pub struct SnpDetector {
    /// Mismatch counts: (genome_id, position, ref_base, alt_base) → count
    mismatch_counts: HashMap<(u32, u32, u8, u8), u32>,
    /// Detected SNPs
    snps: Vec<Snp>,
    /// Whether SNP map has been built
    built: bool,
}

impl SnpDetector {
    /// Create a new empty SNP detector.
    pub fn new() -> Self {
        Self {
            mismatch_counts: HashMap::new(),
            snps: Vec::new(),
            built: false,
        }
    }

    /// Add a mismatch observation.
    ///
    /// # Arguments
    /// * `genome_id` — Genome ID where mismatch occurred
    /// * `position` — Position in the genome (0-indexed)
    /// * `read_base` — Base from the read (alternative)
    /// * `genome_base` — Base from the genome (reference)
    pub fn add_mismatch(&mut self, genome_id: u32, position: u32, read_base: u8, genome_base: u8) {
        let key = (genome_id, position, genome_base, read_base);
        *self.mismatch_counts.entry(key).or_default() += 1;
    }

    /// Build the SNP map from collected mismatches.
    ///
    /// Only mismatches with >= min_support observations are considered SNPs.
    ///
    /// # Arguments
    /// * `min_support` — Minimum number of reads supporting a SNP (default: 3)
    pub fn build_snp_map(&mut self, min_support: u32) {
        self.snps.clear();

        for ((genome_id, position, ref_base, alt_base), count) in &self.mismatch_counts {
            if *count >= min_support {
                self.snps.push(Snp {
                    genome_id: *genome_id,
                    position: *position,
                    ref_base: *ref_base,
                    alt_base: *alt_base,
                    support_count: *count,
                });
            }
        }

        self.snps.sort_by_key(|snp| (snp.genome_id, snp.position));
        self.built = true;
    }

    /// Check if a position is a known SNP.
    ///
    /// # Arguments
    /// * `genome_id` — Genome ID to check
    /// * `position` — Position to check
    ///
    /// # Returns
    /// Some(&Snp) if this position is a known SNP, None otherwise
    pub fn is_known_snp(&self, genome_id: u32, position: u32) -> Option<&Snp> {
        if !self.built {
            return None;
        }

        self.snps.iter().find(|snp| {
            snp.genome_id == genome_id && snp.position == position
        })
    }

    /// Check if a mismatch is on a known SNP position.
    ///
    /// This is the key function for strain resolution:
    /// - If mismatch is on known SNP position → likely from different strain
    /// - If mismatch is random → likely sequencing error
    ///
    /// # Arguments
    /// * `genome_id` — Genome ID where mismatch occurred
    /// * `position` — Position in the genome
    /// * `read_base` — Base from the read
    /// * `genome_base` — Base from the genome
    ///
    /// # Returns
    /// true if this mismatch is on a known SNP position
    pub fn is_snp_mismatch(&self, genome_id: u32, position: u32, read_base: u8, genome_base: u8) -> bool {
        match self.is_known_snp(genome_id, position) {
            Some(snp) => snp.alt_base == read_base && snp.ref_base == genome_base,
            None => false,
        }
    }

    /// Get all SNPs for a specific genome.
    pub fn get_genome_snps(&self, genome_id: u32) -> &[Snp] {
        if !self.built {
            return &[];
        }

        let start = self.snps.partition_point(|snp| snp.genome_id < genome_id);
        let end = self.snps.partition_point(|snp| snp.genome_id <= genome_id);

        &self.snps[start..end]
    }

    /// Get all detected SNPs.
    pub fn get_all_snps(&self) -> &[Snp] {
        &self.snps
    }

    /// Get the number of detected SNPs.
    pub fn snp_count(&self) -> usize {
        self.snps.len()
    }

    /// Check if the SNP map has been built.
    pub fn is_built(&self) -> bool {
        self.built
    }

    /// Export SNP map to JSON format.
    ///
    /// # Returns
    /// JSON string with SNP data
    pub fn to_json(&self) -> String {
        if self.snps.is_empty() {
            return "[]".to_string();
        }

        let mut json = String::from("[\n");
        for (i, snp) in self.snps.iter().enumerate() {
            json.push_str(&format!(
                "  {{\"genome_id\": {}, \"position\": {}, \"ref\": \"{}\", \"alt\": \"{}\", \"support\": {}}}",
                snp.genome_id,
                snp.position,
                base_to_char(snp.ref_base),
                base_to_char(snp.alt_base),
                snp.support_count,
            ));
            if i < self.snps.len() - 1 {
                json.push(',');
            }
            json.push('\n');
        }
        json.push(']');
        json
    }
}

/// Convert 2-bit encoded base to character.
fn base_to_char(base: u8) -> char {
    match base {
        1 => 'A',
        2 => 'C',
        3 => 'G',
        4 => 'T',
        _ => 'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snp_detection_basic() {
        let mut detector = SnpDetector::new();

        // Add 3 identical mismatches (A->G at position 100)
        detector.add_mismatch(0, 100, 3, 1); // G (read) vs A (genome)
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        detector.build_snp_map(3);

        assert_eq!(detector.snp_count(), 1);
        assert!(detector.is_known_snp(0, 100).is_some());
    }

    #[test]
    fn test_snp_detection_below_threshold() {
        let mut detector = SnpDetector::new();

        // Add only 2 mismatches (below threshold of 3)
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        detector.build_snp_map(3);

        assert_eq!(detector.snp_count(), 0);
        assert!(detector.is_known_snp(0, 100).is_none());
    }

    #[test]
    fn test_snp_detection_multiple_positions() {
        let mut detector = SnpDetector::new();

        // SNP at position 100
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        // SNP at position 200
        detector.add_mismatch(0, 200, 4, 2);
        detector.add_mismatch(0, 200, 4, 2);
        detector.add_mismatch(0, 200, 4, 2);

        detector.build_snp_map(3);

        assert_eq!(detector.snp_count(), 2);
        assert!(detector.is_known_snp(0, 100).is_some());
        assert!(detector.is_known_snp(0, 200).is_some());
        assert!(detector.is_known_snp(0, 150).is_none());
    }

    #[test]
    fn test_snp_mismatch_detection() {
        let mut detector = SnpDetector::new();

        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        detector.build_snp_map(3);

        assert!(detector.is_snp_mismatch(0, 100, 3, 1)); // Correct SNP
        assert!(!detector.is_snp_mismatch(0, 100, 4, 1)); // Wrong alt base
        assert!(!detector.is_snp_mismatch(0, 100, 3, 2)); // Wrong ref base
        assert!(!detector.is_snp_mismatch(0, 101, 3, 1)); // Wrong position
    }

    #[test]
    fn test_json_export() {
        let mut detector = SnpDetector::new();

        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        detector.build_snp_map(3);

        let json = detector.to_json();
        assert!(json.contains("\"genome_id\": 0"));
        assert!(json.contains("\"position\": 100"));
        assert!(json.contains("\"ref\": \"A\""));
        assert!(json.contains("\"alt\": \"G\""));
        assert!(json.contains("\"support\": 3"));
    }

    #[test]
    fn test_empty_detector() {
        let detector = SnpDetector::new();
        assert_eq!(detector.snp_count(), 0);
        assert!(!detector.is_built());
        assert!(detector.is_known_snp(0, 100).is_none());
        assert_eq!(detector.to_json(), "[]");
    }

    #[test]
    fn test_get_genome_snps() {
        let mut detector = SnpDetector::new();

        // SNPs for genome 0
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);
        detector.add_mismatch(0, 100, 3, 1);

        // SNPs for genome 1
        detector.add_mismatch(1, 200, 4, 2);
        detector.add_mismatch(1, 200, 4, 2);
        detector.add_mismatch(1, 200, 4, 2);

        detector.build_snp_map(3);

        let genome_0_snps = detector.get_genome_snps(0);
        let genome_1_snps = detector.get_genome_snps(1);
        let genome_2_snps = detector.get_genome_snps(2);

        assert_eq!(genome_0_snps.len(), 1);
        assert_eq!(genome_1_snps.len(), 1);
        assert_eq!(genome_2_snps.len(), 0);
    }
}
