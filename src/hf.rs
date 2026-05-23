use std::collections::HashMap;

/// Homopolymer fingerprint profile for a single genome.
/// Stores normalized frequency of homopolymer runs (base, length).
#[derive(Debug, Clone)]
pub struct HfProfile {
    /// Normalized run frequencies: (base, run_len) -> frequency
    /// base: 0=A, 1=C, 2=G, 3=T
    pub run_freq: HashMap<(u8, usize), f64>,
    /// Total number of runs counted
    pub total_runs: usize,
    /// Minimum run length used during computation
    pub min_run: usize,
}

impl HfProfile {
    /// Compute homopolymer fingerprint from a genome sequence.
    /// Scans the sequence for runs of identical bases >= min_run length,
    /// then normalizes counts to frequencies.
    pub fn compute(sequence: &[u8], min_run: usize) -> Self {
        if sequence.is_empty() || min_run < 2 {
            return Self {
                run_freq: HashMap::new(),
                total_runs: 0,
                min_run,
            };
        }

        let mut counts: HashMap<(u8, usize), usize> = HashMap::new();
        let mut total = 0usize;

        let mut i = 0;
        while i < sequence.len() {
            let base = sequence[i];
            if base > 4 {
                i += 1;
                continue;
            }

            let mut run_len = 1usize;
            while i + run_len < sequence.len() && sequence[i + run_len] == base {
                run_len += 1;
            }

            if run_len >= min_run {
                *counts.entry((base, run_len)).or_default() += 1;
                total += 1;
            }

            i += run_len;
        }

        if total == 0 {
            return Self {
                run_freq: HashMap::new(),
                total_runs: 0,
                min_run,
            };
        }

        let run_freq = counts
            .into_iter()
            .map(|(k, v)| (k, v as f64 / total as f64))
            .collect();

        Self {
            run_freq,
            total_runs: total,
            min_run,
        }
    }

    /// Compute homopolymer fingerprint for a read and return the profile.
    pub fn compute_read(sequence: &[u8], min_run: usize) -> HashMap<(u8, usize), f64> {
        if sequence.is_empty() || min_run < 2 {
            return HashMap::new();
        }

        let mut counts: HashMap<(u8, usize), usize> = HashMap::new();
        let mut total = 0usize;

        let mut i = 0;
        while i < sequence.len() {
            let base = sequence[i];
            if base > 4 {
                i += 1;
                continue;
            }

            let mut run_len = 1usize;
            while i + run_len < sequence.len() && sequence[i + run_len] == base {
                run_len += 1;
            }

            if run_len >= min_run {
                *counts.entry((base, run_len)).or_default() += 1;
                total += 1;
            }

            i += run_len;
        }

        if total == 0 {
            return HashMap::new();
        }

        counts
            .into_iter()
            .map(|(k, v)| (k, v as f64 / total as f64))
            .collect()
    }

    /// Compute similarity score between a read's homopolymer profile and this genome's profile.
    /// Returns value in [0.0, 1.0] where 1.0 = perfect match.
    /// Uses 1 - Jensen-Shannon divergence for comparability.
    pub fn similarity(&self, read_profile: &HashMap<(u8, usize), f64>) -> f64 {
        if self.run_freq.is_empty() || read_profile.is_empty() {
            return 0.5;
        }

        let mut all_keys: std::collections::BTreeSet<(u8, usize)> =
            std::collections::BTreeSet::new();
        for k in self.run_freq.keys() {
            all_keys.insert(*k);
        }
        for k in read_profile.keys() {
            all_keys.insert(*k);
        }

        let mut js_div = 0.0f64;

        for &key in &all_keys {
            let p = *self.run_freq.get(&key).unwrap_or(&0.0);
            let q = *read_profile.get(&key).unwrap_or(&0.0);
            let m = (p + q) / 2.0;

            if m > 0.0 {
                if p > 0.0 {
                    js_div += p * (p / m).ln();
                }
                if q > 0.0 {
                    js_div += q * (q / m).ln();
                }
            }
        }

        js_div /= 2.0;
        let js_distance = js_div.sqrt();
        (1.0 - js_distance).clamp(0.0, 1.0)
    }
}

/// Serialize HfProfile to bytes for persistence.
pub fn serialize_hf_profile(profile: &HfProfile) -> Vec<u8> {
    let mut data = Vec::new();

    data.extend_from_slice(&(profile.total_runs as u64).to_le_bytes());
    data.extend_from_slice(&(profile.min_run as u32).to_le_bytes());
    data.extend_from_slice(&(profile.run_freq.len() as u64).to_le_bytes());

    for (&(base, run_len), &freq) in &profile.run_freq {
        data.push(base);
        data.extend_from_slice(&(run_len as u32).to_le_bytes());
        data.extend_from_slice(&freq.to_le_bytes());
    }

    data
}

/// Deserialize HfProfile from bytes.
pub fn deserialize_hf_profile(data: &[u8]) -> Option<HfProfile> {
    if data.len() < 16 {
        return None;
    }

    let mut pos = 0;
    let total_runs = u64::from_le_bytes(data[pos..pos + 8].try_into().ok()?) as usize;
    pos += 8;

    let min_run = u32::from_le_bytes(data[pos..pos + 4].try_into().ok()?) as usize;
    pos += 4;

    let num_entries = u64::from_le_bytes(data[pos..pos + 8].try_into().ok()?) as usize;
    pos += 8;

    let entry_size = 1 + 4 + 8;
    if pos + num_entries * entry_size > data.len() {
        return None;
    }

    let mut run_freq = HashMap::new();
    for _ in 0..num_entries {
        let base = data[pos];
        pos += 1;
        let run_len = u32::from_le_bytes(data[pos..pos + 4].try_into().ok()?) as usize;
        pos += 4;
        let freq = f64::from_le_bytes(data[pos..pos + 8].try_into().ok()?);
        pos += 8;
        run_freq.insert((base, run_len), freq);
    }

    Some(HfProfile {
        run_freq,
        total_runs,
        min_run,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode_sequence;

    #[test]
    fn test_compute_simple() {
        let seq = encode_sequence("AAAAAACCCGGGTTTTT");
        let profile = HfProfile::compute(&seq, 3);
        assert!(profile.total_runs > 0);
        assert!(profile.run_freq.contains_key(&(1, 6)));
        assert!(profile.run_freq.contains_key(&(2, 3)));
        assert!(profile.run_freq.contains_key(&(3, 3)));
        assert!(profile.run_freq.contains_key(&(4, 5)));
    }

    #[test]
    fn test_compute_no_runs() {
        let seq = encode_sequence("ACGTACGT");
        let profile = HfProfile::compute(&seq, 3);
        assert_eq!(profile.total_runs, 0);
    }

    #[test]
    fn test_similarity_identical() {
        let seq = encode_sequence("AAAAAACCCGGGTTTTTAAAAAACCC");
        let profile = HfProfile::compute(&seq, 3);
        let read = encode_sequence("AAAAAACCCGGGTTTTT");
        let read_prof = HfProfile::compute_read(&read, 3);
        let sim = profile.similarity(&read_prof);
        assert!(sim > 0.3);
    }

    #[test]
    fn test_similarity_different() {
        let seq1 = encode_sequence("AAAAAACCCGGGTTTTT");
        let profile1 = HfProfile::compute(&seq1, 3);
        let seq2 = encode_sequence("GGGGGGGGGGGGGGGGGGGG");
        let profile2 = HfProfile::compute(&seq2, 3);
        let sim = profile1.similarity(&profile2.run_freq);
        assert!(sim < 0.5);
    }

    #[test]
    fn test_serialize_deserialize() {
        let seq = encode_sequence("AAAAAACCCGGGTTTTT");
        let profile = HfProfile::compute(&seq, 3);
        let bytes = serialize_hf_profile(&profile);
        let loaded = deserialize_hf_profile(&bytes).unwrap();
        assert_eq!(loaded.total_runs, profile.total_runs);
        assert_eq!(loaded.min_run, profile.min_run);
        for key in profile.run_freq.keys() {
            assert!(
                (loaded.run_freq.get(key).unwrap() - profile.run_freq.get(key).unwrap()).abs()
                    < 1e-10
            );
        }
    }

    #[test]
    fn test_read_profile_computation() {
        let read = encode_sequence("AAAAACCCGGGGG");
        let profile = HfProfile::compute_read(&read, 3);
        assert!(profile.contains_key(&(1, 5)));
        assert!(profile.contains_key(&(2, 3)));
        assert!(profile.contains_key(&(3, 5)));
        let sum: f64 = profile.values().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }
}
