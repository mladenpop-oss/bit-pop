//! EM Algorithm for soft-assignment read classification.
//!
//! Improves strain-level resolution for highly similar genomes (evo_* strains)
//! by using probabilistic assignments instead of hard classification.
//!
//! Algorithm:
//!   Phase 1: Soft assignment - each read gets probabilities for ALL genomes
//!   Phase 2: EM iterations - E-step (compute probabilities) + M-step (update abundances)
//!   Phase 3: Back-propagate - hard assignment from final probabilities
//!
//! Usage:
//!   let em = EMClassifier::new(0.001, 50, 1e-6);
//!   let results = em.classify(&mappings);

use std::collections::{HashMap, HashSet};

/// Configuration for EM algorithm.
#[derive(Debug, Clone)]
pub struct EMConfig {
    /// Convergence threshold (KL divergence).
    pub convergence_threshold: f64,
    /// Maximum number of EM iterations.
    pub max_iterations: usize,
    /// Minimum abundance to keep a genome.
    pub min_abundance: f64,
    /// Softmax temperature for initial soft assignment.
    pub temperature: f64,
    /// Top-K genomes per read for EM.
    pub top_k: usize,
}

impl Default for EMConfig {
    fn default() -> Self {
        Self {
            convergence_threshold: 0.001,
            max_iterations: 50,
            min_abundance: 1e-6,
            temperature: 0.1,
            top_k: 10,
        }
    }
}

/// Soft assignment for a single read: genome -> probability.
pub type SoftAssignment = HashMap<String, f64>;

/// Mapping input: read name -> list of (genome_name, score).
pub type ReadMappings = Vec<(String, String, f64)>;

/// EM classifier for soft-assignment read classification.
pub struct EMClassifier {
    config: EMConfig,

    /// All genomes seen across all reads.
    all_genomes: HashSet<String>,

    /// Read ID -> soft assignment (genome -> probability).
    soft_assignments: Vec<(String, SoftAssignment)>,

    /// Genome abundances (theta).
    abundances: HashMap<String, f64>,

    /// Statistics: number of iterations run.
    pub iterations_run: usize,

    /// Final KL divergence.
    pub final_kl: f64,
}

impl EMClassifier {
    /// Create a new EM classifier with the given config.
    pub fn new(config: EMConfig) -> Self {
        Self {
            config,
            all_genomes: HashSet::new(),
            soft_assignments: Vec::new(),
            abundances: HashMap::new(),
            iterations_run: 0,
            final_kl: 0.0,
        }
    }

    /// Create with default config.
    pub fn default_config() -> Self {
        Self::new(EMConfig::default())
    }

    /// Run EM classification on the given mappings.
    ///
    /// # Arguments
    /// * `mappings` — list of (read_name, genome_name, score) tuples
    ///
    /// # Returns
    /// List of (read_name, best_genome) hard assignments.
    pub fn classify(&mut self, mappings: &ReadMappings) -> Vec<(String, Option<String>)> {
        // Phase 1: Initialize soft assignments
        self.initialize(mappings);

        // Phase 2: EM iterations
        let _stats = self.run_em();

        // Phase 3: Back-propagate to hard assignments
        self.back_propagate()
    }

    /// Initialize soft assignments from raw mappings.
    fn initialize(&mut self, mappings: &ReadMappings) {
        // Collect all genomes
        for (_, genome_name, _) in mappings {
            self.all_genomes.insert(genome_name.clone());
        }

        // Group mappings by read name
        let mut read_map: HashMap<String, Vec<(String, f64)>> = HashMap::new();
        for (read_name, genome_name, score) in mappings {
            read_map
                .entry(read_name.clone())
                .or_default()
                .push((genome_name.clone(), *score));
        }

        // Initialize soft assignments using softmax
        for (read_name, genome_scores) in read_map {
            let probs = normalize_to_probabilities(
                &genome_scores,
                self.config.temperature,
                self.config.top_k,
            );
            if !probs.is_empty() {
                self.soft_assignments.push((read_name, probs));
            }
        }

        // Initialize uniform abundances
        let n_genomes = self.all_genomes.len();
        if n_genomes > 0 {
            let uniform = 1.0 / n_genomes as f64;
            for genome in &self.all_genomes {
                self.abundances.insert(genome.clone(), uniform);
            }
        }
    }

    /// Run EM iterations (E-step + M-step).
    fn run_em(&mut self) -> EMStats {
        let mut prev_abundances: HashMap<String, f64> = self.abundances.clone();
        let mut kl_divergence = f64::INFINITY;

        for iteration in 1..=self.config.max_iterations {
            // E-step
            self.e_step();

            // M-step
            self.m_step();

            // Check convergence
            if iteration > 1 {
                kl_divergence = compute_kl_divergence(
                    &prev_abundances,
                    &self.abundances,
                    &self.all_genomes,
                    self.config.min_abundance,
                );
            }

            prev_abundances = self.abundances.clone();

            if kl_divergence < self.config.convergence_threshold {
                self.iterations_run = iteration;
                self.final_kl = kl_divergence;
                break;
            }

            if iteration == self.config.max_iterations {
                self.iterations_run = iteration;
                self.final_kl = kl_divergence;
            }
        }

        EMStats {
            iterations: self.iterations_run,
            final_kl: self.final_kl,
            active_genomes: self
                .abundances
                .values()
                .filter(|&&a| a > self.config.min_abundance)
                .count(),
        }
    }

    /// E-step: P(genome_i | read) = P(read | genome_i) * theta_i / Z
    fn e_step(&mut self) {
        for (_, probs) in &mut self.soft_assignments {
            if probs.is_empty() {
                continue;
            }

            // Apply abundance priors
            let weighted: HashMap<String, f64> = probs
                .iter()
                .map(|(genome, prob)| {
                    let abundance = self
                        .abundances
                        .get(genome)
                        .copied()
                        .unwrap_or(self.config.min_abundance);
                    (genome.clone(), prob * abundance)
                })
                .collect();

            // Normalize
            let total: f64 = weighted.values().sum();
            if total > 0.0 {
                *probs = weighted.into_iter().map(|(g, p)| (g, p / total)).collect();
            } else {
                // Fallback to uniform
                let n = probs.len();
                *probs = probs.keys().cloned().map(|g| (g, 1.0 / n as f64)).collect();
            }
        }
    }

    /// M-step: Update theta_i = sum(P(genome_i | read_r)) / total_reads
    fn m_step(&mut self) {
        // Sum probabilities for each genome
        let mut sums: HashMap<String, f64> = HashMap::new();
        let mut total_reads = 0usize;

        for (_, probs) in &self.soft_assignments {
            for (genome, prob) in probs {
                *sums.entry(genome.clone()).or_insert(0.0) += prob;
            }
            total_reads += 1;
        }

        if total_reads == 0 {
            return;
        }

        // Normalize to get new abundances
        let new_abundances: HashMap<String, f64> = self
            .all_genomes
            .iter()
            .filter_map(|genome| {
                let abundance = sums.get(genome).copied().unwrap_or(0.0) / total_reads as f64;
                if abundance >= self.config.min_abundance {
                    Some((genome.clone(), abundance))
                } else {
                    None
                }
            })
            .collect();

        // Renormalize
        let total: f64 = new_abundances.values().sum();
        if total > 0.0 {
            self.abundances = new_abundances
                .into_iter()
                .map(|(g, a)| (g, a / total))
                .collect();
        }
    }

    /// Convert soft assignments back to hard classifications.
    fn back_propagate(&self) -> Vec<(String, Option<String>)> {
        self.soft_assignments
            .iter()
            .map(|(read_name, probs)| {
                let best = probs
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
                (read_name.clone(), best.map(|(g, _)| g.clone()))
            })
            .collect()
    }

    /// Get sorted abundance report.
    pub fn get_abundance_report(&self) -> Vec<(String, f64)> {
        let mut sorted: Vec<_> = self
            .abundances
            .iter()
            .map(|(g, a)| (g.clone(), *a))
            .collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        sorted
    }
}

/// Statistics from EM run.
#[derive(Debug)]
pub struct EMStats {
    pub iterations: usize,
    pub final_kl: f64,
    pub active_genomes: usize,
}

/// Convert genome scores to probability distribution using softmax with temperature.
fn normalize_to_probabilities(
    genome_scores: &[(String, f64)],
    temperature: f64,
    top_k: usize,
) -> HashMap<String, f64> {
    if genome_scores.is_empty() {
        return HashMap::new();
    }

    // Sort by score descending
    let mut sorted: Vec<_> = genome_scores.iter().collect();
    sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Take top-K
    let top: Vec<_> = sorted.iter().take(top_k).copied().collect();

    // Filter by minimum score
    let top: Vec<_> = top.into_iter().filter(|&(_, s)| *s >= 0.0).collect();

    if top.is_empty() {
        return HashMap::new();
    }

    // Softmax with temperature
    let max_logit = top
        .iter()
        .map(|&(_, s)| s)
        .fold(f64::NEG_INFINITY, |a, b| f64::max(a, *b));

    let mut exps: HashMap<String, f64> = HashMap::new();
    let mut total = 0.0;

    for &(genome, score) in &top {
        let adjusted = (score - max_logit) / temperature;
        let exp_val = (adjusted.max(-50.0).min(50.0)).exp();
        total += exp_val;
        exps.insert(genome.clone(), exp_val);
    }

    if total == 0.0 {
        // Uniform distribution as fallback
        let n = top.len();
        exps.into_iter().map(|(g, _)| (g, 1.0 / n as f64)).collect()
    } else {
        exps.into_iter().map(|(g, e)| (g, e / total)).collect()
    }
}

/// Compute KL divergence: KL(p || q) = sum(p * log(p/q))
fn compute_kl_divergence(
    p: &HashMap<String, f64>,
    q: &HashMap<String, f64>,
    all_genomes: &HashSet<String>,
    min_abundance: f64,
) -> f64 {
    let mut kl = 0.0;
    for genome in all_genomes {
        let p_val = p.get(genome).copied().unwrap_or(min_abundance);
        let q_val = q.get(genome).copied().unwrap_or(min_abundance);
        if p_val > 0.0 {
            kl += p_val * (p_val / q_val).ln();
        }
    }
    kl
}
