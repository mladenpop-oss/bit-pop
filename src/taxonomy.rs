//! NCBI Taxonomy tree parsing and LCA (Lowest Common Ancestor) algorithm.
//!
//! Parses NCBI taxonomy dump files (nodes.dmp, names.dmp) to build a hierarchical
//! taxonomy tree, then uses LCA to aggregate genome-level classifications into
//! taxonomic rank profiles (species, genus, family, phylum, etc.).
//!
//! # NCBI taxonomy file formats
//!
//! **nodes.dmp** (tab-separated):
//!   taxon_id \t parent_id \t rank \t division_id \t ...
//!
//! **names.dmp** (tab-separated):
//!   taxon_id \t name_name \t name_class \t ...

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;

/// A single node in the taxonomy tree.
#[derive(Debug, Clone)]
pub struct TaxNode {
    /// NCBI taxonomy ID
    pub tax_id: u32,
    /// Parent taxonomy ID (0 = root)
    pub parent_id: u32,
    /// Taxonomic rank (e.g., "species", "genus", "phylum")
    pub rank: String,
    /// Scientific name (from names.dmp, "scientific name" class)
    pub name: String,
}

/// Taxonomy tree: tax_id -> TaxNode mapping + cache for LCA queries.
pub struct TaxonomyTree {
    /// All taxonomy nodes
    nodes: HashMap<u32, TaxNode>,
    /// Reverse index: name -> tax_id
    name_to_id: HashMap<String, u32>,
    /// Genome name -> tax_id mapping (set by user)
    genome_tax_ids: HashMap<String, u32>,
    /// Standard taxonomic ranks in order from root to leaf
    pub standard_ranks: Vec<String>,
}

impl Default for TaxonomyTree {
    fn default() -> Self {
        Self::new()
    }
}

impl TaxonomyTree {
    /// Create an empty taxonomy tree.
    pub fn new() -> Self {
        Self {
            nodes: HashMap::new(),
            name_to_id: HashMap::new(),
            genome_tax_ids: HashMap::new(),
            standard_ranks: vec![
                "no rank".to_string(),
                "superkingdom".to_string(),
                "phylum".to_string(),
                "class".to_string(),
                "order".to_string(),
                "family".to_string(),
                "genus".to_string(),
                "species".to_string(),
            ],
        }
    }

    /// Load taxonomy from NCBI dump files.
    ///
    /// # Arguments
    /// * `nodes_path` — Path to nodes.dmp
    /// * `names_path` — Path to names.dmp
    pub fn load(nodes_path: &str, names_path: &str) -> Result<Self, TaxError> {
        let mut self_ = Self::new();
        self_.parse_nodes(nodes_path)?;
        self_.parse_names(names_path)?;
        Ok(self_)
    }

    /// Parse nodes.dmp file.
    fn parse_nodes(&mut self, path: &str) -> Result<(), TaxError> {
        let content = fs::read_to_string(path).map_err(TaxError::Io)?;

        for line in content.lines() {
            // Format: taxon_id | parent_id | rank | ...
            // Fields are separated by '\t|\t'
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                continue;
            }

            let tax_id: u32 = parts[0]
                .trim()
                .parse()
                .map_err(|_| TaxError::Parse(format!("Invalid tax_id: {}", parts[0].trim())))?;
            let parent_id: u32 = parts[1]
                .trim()
                .parse()
                .map_err(|_| TaxError::Parse(format!("Invalid parent_id: {}", parts[1].trim())))?;
            let rank = parts[2].trim().to_string();

            // Remove trailing '|' artifacts
            let rank = rank.trim_end_matches('|').trim().to_string();

            let node = TaxNode {
                tax_id,
                parent_id,
                rank,
                name: String::new(), // Will be filled from names.dmp
            };

            self.nodes.insert(tax_id, node);
        }

        Ok(())
    }

    /// Parse names.dmp file, filling in scientific names.
    fn parse_names(&mut self, path: &str) -> Result<(), TaxError> {
        let content = fs::read_to_string(path).map_err(TaxError::Io)?;

        for line in content.lines() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                continue;
            }

            let tax_id: u32 = parts[0].trim().parse().map_err(|_| {
                TaxError::Parse(format!("Invalid tax_id in names: {}", parts[0].trim()))
            })?;
            let name_name = parts[1].trim().to_string();
            let name_class = parts[2].trim().to_string();
            let name_class = name_class.trim_end_matches('|').trim().to_string();

            // Only use "scientific name" class
            if name_class == "scientific name" {
                if let Some(node) = self.nodes.get_mut(&tax_id) {
                    node.name = name_name.clone();
                    self.name_to_id.insert(name_name, tax_id);
                }
            }
        }

        Ok(())
    }

    /// Map a genome name to a taxonomy ID.
    ///
    /// Tries exact match first, then substring match against known taxonomy names.
    ///
    /// # Arguments
    /// * `genome_name` — Name of the genome (e.g., "Escherichia coli")
    pub fn map_genome_name(&mut self, genome_name: &str) -> Option<u32> {
        // Try exact match
        if let Some(&tax_id) = self.name_to_id.get(genome_name) {
            self.genome_tax_ids.insert(genome_name.to_string(), tax_id);
            return Some(tax_id);
        }

        // Try case-insensitive exact match
        let lower = genome_name.to_lowercase();
        for (name, &tax_id) in &self.name_to_id {
            if name.to_lowercase() == lower {
                self.genome_tax_ids.insert(genome_name.to_string(), tax_id);
                return Some(tax_id);
            }
        }

        // Try substring match (genome name contains tax name or vice versa)
        let best = self.name_to_id.iter().find(|(name, _)| {
            let n_lower = name.to_lowercase();
            n_lower.contains(&lower) || lower.contains(&n_lower)
        });

        if let Some((_, &tax_id)) = best {
            self.genome_tax_ids.insert(genome_name.to_string(), tax_id);
            Some(tax_id)
        } else {
            None
        }
    }

    /// Set taxonomy ID for a genome directly.
    pub fn set_genome_tax_id(&mut self, genome_name: &str, tax_id: u32) {
        self.genome_tax_ids.insert(genome_name.to_string(), tax_id);
    }

    /// Get taxonomy ID for a genome.
    pub fn genome_tax_id(&self, genome_name: &str) -> Option<u32> {
        self.genome_tax_ids.get(genome_name).copied()
    }

    /// Get all mapped genome taxonomy IDs.
    pub fn genome_tax_ids(&self) -> &HashMap<String, u32> {
        &self.genome_tax_ids
    }

    /// Get a taxonomy node by ID.
    pub fn node(&self, tax_id: u32) -> Option<&TaxNode> {
        self.nodes.get(&tax_id)
    }

    /// Get the parent of a taxonomy node.
    pub fn parent(&self, tax_id: u32) -> Option<&TaxNode> {
        let node = self.nodes.get(&tax_id)?;
        if node.parent_id == 0 {
            return None;
        }
        self.nodes.get(&node.parent_id)
    }

    /// Walk from a tax_id to the root, collecting all ancestor IDs (inclusive).
    pub fn ancestors(&self, tax_id: u32) -> Vec<u32> {
        let mut result = Vec::new();
        let mut current = tax_id;
        while let Some(node) = self.nodes.get(&current) {
            result.push(current);
            if node.parent_id == 0 {
                break;
            }
            current = node.parent_id;
        }
        result
    }

    /// Compute the Lowest Common Ancestor of two taxonomy IDs.
    pub fn lca(&self, tax_id_a: u32, tax_id_b: u32) -> Option<u32> {
        let ancestors_a: HashSet<u32> = self.ancestors(tax_id_a).into_iter().collect();
        let ancestors_b = self.ancestors(tax_id_b);

        // Find the deepest common ancestor (smallest tax_id with highest depth)
        let mut best: Option<u32> = None;
        let mut best_depth: usize = 0;

        for &id in &ancestors_b {
            if ancestors_a.contains(&id) {
                let depth = self.ancestors(id).len();
                if best.is_none() || depth > best_depth {
                    best = Some(id);
                    best_depth = depth;
                }
            }
        }

        best
    }

    /// Get the taxonomic path from root to a given tax_id.
    pub fn tax_path(&self, tax_id: u32) -> Vec<&TaxNode> {
        let ancestor_ids = self.ancestors(tax_id);
        ancestor_ids
            .into_iter()
            .rev()
            .filter_map(|id| self.nodes.get(&id))
            .collect()
    }

    /// Get the species name for a tax_id (walks up to nearest species rank).
    pub fn species_name(&self, tax_id: u32) -> Option<String> {
        let mut current = tax_id;
        while let Some(node) = self.nodes.get(&current) {
            if node.rank == "species" {
                return Some(node.name.clone());
            }
            if node.parent_id == 0 {
                break;
            }
            current = node.parent_id;
        }
        None
    }

    /// Get the genus name for a tax_id.
    pub fn genus_name(&self, tax_id: u32) -> Option<String> {
        let mut current = tax_id;
        while let Some(node) = self.nodes.get(&current) {
            if node.rank == "genus" {
                return Some(node.name.clone());
            }
            if node.parent_id == 0 {
                break;
            }
            current = node.parent_id;
        }
        None
    }

    /// Get the rank name for a tax_id.
    pub fn rank_name(&self, tax_id: u32) -> Option<&str> {
        self.nodes.get(&tax_id).map(|n| n.rank.as_str())
    }

    /// Get total number of nodes in the tree.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

/// Taxonomic abundance report: rank -> (name, count).
#[derive(Debug, Clone)]
pub struct TaxReport {
    /// Total reads classified
    pub total_reads: usize,
    /// Per-rank breakdown: rank -> (tax_name, read_count, percentage)
    pub ranks: BTreeMap<String, Vec<(String, usize, f64)>>,
    /// Per-genome taxonomic path with counts
    pub genomes: Vec<(String, String, String, usize)>,
}

/// Compute taxonomic report from genome-level read counts.
///
/// # Arguments
/// * `tree` — Taxonomy tree with genome mappings
/// * `genome_counts` — genome_name -> read_count
pub fn compute_tax_report(
    tree: &TaxonomyTree,
    genome_counts: &HashMap<String, usize>,
) -> TaxReport {
    let total_reads: usize = genome_counts.values().sum();

    // Aggregate reads by taxonomy node using LCA
    let mut tax_counts: HashMap<u32, usize> = HashMap::new();

    for (genome_name, &count) in genome_counts {
        if let Some(&tax_id) = tree.genome_tax_ids().get(genome_name) {
            // Add count to all ancestors
            for &ancestor_id in &tree.ancestors(tax_id) {
                *tax_counts.entry(ancestor_id).or_insert(0) += count;
            }
        }
    }

    // Build per-rank breakdown
    let mut ranks: BTreeMap<String, Vec<(String, usize, f64)>> = BTreeMap::new();

    for (&tax_id, &count) in &tax_counts {
        if let Some(node) = tree.node(tax_id) {
            let pct = if total_reads > 0 {
                (count as f64 / total_reads as f64) * 100.0
            } else {
                0.0
            };
            ranks
                .entry(node.rank.clone())
                .or_default()
                .push((node.name.clone(), count, pct));
        }
    }

    // Sort each rank by count descending
    for entry in ranks.values_mut() {
        entry.sort_by_key(|b| std::cmp::Reverse(b.1));
    }

    // Build genome-level taxonomic paths
    let mut genomes: Vec<(String, String, String, usize)> = Vec::new();
    for (genome_name, &count) in genome_counts {
        if let Some(&tax_id) = tree.genome_tax_ids().get(genome_name) {
            let genus = tree.genus_name(tax_id).unwrap_or_else(|| "N/A".to_string());
            let species = tree
                .species_name(tax_id)
                .unwrap_or_else(|| genome_name.clone());
            genomes.push((genome_name.clone(), genus, species, count));
        }
    }
    genomes.sort_by_key(|b| std::cmp::Reverse(b.3));

    TaxReport {
        total_reads,
        ranks,
        genomes,
    }
}

/// Format taxonomic report as a human-readable string.
pub fn format_tax_report(report: &TaxReport, top_n: usize) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "=== Taxonomic Report ({} total reads) ===\n\n",
        report.total_reads
    ));

    for (rank, entries) in &report.ranks {
        // Skip "no rank" entries
        if rank == "no rank" {
            continue;
        }
        out.push_str(&format!("{}\n", rank.to_uppercase()));
        out.push_str(&format!("{:<40} {:>10} {:>10}\n", "Name", "Reads", "%"));
        out.push_str(&format!("{:-<64}\n", ""));

        for (name, count, pct) in entries.iter().take(top_n) {
            out.push_str(&format!("{:<40} {:>10} {:>9.1}%\n", name, count, pct));
        }
        out.push('\n');
    }

    // Genome-level detail
    out.push_str("GENOME DETAIL\n");
    out.push_str(&format!(
        "{:<35} {:<25} {:<25} {:>10}\n",
        "Genome", "Genus", "Species", "Reads"
    ));
    out.push_str(&format!("{:-<96}\n", ""));

    for (genome, genus, species, count) in &report.genomes {
        out.push_str(&format!(
            "{:<35} {:<25} {:<25} {:>10}\n",
            genome, genus, species, count
        ));
    }

    out
}

/// Aggregate genome-level EM abundances to taxonomy level using LCA.
///
/// Takes genome -> probability mappings and aggregates them to the
/// lowest common ancestor for each read across all candidate genomes.
///
/// # Arguments
/// * `tree` — Taxonomy tree with genome mappings
/// * `read_candidates` — For each read: list of (genome_name, probability)
/// * `top_k` — Number of top candidates per read to consider
///
/// # Returns
/// tax_id -> aggregated probability
pub fn lca_aggregate(
    tree: &TaxonomyTree,
    read_candidates: &[(String, Vec<(String, f64)>)],
    top_k: usize,
) -> HashMap<u32, f64> {
    let mut tax_probs: HashMap<u32, f64> = HashMap::new();

    for (_read_name, candidates) in read_candidates {
        // Sort by probability descending, take top-K
        let mut sorted: Vec<_> = candidates.iter().collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let top: Vec<_> = sorted.iter().take(top_k).collect();

        if top.is_empty() {
            continue;
        }

        // Get tax IDs for top candidates
        let tax_ids: Vec<u32> = top
            .iter()
            .filter_map(|(genome, _)| tree.genome_tax_ids().get(genome).copied())
            .collect();

        if tax_ids.is_empty() {
            continue;
        }

        // Compute LCA of all top candidates
        let mut lca_id = tax_ids[0];
        for &tid in &tax_ids[1..] {
            if let Some(new_lca) = tree.lca(lca_id, tid) {
                lca_id = new_lca;
            }
        }

        // Sum the probabilities of top candidates
        let total_prob: f64 = top.iter().map(|(_, p)| *p).sum();
        *tax_probs.entry(lca_id).or_insert(0.0) += total_prob;
    }

    tax_probs
}

#[derive(Debug, thiserror::Error)]
pub enum TaxError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_tree() -> TaxonomyTree {
        // Build a small test taxonomy tree manually
        let mut tree = TaxonomyTree::new();

        // Root -> superkingdom -> phylum -> class -> order -> family -> genus -> species
        let nodes = vec![
            (1, 0, "no rank", "root"),
            (2, 1, "superkingdom", "Bacteria"),
            (3, 2, "phylum", "Proteobacteria"),
            (4, 3, "class", "Gammaproteobacteria"),
            (5, 4, "order", "Enterobacterales"),
            (6, 5, "family", "Enterobacteriaceae"),
            (7, 6, "genus", "Escherichia"),
            (8, 7, "species", "Escherichia coli"),
            // Second species in same genus
            (9, 7, "species", "Escherichia fergusonii"),
            // Different genus, same family
            (10, 6, "genus", "Salmonella"),
            (11, 10, "species", "Salmonella enterica"),
        ];

        for (tax_id, parent_id, rank, name) in nodes {
            tree.nodes.insert(
                tax_id,
                TaxNode {
                    tax_id,
                    parent_id,
                    rank: rank.to_string(),
                    name: name.to_string(),
                },
            );
            tree.name_to_id.insert(name.to_string(), tax_id);
        }

        tree
    }

    #[test]
    fn test_ancestors() {
        let tree = create_test_tree();
        let ancestors = tree.ancestors(8); // E. coli
        assert_eq!(ancestors, vec![8, 7, 6, 5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_lca_same_genus() {
        let tree = create_test_tree();
        // E. coli (8) and E. fergusonii (9) should have LCA = Escherichia (7)
        let lca = tree.lca(8, 9);
        assert_eq!(lca, Some(7));
    }

    #[test]
    fn test_lca_different_genus() {
        let tree = create_test_tree();
        // E. coli (8) and S. enterica (11) should have LCA = Enterobacteriaceae (6)
        let lca = tree.lca(8, 11);
        assert_eq!(lca, Some(6));
    }

    #[test]
    fn test_lca_same_node() {
        let tree = create_test_tree();
        let lca = tree.lca(8, 8);
        assert_eq!(lca, Some(8));
    }

    #[test]
    fn test_genus_name() {
        let tree = create_test_tree();
        assert_eq!(tree.genus_name(8), Some("Escherichia".to_string()));
        assert_eq!(tree.genus_name(11), Some("Salmonella".to_string()));
    }

    #[test]
    fn test_species_name() {
        let tree = create_test_tree();
        assert_eq!(tree.species_name(8), Some("Escherichia coli".to_string()));
    }

    #[test]
    fn test_map_genome_name() {
        let mut tree = create_test_tree();
        assert_eq!(tree.map_genome_name("Escherichia coli"), Some(8));
        assert_eq!(tree.map_genome_name("Salmonella enterica"), Some(11));
        assert_eq!(tree.map_genome_name("Unknown organism"), None);
    }

    #[test]
    fn test_tax_report() {
        let mut tree = create_test_tree();
        tree.map_genome_name("Escherichia coli");
        tree.map_genome_name("Salmonella enterica");

        let mut counts = HashMap::new();
        counts.insert("Escherichia coli".to_string(), 100);
        counts.insert("Salmonella enterica".to_string(), 50);

        let report = compute_tax_report(&tree, &counts);
        assert_eq!(report.total_reads, 150);
        assert!(report.ranks.contains_key("family"));
    }
}
