use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use rayon::prelude::*;

use crate::fastq::ReadsFormat;
use crate::report_atomic_progress;
use crate::BitPop;

/// Consensus voting strategy for multi-k mapping.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ConsensusStrategy {
    Majority,
    #[default]
    WeightedScore,
    /// Union mode: take best score from any k, no voting
    BestScore,
    /// Raw score sum (no k-weight), like JNI Android
    BaseScore,
}

/// Per-k mapping result for a single read.
#[derive(Debug, Clone)]
pub struct KResult {
    pub k: usize,
    pub genome_id: u32,
    pub genome_name: String,
    pub score: f64,
    pub rarity: f64,
    pub cigar: String,
    pub position: u64,
    pub is_reverse: bool,
}

/// Consensus result for a single read across multiple k-values.
#[derive(Debug, Clone)]
pub struct ConsensusResult {
    pub genome_id: u32,
    pub genome_name: String,
    pub vote_count: usize,
    pub k_count: usize,
    pub consensus_score: f64,
    pub k_results: Vec<KResult>,
}

/// Multi-k consensus mapper.
pub struct MultiKConsensus {
    pub indexes: HashMap<usize, BitPop>,
    pub k_values: Vec<usize>,
    pub strategy: ConsensusStrategy,
    pub min_agreement: usize,
    pub min_k_mappings: usize,
    pub k_weights: HashMap<usize, f64>,
    pub min_score: f64,
    pub context_window: usize,
    pub chunk_size: usize,
    pub chunk_pct: f64,
    pub chunk_min: usize,
    pub chunk_max: usize,
    pub enable_snp_detect: bool,
    pub snp_min_support: u32,
    pub snp_penalty: f64,
    pub top_n: usize,
    /// top_n for each BitPop index (controls how many rare k-mers are used as anchors)
    pub map_top_n: usize,
}

impl MultiKConsensus {
    pub fn from_paths(index_paths: &[PathBuf], min_score: f64) -> Result<Self, String> {
        let mut indexes = HashMap::new();
        let mut k_values = Vec::new();

        for path in index_paths {
            let bp = BitPop::deserialize_from_file(path.to_str().unwrap())
                .map_err(|e| format!("Failed to load index {}: {}", path.display(), e))?;
            let k = bp.k();
            println!("  Loading k={} index: {}", k, path.display());
            indexes.insert(k, bp);
            k_values.push(k);
        }

        k_values.sort();

        let mut k_weights = HashMap::new();
        let min_k = *k_values.iter().min().unwrap_or(&10);
        for &k in &k_values {
            k_weights.insert(k, k as f64 / min_k as f64);
        }

        Ok(Self {
            indexes,
            k_values,
            strategy: ConsensusStrategy::default(),
            min_agreement: 0,
            min_k_mappings: 1,
            k_weights,
            min_score,
            context_window: 0,
            chunk_size: 0,
            chunk_pct: 0.0,
            chunk_min: 20,
            chunk_max: 500,
            enable_snp_detect: false,
            snp_min_support: 3,
            snp_penalty: 0.1,
            top_n: 1,
            map_top_n: 1,
        })
    }

    pub fn map_read(&self, read_seq: &str, _read_qual: &str) -> Vec<ConsensusResult> {
        let mut all_k_results: Vec<KResult> = Vec::new();

        for &k in &self.k_values {
            let bp = self.indexes.get(&k).unwrap();

            let results = bp.map_read(read_seq, self.context_window);

            if let Some(best) = results.first() {
                let genome_name = bp.genome_name(best.genome_id).unwrap_or("?").to_string();
                all_k_results.push(KResult {
                    k,
                    genome_id: best.genome_id,
                    genome_name,
                    score: best.score,
                    rarity: best.rarity,
                    cigar: best.cigar.clone(),
                    position: best.position,
                    is_reverse: best.is_reverse,
                });
            }
        }

        if all_k_results.is_empty() {
            return Vec::new();
        }

        // BestScore: take best from any k, no voting
        if self.strategy == ConsensusStrategy::BestScore {
            let mut best_idx = 0;
            for i in 1..all_k_results.len() {
                if all_k_results[i].score > all_k_results[best_idx].score {
                    best_idx = i;
                }
            }
            let kr = &all_k_results[best_idx];
            return vec![ConsensusResult {
                genome_id: kr.genome_id,
                genome_name: kr.genome_name.clone(),
                vote_count: 1,
                k_count: all_k_results.len(),
                consensus_score: kr.score,
                k_results: vec![kr.clone()],
            }];
        }

        // Require at least min_k_mappings k-values to find a mapping
        if all_k_results.len() < self.min_k_mappings {
            return Vec::new();
        }

        if self.min_score > 0.0 {
            all_k_results.retain(|r| r.score >= self.min_score);
        }

        if all_k_results.is_empty() {
            return Vec::new();
        }

        let candidates = match self.strategy {
            ConsensusStrategy::Majority => self.vote_majority_multi(&all_k_results),
            ConsensusStrategy::WeightedScore | ConsensusStrategy::BestScore => {
                self.vote_weighted_multi(&all_k_results)
            }
            ConsensusStrategy::BaseScore => self.vote_base_score(&all_k_results),
        };

        let limit = if self.top_n > 1 { self.top_n } else { 1 };
        candidates.into_iter().take(limit).collect()
    }

    fn vote_majority_multi(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
        let mut genome_votes: HashMap<u32, (usize, f64, String, Vec<KResult>)> = HashMap::new();

        for kr in k_results {
            let entry = genome_votes.entry(kr.genome_id).or_insert((
                0,
                0.0,
                kr.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += kr.score;
            entry.3.push(kr.clone());
        }

        let mut candidates: Vec<_> = genome_votes
            .into_iter()
            .map(|(id, (vote_count, total_score, name, kr))| {
                (id, vote_count, total_score, name, kr)
            })
            .collect();
        candidates.sort_by(|a, b| {
            a.1.cmp(&b.1)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });

        candidates
            .into_iter()
            .map(|(id, vote_count, total_score, name, kr)| ConsensusResult {
                genome_id: id,
                genome_name: name,
                vote_count,
                k_count: k_results.len(),
                consensus_score: total_score / k_results.len() as f64,
                k_results: kr,
            })
            .collect()
    }

    fn vote_weighted_multi(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
        let mut genome_scores: HashMap<u32, (usize, f64, String, Vec<KResult>)> = HashMap::new();

        for kr in k_results {
            let weight = self.k_weights.get(&kr.k).copied().unwrap_or(1.0);
            let weighted_score = kr.score * weight;
            let entry = genome_scores.entry(kr.genome_id).or_insert((
                0,
                0.0,
                kr.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += weighted_score;
            entry.3.push(kr.clone());
        }

        let mut candidates: Vec<_> = genome_scores
            .into_iter()
            .map(|(id, (vote_count, total_score, name, kr))| {
                (id, vote_count, total_score, name, kr)
            })
            .collect();
        candidates.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        candidates
            .into_iter()
            .map(|(_id, vote_count, total_score, name, kr)| ConsensusResult {
                genome_id: _id,
                genome_name: name,
                vote_count,
                k_count: k_results.len(),
                consensus_score: total_score / k_results.len() as f64,
                k_results: kr,
            })
            .collect()
    }

    fn vote_base_score(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
        let mut genome_scores: HashMap<u32, (usize, f64, String, Vec<KResult>)> = HashMap::new();

        for kr in k_results {
            let entry = genome_scores.entry(kr.genome_id).or_insert((
                0,
                0.0,
                kr.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += kr.score;
            entry.3.push(kr.clone());
        }

        let mut candidates: Vec<_> = genome_scores
            .into_iter()
            .map(|(id, (vote_count, total_score, name, kr))| {
                (id, vote_count, total_score, name, kr)
            })
            .collect();
        candidates.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        candidates
            .into_iter()
            .map(|(id, vote_count, total_score, name, kr)| ConsensusResult {
                genome_id: id,
                genome_name: name,
                vote_count,
                k_count: k_results.len(),
                consensus_score: total_score / k_results.len() as f64,
                k_results: kr,
            })
            .collect()
    }

    pub fn map_reads_to_sam(
        &self,
        reads_path: &std::path::Path,
        output_path: &std::path::Path,
        threads: usize,
    ) -> Result<(usize, usize), String> {
        let reads = ReadsFormat::Fastq(
            crate::fastq::parse_fastq(reads_path.to_str().unwrap())
                .map_err(|e| format!("Failed to parse reads: {}", e))?,
        );

        let total = reads.count();
        println!("  Loaded {} reads", total);

        // Collect all genome names for SAM header
        let mut all_genomes: Vec<(String, u64)> = Vec::new();

        // Just collect from first index (all indexes have same genomes)
        if let Some(first_bp) = self.indexes.values().next() {
            let names = first_bp.genome_names_ordered();
            for (gid, name) in names.iter().enumerate() {
                let len = first_bp.genome_seq_len(gid as u32).unwrap_or(0) as u64;
                all_genomes.push((name.clone(), len));
            }
        }

        let file = File::create(output_path)
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        for (name, len) in &all_genomes {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, len)
                .map_err(|e| format!("Failed to write SQ line: {}", e))?;
        }

        println!("\n[1/2] Mapping reads with multi-k consensus...");
        let pb = indicatif::ProgressBar::new(total as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
                .progress_chars("#>-"),
        );

        // Convert to vector for parallel iteration
        let reads_vec: Vec<_> = match &reads {
            ReadsFormat::Fastq(r) => r
                .iter()
                .map(|(name, seq, qual)| (name.clone(), seq.clone(), qual.to_vec()))
                .collect(),
            ReadsFormat::Fasta(r) => r
                .iter()
                .map(|(name, seq)| (name.clone(), seq.clone(), Vec::new()))
                .collect(),
        };

        let mut mapped = 0usize;

        if threads > 1 {
            let results: Vec<_> = reads_vec
                .into_par_iter()
                .map(|(name, seq, _qual)| {
                    let candidates = self.map_read(&seq, "");
                    if !candidates.is_empty() {
                        Some((name, seq, candidates))
                    } else {
                        None
                    }
                })
                .collect();

            for result in results {
                if let Some((name, seq, candidates)) = result {
                    for (i, cr) in candidates.iter().enumerate() {
                        self.write_consensus_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                        mapped += 1;
                    }
                }
                pb.inc(1);
                report_atomic_progress(pb.position(), total as u64);
            }
        } else {
            for (name, seq, _qual) in &reads_vec {
                let candidates = self.map_read(seq, "");
                if !candidates.is_empty() {
                    for (i, cr) in candidates.iter().enumerate() {
                        self.write_consensus_sam_line(&mut writer, name, seq, cr, i == 0)?;
                        mapped += 1;
                    }
                }
                pb.inc(1);
                report_atomic_progress(pb.position(), total as u64);
            }
        }

        pb.finish();
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output: {}", e))?;

        Ok((mapped, total))
    }

    /// Stream-based mapping: process reads in chunks to limit memory usage.
    pub fn map_reads_to_sam_stream(
        &self,
        reads_path: &std::path::Path,
        output_path: &std::path::Path,
        threads: usize,
        chunk_size: usize,
    ) -> Result<(usize, usize), String> {
        use crate::fastq::FastqChunkParser;

        let total = FastqChunkParser::count_reads(reads_path.to_str().unwrap())
            .map_err(|e| format!("Failed to count reads: {}", e))?;
        println!("  Total reads: {} (streaming, chunk={})", total, chunk_size);

        let mut parser = FastqChunkParser::new(reads_path.to_str().unwrap(), chunk_size)
            .map_err(|e| format!("Failed to open FASTQ: {}", e))?;

        // Collect genomes from first index
        let mut all_genomes: Vec<(String, u64)> = Vec::new();
        if let Some(first_bp) = self.indexes.values().next() {
            let names = first_bp.genome_names_ordered();
            for (gid, name) in names.iter().enumerate() {
                let len = first_bp.genome_seq_len(gid as u32).unwrap_or(0) as u64;
                all_genomes.push((name.clone(), len));
            }
        }

        let file = File::create(output_path)
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        for (name, len) in &all_genomes {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, len)
                .map_err(|e| format!("Failed to write SQ line: {}", e))?;
        }

        println!("\n[1/2] Mapping reads with multi-k consensus (streaming)...");
        let pb = indicatif::ProgressBar::new(total as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
                .progress_chars("#>-"),
        );

        let mut mapped = 0usize;
        let mut chunk_num = 0usize;

        while let Some(chunk) = parser
            .next_chunk()
            .map_err(|e| format!("Parse error: {}", e))?
        {
            chunk_num += 1;
            let chunk_len = chunk.len();
            if chunk_len == 0 {
                break;
            }

            if threads > 1 {
                let results: Vec<_> = chunk
                    .into_par_iter()
                    .map(|(name, seq, _qual)| {
                        let candidates = self.map_read(&seq, "");
                        if !candidates.is_empty() {
                            Some((name, seq, candidates))
                        } else {
                            None
                        }
                    })
                    .collect();

                for result in results {
                    if let Some((name, seq, candidates)) = result {
                        for (i, cr) in candidates.iter().enumerate() {
                            self.write_consensus_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                            mapped += 1;
                        }
                    }
                    pb.inc(1);
                    report_atomic_progress(pb.position(), total as u64);
                }
            } else {
                for (name, seq, _qual) in chunk {
                    let candidates = self.map_read(&seq, "");
                    if !candidates.is_empty() {
                        for (i, cr) in candidates.iter().enumerate() {
                            self.write_consensus_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                            mapped += 1;
                        }
                    }
                    pb.inc(1);
                    report_atomic_progress(pb.position(), total as u64);
                }
            }

            println!("  Chunk {}: {} reads", chunk_num, chunk_len);
        }

        pb.finish();
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output: {}", e))?;

        Ok((mapped, total))
    }

    fn write_consensus_sam_line(
        &self,
        writer: &mut BufWriter<File>,
        read_name: &str,
        read_seq: &str,
        cr: &ConsensusResult,
        is_primary: bool,
    ) -> Result<(), String> {
        let mut flag: u16 = if cr.k_results.is_empty() {
            0x4
        } else if cr.k_results[0].is_reverse {
            0x10
        } else {
            0
        };

        if !is_primary {
            flag |= 0x800;
        }

        let mapq = (cr.consensus_score * 60.0) as u16;
        let pos = if cr.k_results.is_empty() {
            0
        } else {
            (cr.k_results[0].position + 1) as u32
        };
        let cigar = if cr.k_results.is_empty() {
            "*".to_string()
        } else {
            cr.k_results[0].cigar.clone()
        };

        let nm = if cr.k_results.is_empty() {
            0u32
        } else {
            self.compute_nm_from_cigar(&cr.k_results[0].cigar)
        };

        let mut rk_tags = String::new();
        for kr in &cr.k_results {
            rk_tags.push_str(&format!("\tRK{}:Z:{}", kr.k, kr.genome_name));
        }

        let kv_tag = format!("\tKV:Z:{}", cr.genome_name);
        let kc_tag = format!("\tKC:i:{}", cr.vote_count);
        let kk_tag = format!("\tKK:i:{}", cr.k_count);
        let as_tag = format!("\tAS:f:{:.4}", cr.consensus_score);
        let xs_tag = if is_primary {
            String::new()
        } else {
            format!("\tXS:f:{:.4}", cr.consensus_score)
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tNM:i:{}\tRK:Z:{}{}{}{}{}{}",
            read_name,
            flag,
            cr.genome_name,
            pos,
            mapq,
            cigar,
            read_seq,
            nm,
            rk_tags,
            kv_tag,
            kc_tag,
            kk_tag,
            as_tag,
            xs_tag,
        )
        .map_err(|e| format!("Failed to write SAM line: {}", e))?;
        Ok(())
    }

    /// Two-pass mapping: map each k separately (like Python script), then combine.
    /// Much faster because each k's mapping is fully parallelized.
    pub fn map_reads_to_sam_two_pass(
        &self,
        reads_path: &std::path::Path,
        output_path: &std::path::Path,
        threads: usize,
        _use_temp_files: bool,
    ) -> Result<(usize, usize), String> {
        use rayon::prelude::*;

        let reads = ReadsFormat::Fastq(
            crate::fastq::parse_fastq(reads_path.to_str().unwrap())
                .map_err(|e| format!("Failed to parse reads: {}", e))?,
        );

        let total = reads.count();
        println!("  Loaded {} reads", total);

        // Convert to vector
        let reads_vec: Vec<_> = match &reads {
            ReadsFormat::Fastq(r) => r
                .iter()
                .map(|(name, seq, qual)| (name.clone(), seq.clone(), qual.to_vec()))
                .collect(),
            ReadsFormat::Fasta(r) => r
                .iter()
                .map(|(name, seq)| (name.clone(), seq.clone(), Vec::new()))
                .collect(),
        };

        // Collect genomes from first index for SAM header
        let mut all_genomes: Vec<(String, u64)> = Vec::new();
        if let Some(first_bp) = self.indexes.values().next() {
            let names = first_bp.genome_names_ordered();
            for (gid, name) in names.iter().enumerate() {
                let len = first_bp.genome_seq_len(gid as u32).unwrap_or(0) as u64;
                all_genomes.push((name.clone(), len));
            }
        }

        // Phase 1: Map each k value separately (fully parallelized per k)
        println!(
            "\n[1/{}] Mapping each k separately...",
            self.k_values.len() + 1
        );

        // Store results: Vec<(read_index, k, KResult)>
        let mut all_results: Vec<(usize, KResult)> = Vec::new();

        for &k in &self.k_values {
            let bp = self.indexes.get(&k).unwrap();
            println!(
                "  k={}: mapping {} reads with {} threads...",
                k, total, threads
            );
            let pb = indicatif::ProgressBar::new(total as u64);
            pb.set_style(
                indicatif::ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} k={k}")
                    .unwrap()
                    .progress_chars("#>-"),
            );

            // Map all reads on this k in parallel with progress
            let k_results: Vec<_> = reads_vec
                .par_iter()
                .enumerate()
                .filter_map(|(idx, (_name, seq, _qual))| {
                    pb.inc(1);
                    report_atomic_progress(pb.position(), total as u64);
                    let results = bp.map_read(seq, self.context_window);
                    if let Some(best) = results.first() {
                        let genome_name = bp.genome_name(best.genome_id).unwrap_or("?").to_string();
                        Some((
                            idx,
                            KResult {
                                k,
                                genome_id: best.genome_id,
                                genome_name,
                                score: best.score,
                                rarity: best.rarity,
                                cigar: best.cigar.clone(),
                                position: best.position,
                                is_reverse: best.is_reverse,
                            },
                        ))
                    } else {
                        None
                    }
                })
                .collect();

            pb.finish_with_message(format!("  k={}: {} reads mapped", k, k_results.len()));
            all_results.extend(k_results);
        }

        // Phase 2: Combine results per read using consensus voting
        println!(
            "\n[2/{}] Combining results with {}...",
            self.k_values.len() + 1,
            match self.strategy {
                ConsensusStrategy::Majority => "majority",
                ConsensusStrategy::WeightedScore => "weighted_score",
                ConsensusStrategy::BestScore => "best_score",
                ConsensusStrategy::BaseScore => "base_score",
            }
        );

        // Group results by read index
        let mut read_results: HashMap<usize, Vec<KResult>> = HashMap::new();
        for (read_idx, kr) in all_results {
            read_results.entry(read_idx).or_default().push(kr);
        }

        // Write SAM output
        let file = File::create(output_path)
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        for (name, len) in &all_genomes {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, len)
                .map_err(|e| format!("Failed to write SQ line: {}", e))?;
        }

        let pb = indicatif::ProgressBar::new(total as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} combining")
                .unwrap()
                .progress_chars("#>-"),
        );

        let mut mapped = 0usize;

        for (read_idx, k_results) in read_results.iter().take(total) {
            let (name, seq, _qual) = &reads_vec[*read_idx];

            // Apply min_k_mappings filter
            if k_results.len() < self.min_k_mappings {
                pb.inc(1);
                report_atomic_progress(pb.position(), total as u64);
                continue;
            }

            // Apply min_score filter
            let filtered: Vec<_> = if self.min_score > 0.0 {
                k_results
                    .iter()
                    .filter(|r| r.score >= self.min_score)
                    .cloned()
                    .collect()
            } else {
                k_results.clone()
            };

            if filtered.is_empty() {
                pb.inc(1);
                report_atomic_progress(pb.position(), total as u64);
                continue;
            }

            // Consensus voting
            let candidates = match self.strategy {
                ConsensusStrategy::BestScore => {
                    let best_idx = filtered
                        .iter()
                        .enumerate()
                        .max_by(|a, b| a.1.score.partial_cmp(&b.1.score).unwrap())
                        .map(|(i, _)| i);
                    if let Some(bi) = best_idx {
                        let kr = &filtered[bi];
                        vec![ConsensusResult {
                            genome_id: kr.genome_id,
                            genome_name: kr.genome_name.clone(),
                            vote_count: 1,
                            k_count: filtered.len(),
                            consensus_score: kr.score,
                            k_results: vec![kr.clone()],
                        }]
                    } else {
                        Vec::new()
                    }
                }
                ConsensusStrategy::Majority => self.vote_majority_multi(&filtered),
                ConsensusStrategy::WeightedScore => self.vote_weighted_multi(&filtered),
                ConsensusStrategy::BaseScore => self.vote_base_score(&filtered),
            };

            let limit = if self.top_n > 1 { self.top_n } else { 1 };
            for (i, cr) in candidates.iter().take(limit).enumerate() {
                self.write_consensus_sam_line(&mut writer, name, seq, cr, i == 0)?;
                mapped += 1;
            }

            pb.inc(1);
            report_atomic_progress(pb.position(), total as u64);
        }

        pb.finish();
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output: {}", e))?;

        Ok((mapped, total))
    }

    fn compute_nm_from_cigar(&self, cigar: &str) -> u32 {
        let mut nm = 0u32;
        let mut num = String::new();
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                num.push(ch);
            } else {
                if !num.is_empty() {
                    let len: u32 = num.parse().unwrap_or(0);
                    if "IDX".contains(ch) {
                        nm += len;
                    }
                    num.clear();
                }
            }
        }
        nm
    }
}
