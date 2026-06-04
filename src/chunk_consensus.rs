use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use rayon::prelude::*;

use crate::fastq::ReadsFormat;
use crate::report_atomic_progress;
use crate::BitPop;

/// Per-config mapping result for a single read.
#[derive(Debug, Clone)]
pub struct ChunkResult {
    pub chunk_pct: f64,
    pub genome_id: u32,
    pub genome_name: String,
    pub score: f64,
    pub rarity: f64,
    pub cigar: String,
    pub position: u64,
    pub is_reverse: bool,
}

/// Consensus result for a single read across multiple chunk-% configs.
#[derive(Debug, Clone)]
pub struct ChunkConsensusResult {
    pub genome_id: u32,
    pub genome_name: String,
    pub vote_count: usize,
    pub config_count: usize,
    pub consensus_score: f64,
    pub chunk_results: Vec<ChunkResult>,
}

/// Multi chunk-% consensus mapper.
///
/// Loads the same index N times, each with a different chunk_pct configuration.
/// For each read, maps against all N configs and requires agreement from
/// at least `min_agreement` configs to accept the mapping.
pub struct MultiChunkConsensus {
    pub bp_instances: Vec<BitPop>,
    pub chunk_pcts: Vec<f64>,
    pub strategy: ConsensusStrategy,
    pub min_agreement: usize,
    pub min_score: f64,
    pub context_window: usize,
    pub chunk_min: usize,
    pub chunk_max: usize,
    pub top_n: usize,
}

/// Consensus voting strategy.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ConsensusStrategy {
    #[default]
    Majority,
    WeightedScore,
    /// Raw score sum (no chunk-weight), like JNI Android
    BaseScore,
}

/// Load a single BitPop instance from file with the given chunk_pct.
fn load_instance(
    path: &str,
    chunk_pct: f64,
    chunk_min: usize,
    chunk_max: usize,
    anchor_min_score: f64,
) -> Result<BitPop, String> {
    let mut bp = BitPop::deserialize_from_file(path)
        .map_err(|e| format!("Failed to load index {}: {}", path, e))?;
    bp.set_chunk_pct(chunk_pct);
    bp.set_chunk_min(chunk_min);
    bp.set_chunk_max(chunk_max);
    bp.set_chunk_anchor_min_score(anchor_min_score);
    Ok(bp)
}

impl MultiChunkConsensus {
    /// Load N instances of the same index, each with a different chunk_pct.
    pub fn from_path(
        index_path: &str,
        chunk_pcts: &[f64],
        chunk_min: usize,
        chunk_max: usize,
        min_score: f64,
        anchor_min_score: f64,
    ) -> Result<Self, String> {
        let mut bp_instances = Vec::new();

        for (i, &pct) in chunk_pcts.iter().enumerate() {
            println!(
                "  Loading config {} (chunk_pct={:.0}%): {}",
                i + 1,
                pct * 100.0,
                index_path
            );
            let bp = load_instance(index_path, pct, chunk_min, chunk_max, anchor_min_score)?;
            bp_instances.push(bp);
        }

        Ok(Self {
            bp_instances,
            chunk_pcts: chunk_pcts.to_vec(),
            strategy: ConsensusStrategy::default(),
            min_agreement: (chunk_pcts.len() / 2) + 1,
            min_score,
            context_window: 10,
            chunk_min,
            chunk_max,
            top_n: 1,
        })
    }

    /// Map a single read across all chunk-% configs and return consensus results.
    pub fn map_read(&self, read_seq: &str) -> Vec<ChunkConsensusResult> {
        let mut all_chunk_results: Vec<ChunkResult> = Vec::new();

        for (i, bp) in self.bp_instances.iter().enumerate() {
            let pct = self.chunk_pcts[i];

            let results = if pct > 0.0 && read_seq.len() >= bp.k() * 3 {
                bp.map_read_with_chunking(read_seq, self.context_window)
            } else {
                bp.map_read(read_seq, self.context_window)
            };

            if let Some(best) = results.first() {
                let genome_name = bp.genome_name(best.genome_id).unwrap_or("?").to_string();
                all_chunk_results.push(ChunkResult {
                    chunk_pct: pct,
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

        if all_chunk_results.is_empty() {
            return Vec::new();
        }

        // Require at least min_agreement configs to find a mapping
        if all_chunk_results.len() < self.min_agreement {
            return Vec::new();
        }

        all_chunk_results.retain(|r| r.score >= self.min_score);

        if all_chunk_results.is_empty() {
            return Vec::new();
        }

        let candidates = match self.strategy {
            ConsensusStrategy::Majority => self.vote_majority(&all_chunk_results),
            ConsensusStrategy::WeightedScore => self.vote_weighted(&all_chunk_results),
            ConsensusStrategy::BaseScore => self.vote_base_score(&all_chunk_results),
        };

        // Filter: only keep candidates that meet min_agreement threshold
        let filtered: Vec<_> = candidates
            .into_iter()
            .filter(|cr| cr.vote_count >= self.min_agreement)
            .collect();

        let limit = if self.top_n > 1 { self.top_n } else { 1 };
        filtered.into_iter().take(limit).collect()
    }

    fn vote_majority(&self, results: &[ChunkResult]) -> Vec<ChunkConsensusResult> {
        let mut genome_votes: HashMap<u32, (usize, f64, String, Vec<ChunkResult>)> = HashMap::new();

        for r in results {
            let entry = genome_votes.entry(r.genome_id).or_insert((
                0,
                0.0,
                r.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += r.score;
            entry.3.push(r.clone());
        }

        let mut candidates: Vec<_> = genome_votes
            .into_iter()
            .map(|(id, (count, score, name, rs))| (id, count, score, name, rs))
            .collect();
        candidates.sort_by(|a, b| {
            b.1.cmp(&a.1)
                .then_with(|| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
        });

        candidates
            .into_iter()
            .map(|(id, count, total_score, name, rs)| ChunkConsensusResult {
                genome_id: id,
                genome_name: name,
                vote_count: count,
                config_count: results.len(),
                consensus_score: total_score / results.len() as f64,
                chunk_results: rs,
            })
            .collect()
    }

    fn vote_weighted(&self, results: &[ChunkResult]) -> Vec<ChunkConsensusResult> {
        // Weight: lower chunk_pct → higher weight (smaller chunks = more sensitive)
        // Weight = 1.0 / chunk_pct for pct > 0, else 1.0
        let mut genome_scores: HashMap<u32, (usize, f64, String, Vec<ChunkResult>)> =
            HashMap::new();

        for r in results {
            let weight = if r.chunk_pct > 0.0 {
                1.0 / r.chunk_pct
            } else {
                1.0
            };
            let weighted = r.score * weight;
            let entry = genome_scores.entry(r.genome_id).or_insert((
                0,
                0.0,
                r.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += weighted;
            entry.3.push(r.clone());
        }

        let mut candidates: Vec<_> = genome_scores
            .into_iter()
            .map(|(id, (count, score, name, rs))| (id, count, score, name, rs))
            .collect();
        candidates.sort_by(|a, b| {
            a.2.partial_cmp(&b.2)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.1.cmp(&a.1))
        });

        candidates
            .into_iter()
            .map(|(id, count, total_score, name, rs)| ChunkConsensusResult {
                genome_id: id,
                genome_name: name,
                vote_count: count,
                config_count: results.len(),
                consensus_score: total_score / results.len() as f64,
                chunk_results: rs,
            })
            .collect()
    }

    fn vote_base_score(&self, results: &[ChunkResult]) -> Vec<ChunkConsensusResult> {
        let mut genome_scores: HashMap<u32, (usize, f64, String, Vec<ChunkResult>)> =
            HashMap::new();

        for r in results {
            let entry = genome_scores.entry(r.genome_id).or_insert((
                0,
                0.0,
                r.genome_name.clone(),
                Vec::new(),
            ));
            entry.0 += 1;
            entry.1 += r.score;
            entry.3.push(r.clone());
        }

        let mut candidates: Vec<_> = genome_scores
            .into_iter()
            .map(|(id, (count, score, name, rs))| (id, count, score, name, rs))
            .collect();
        candidates.sort_by(|a, b| {
            b.2.partial_cmp(&a.2)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.1.cmp(&a.1))
        });

        candidates
            .into_iter()
            .map(|(id, count, total_score, name, rs)| ChunkConsensusResult {
                genome_id: id,
                genome_name: name,
                vote_count: count,
                config_count: results.len(),
                consensus_score: total_score / results.len() as f64,
                chunk_results: rs,
            })
            .collect()
    }

    /// Map all reads and write SAM output.
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

        // Collect genome names from first instance
        let mut all_genomes: Vec<(String, u64)> = Vec::new();
        if let Some(first_bp) = self.bp_instances.first() {
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

        println!("\n[1/2] Mapping reads with chunk-consensus...");
        println!(
            "  Configs: {:?}",
            self.chunk_pcts
                .iter()
                .map(|p| format!("{:.0}%", p * 100.0))
                .collect::<Vec<_>>()
        );
        println!(
            "  Min agreement: {}/{}",
            self.min_agreement,
            self.bp_instances.len()
        );

        let pb = indicatif::ProgressBar::new(total as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
                .progress_chars("#>-"),
        );

        let reads_vec: Vec<_> = match &reads {
            ReadsFormat::Fastq(r) => r
                .iter()
                .map(|(n, s, q)| (n.clone(), s.clone(), q.to_vec()))
                .collect(),
            ReadsFormat::Fasta(r) => r
                .iter()
                .map(|(n, s)| (n.clone(), s.clone(), Vec::new()))
                .collect(),
        };

        let mut mapped = 0usize;

        if threads > 1 {
            let results: Vec<_> = reads_vec
                .into_par_iter()
                .map(|(name, seq, _qual)| {
                    let candidates = self.map_read(&seq);
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
                        self.write_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                        mapped += 1;
                    }
                }
                pb.inc(1);
                report_atomic_progress(pb.position(), total as u64);
            }
        } else {
            for (name, seq, _qual) in &reads_vec {
                let candidates = self.map_read(seq);
                if !candidates.is_empty() {
                    for (i, cr) in candidates.iter().enumerate() {
                        self.write_sam_line(&mut writer, name, seq, cr, i == 0)?;
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

    /// Stream version: process reads in chunks to limit memory usage.
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

        // Collect genome names from first instance
        let mut all_genomes: Vec<(String, u64)> = Vec::new();
        if let Some(first_bp) = self.bp_instances.first() {
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

        println!("\n[1/2] Mapping reads with chunk-consensus (streaming)...");

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
                        let candidates = self.map_read(&seq);
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
                            self.write_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                            mapped += 1;
                        }
                    }
                    pb.inc(1);
                    report_atomic_progress(pb.position(), total as u64);
                }
            } else {
                for (name, seq, _qual) in chunk {
                    let candidates = self.map_read(&seq);
                    if !candidates.is_empty() {
                        for (i, cr) in candidates.iter().enumerate() {
                            self.write_sam_line(&mut writer, &name, &seq, cr, i == 0)?;
                            mapped += 1;
                        }
                    }
                    pb.inc(1);
                    report_atomic_progress(pb.position(), total as u64);
                }
            }

            println!("  Chunk {} done: {} reads", chunk_num, chunk_len);
        }

        pb.finish();
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output: {}", e))?;

        Ok((mapped, total))
    }

    fn write_sam_line(
        &self,
        writer: &mut BufWriter<File>,
        read_name: &str,
        read_seq: &str,
        cr: &ChunkConsensusResult,
        is_primary: bool,
    ) -> Result<(), String> {
        let mut flag: u16 = if cr.chunk_results.is_empty() {
            0x4
        } else if cr.chunk_results[0].is_reverse {
            0x10
        } else {
            0
        };

        if !is_primary {
            flag |= 0x800;
        }

        let mapq = (cr.consensus_score * 60.0) as u16;
        let pos = if cr.chunk_results.is_empty() {
            0
        } else {
            (cr.chunk_results[0].position + 1) as u32
        };
        let cigar = if cr.chunk_results.is_empty() {
            "*".to_string()
        } else {
            cr.chunk_results[0].cigar.clone()
        };

        let nm = if cr.chunk_results.is_empty() {
            0u32
        } else {
            self.compute_nm_from_cigar(&cr.chunk_results[0].cigar)
        };

        // Per-config tags: CP1:Z:genome/score, CP10:Z:genome/score, etc.
        let mut cp_tags = String::new();
        for cr_item in &cr.chunk_results {
            let pct_label = (cr_item.chunk_pct * 100.0) as u32;
            cp_tags.push_str(&format!(
                "\tCP{}:Z:{}/{:.4}",
                pct_label, cr_item.genome_name, cr_item.score
            ));
        }

        let cv_tag = format!("\tCV:Z:{}", cr.genome_name);
        let cc_tag = format!("\tCC:i:{}", cr.config_count);
        let vc_tag = format!("\tVC:i:{}", cr.vote_count);
        let as_tag = format!("\tAS:f:{:.4}", cr.consensus_score);
        let xs_tag = if is_primary {
            String::new()
        } else {
            format!("\tXS:f:{:.4}", cr.consensus_score)
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tNM:i:{}{}{}{}{}{}{}",
            read_name,
            flag,
            cr.genome_name,
            pos,
            mapq,
            cigar,
            read_seq,
            nm,
            cp_tags,
            cv_tag,
            cc_tag,
            vc_tag,
            as_tag,
            xs_tag,
        )
        .map_err(|e| format!("Failed to write SAM line: {}", e))?;
        Ok(())
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
