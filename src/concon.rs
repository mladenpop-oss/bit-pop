use crate::BitPop;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Per-index mapping result for a single read.
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

/// Consensus result for a single read across multiple indexes.
#[derive(Debug, Clone)]
pub struct ConsensusResult {
    pub genome_id: u32,
    pub genome_name: String,
    pub vote_count: usize,
    pub k_count: usize,
    pub consensus_score: f64,
    pub k_results: Vec<KResult>,
}

/// Consensus voting strategy.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ConsensusStrategy {
    Majority,
    #[default]
    WeightedScore,
    BestScore,
    /// Raw score sum (no k-weight), like JNI Android
    BaseScore,
}

/// Fast consensus: runs `bit-pop map` subprocess for each index, then combines.
pub struct ConCon {
    /// Index paths (k-value read from each index file)
    pub indexes: Vec<PathBuf>,
    /// k-values loaded from index files (same order as indexes)
    pub k_values: Vec<usize>,
    pub strategy: ConsensusStrategy,
    pub min_score: f64,
    pub min_k_mappings: usize,
    pub top_n: usize,
    pub map_top_n: usize,
    pub k_weights: HashMap<usize, f64>,
    pub bit_pop_exe: PathBuf,
    pub chunk_size: usize,
    pub chunk_pct: f64,
    pub chunk_min: usize,
    pub chunk_max: usize,
    pub context_window: usize,
    pub anchor_min_score: f64,
    pub anchor_filter: bool,
}

impl ConCon {
    pub fn new(
        indexes: Vec<PathBuf>,
        bit_pop_exe: PathBuf,
        min_score: f64,
        context_window: usize,
        anchor_min_score: f64,
        anchor_filter: bool,
    ) -> Result<Self, String> {
        // Load k-values from index files
        let mut k_values = Vec::new();
        for path in &indexes {
            let bp = BitPop::deserialize_from_file(path.to_str().unwrap())
                .map_err(|e| format!("Failed to load index {}: {}", path.display(), e))?;
            k_values.push(bp.k());
        }

        let min_k = k_values.iter().min().ok_or("No indexes provided")?;

        let mut k_weights = HashMap::new();
        for &k in &k_values {
            k_weights.insert(k, k as f64 / *min_k as f64);
        }

        Ok(Self {
            indexes,
            k_values,
            strategy: ConsensusStrategy::default(),
            min_score,
            min_k_mappings: 1,
            top_n: 1,
            map_top_n: 1,
            k_weights,
            bit_pop_exe,
            chunk_size: 0,
            chunk_pct: 0.0,
            chunk_min: 20,
            chunk_max: 500,
            context_window,
            anchor_min_score,
            anchor_filter,
        })
    }

    /// Derive temp SAM path from index path: index.bitpop -> index.concon.sam
    fn temp_sam_path(index_path: &Path) -> PathBuf {
        let mut stem = index_path
            .file_stem()
            .unwrap_or(index_path.as_os_str())
            .to_string_lossy()
            .to_string();
        if let Some(ext) = index_path.extension() {
            stem.push_str(&format!(".{}", ext.to_string_lossy()));
        }
        stem.push_str(".concon.sam");
        index_path
            .parent()
            .unwrap_or(std::path::Path::new("."))
            .join(stem)
    }

    // Phase 1: Run `bit-pop map` for each index
    fn phase1_map(&self, reads_path: &Path, threads: usize) -> Result<Vec<PathBuf>, String> {
        let mut sam_paths = Vec::new();

        for index_path in &self.indexes {
            let sam_path = Self::temp_sam_path(index_path);
            println!(
                "  Mapping: {} -> {}",
                index_path.display(),
                sam_path.display()
            );

            let mut args: Vec<String> = vec![
                "map".into(),
                "-i".into(),
                index_path.to_string_lossy().to_string(),
                "-r".into(),
                reads_path.to_string_lossy().to_string(),
                "-o".into(),
                sam_path.to_string_lossy().to_string(),
                "-a".into(),
                "xor".into(),
                "--top-n".into(),
                self.map_top_n.to_string(),
                "-t".into(),
                threads.to_string(),
            ];

            if self.chunk_pct > 0.0 {
                args.push("--chunk-pct".into());
                args.push(self.chunk_pct.to_string());
            } else if self.chunk_size > 0 {
                args.push("--chunk-size".into());
                args.push(self.chunk_size.to_string());
            }

            if self.chunk_min != 20 {
                args.push("--chunk-min".into());
                args.push(self.chunk_min.to_string());
            }

            if self.chunk_max != 500 {
                args.push("--chunk-max".into());
                args.push(self.chunk_max.to_string());
            }

            if self.min_score < 0.7 {
                args.push("--min-score".into());
                args.push(self.min_score.to_string());
            }

            if self.context_window != 50 {
                args.push("--context-window".into());
                args.push(self.context_window.to_string());
            }

            if self.anchor_min_score != 0.5 {
                args.push("--anchor-min-score".into());
                args.push(self.anchor_min_score.to_string());
            }

            if self.anchor_filter {
                args.push("--anchor-filter".into());
            }

            let mut child = Command::new(&self.bit_pop_exe)
                .args(&args)
                .stdout(std::process::Stdio::inherit())
                .stderr(std::process::Stdio::inherit())
                .spawn()
                .map_err(|e| format!("Failed to run bit-pop map: {}", e))?;

            let status = child
                .wait()
                .map_err(|e| format!("Failed to wait for bit-pop map: {}", e))?;

            if !status.success() {
                return Err(format!("bit-pop map failed with status: {}", status));
            }

            sam_paths.push(sam_path);
        }

        Ok(sam_paths)
    }

    // Phase 2: Parse all temp SAMs
    fn phase2_parse(
        &self,
        sam_paths: &[PathBuf],
    ) -> Result<(HashMap<String, Vec<KResult>>, HashMap<String, u32>), String> {
        // Build name -> id from first SAM
        let mut name_to_id: HashMap<String, u32> = HashMap::new();
        if !sam_paths.is_empty() {
            name_to_id = Self::build_name_to_id(&sam_paths[0])?;
        }

        let mut read_results: HashMap<String, Vec<KResult>> = HashMap::new();

        for (i, sam_path) in sam_paths.iter().enumerate() {
            let k = self.k_values[i];
            println!("  Parsing: {}", sam_path.display());

            let file = File::open(sam_path)
                .map_err(|e| format!("Failed to open SAM {}: {}", sam_path.display(), e))?;
            let reader = BufReader::new(file);

            for line in reader.lines() {
                let line = line.map_err(|e| format!("Failed to read SAM: {}", e))?;
                if line.starts_with('@') {
                    continue;
                }

                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 11 {
                    continue;
                }

                let read_name = fields[0].to_string();
                let genome_name = fields[2].to_string();
                if genome_name == "*" {
                    continue;
                }

                let pos: u64 = fields[3].parse().unwrap_or(1) - 1;
                let cigar = fields[5].to_string();
                let flag: i32 = fields[1].parse().unwrap_or(0);
                let is_reverse = flag & 0x10 != 0;

                let mut score = 0.0;
                let mut rarity = 0.0;

                for f in &fields[11..] {
                    if let Some(val) = f.strip_prefix("AS:f:") {
                        score = val.parse().unwrap_or(0.0);
                    } else if let Some(val) = f.strip_prefix("RK:f:") {
                        rarity = val.parse().unwrap_or(0.0);
                    }
                }

                let genome_id = *name_to_id.get(&genome_name).unwrap_or(&0);

                read_results.entry(read_name).or_default().push(KResult {
                    k,
                    genome_id,
                    genome_name,
                    score,
                    rarity,
                    cigar,
                    position: pos,
                    is_reverse,
                });
            }
        }

        Ok((read_results, name_to_id))
    }

    fn build_name_to_id(sam_path: &PathBuf) -> Result<HashMap<String, u32>, String> {
        let mut name_to_id: HashMap<String, u32> = HashMap::new();
        let mut gid: u32 = 0;

        let file = File::open(sam_path)
            .map_err(|e| format!("Failed to open SAM {}: {}", sam_path.display(), e))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read SAM: {}", e))?;
            if let Some(rest) = line.strip_prefix("@SQ") {
                for field in rest.split('\t') {
                    if let Some(name) = field.strip_prefix("SN:") {
                        name_to_id.insert(name.to_string(), gid);
                        gid += 1;
                        break;
                    }
                }
            }
        }

        Ok(name_to_id)
    }

    // Phase 3: Consensus voting + write output SAM
    fn phase3_combine(
        &self,
        read_results: &HashMap<String, Vec<KResult>>,
        _name_to_id: &HashMap<String, u32>,
        output_path: &PathBuf,
        sam_paths: &[PathBuf],
    ) -> Result<(usize, usize), String> {
        // Collect genome names and lengths from first SAM
        let mut all_genomes: Vec<(String, u64)> = Vec::new();
        if !sam_paths.is_empty() {
            let file =
                File::open(&sam_paths[0]).map_err(|e| format!("Failed to open SAM: {}", e))?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line.map_err(|e| format!("Failed to read SAM: {}", e))?;
                if let Some(rest) = line.strip_prefix("@SQ") {
                    let mut name = String::new();
                    let mut len: u64 = 0;
                    for field in rest.split('\t') {
                        if let Some(n) = field.strip_prefix("SN:") {
                            name = n.to_string();
                        } else if let Some(l) = field.strip_prefix("LN:") {
                            len = l.parse().unwrap_or(0);
                        }
                    }
                    if !name.is_empty() {
                        all_genomes.push((name, len));
                    }
                }
            }
        }

        // Count total reads from first SAM (non-header lines)
        let total = {
            let file =
                File::open(&sam_paths[0]).map_err(|e| format!("Failed to open SAM: {}", e))?;
            let reader = BufReader::new(file);
            reader
                .lines()
                .filter(|l| {
                    if let Ok(l) = l {
                        !l.starts_with('@')
                    } else {
                        false
                    }
                })
                .count()
        };

        let file = File::create(output_path)
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")
            .map_err(|e| format!("Failed to write header: {}", e))?;
        for (name, len) in &all_genomes {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, len)
                .map_err(|e| format!("Failed to write SQ line: {}", e))?;
        }

        println!("  Combining results (total reads in first SAM: {})", total);

        // Collect all read names from all SAMs to get unique set
        let mut all_read_names: Vec<String> = Vec::new();
        for sam_path in sam_paths {
            let file = File::open(sam_path).map_err(|e| format!("Failed to open SAM: {}", e))?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line.map_err(|e| format!("Failed to read SAM: {}", e))?;
                if !line.starts_with('@') {
                    if let Some(name) = line.split('\t').next() {
                        all_read_names.push(name.to_string());
                    }
                }
            }
        }
        all_read_names.sort();
        all_read_names.dedup();
        let total_unique = all_read_names.len();

        let pb = indicatif::ProgressBar::new(total_unique as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} combining")
                .unwrap()
                .progress_chars("#>-"),
        );

        let mut mapped_reads = 0usize;

        // Build read sequence lookup from first SAM
        let mut read_seqs: HashMap<String, String> = HashMap::new();
        for sam_path in sam_paths {
            let file = File::open(sam_path).map_err(|e| format!("Failed to open SAM: {}", e))?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line.map_err(|e| format!("Failed to read SAM: {}", e))?;
                if !line.starts_with('@') {
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() >= 10 {
                        let name = fields[0].to_string();
                        let seq = fields[9].to_string();
                        read_seqs.entry(name).or_insert(seq);
                    }
                }
            }
        }

        for read_name in &all_read_names {
            let seq = read_seqs
                .get(read_name)
                .cloned()
                .unwrap_or_else(|| "*".to_string());

            let k_results = match read_results.get(read_name) {
                Some(r) => r,
                None => {
                    // Read not found in any index - write as unmapped
                    self.write_unmapped_line(&mut writer, read_name, &seq)?;
                    pb.inc(1);
                    continue;
                }
            };

            if k_results.len() < self.min_k_mappings {
                // Not enough mappings across indexes - write as unmapped
                self.write_unmapped_line(&mut writer, read_name, &seq)?;
                pb.inc(1);
                continue;
            }

            let filtered: Vec<KResult> = if self.min_score > 0.0 {
                k_results
                    .iter()
                    .filter(|r| r.score >= self.min_score)
                    .cloned()
                    .collect()
            } else {
                k_results.clone()
            };

            if filtered.is_empty() {
                // All mappings below min score - write as unmapped
                self.write_unmapped_line(&mut writer, read_name, &seq)?;
                pb.inc(1);
                continue;
            }

            let candidates = match self.strategy {
                ConsensusStrategy::BestScore => self.vote_best_score(&filtered),
                ConsensusStrategy::Majority => self.vote_majority(&filtered),
                ConsensusStrategy::WeightedScore => self.vote_weighted(&filtered),
                ConsensusStrategy::BaseScore => self.vote_base_score(&filtered),
            };

            let limit = if self.top_n > 1 { self.top_n } else { 1 };

            for (i, cr) in candidates.iter().take(limit).enumerate() {
                self.write_sam_line(&mut writer, read_name, &seq, cr, i == 0)?;
            }
            mapped_reads += 1;

            pb.inc(1);
        }

        pb.finish();
        writer
            .flush()
            .map_err(|e| format!("Failed to flush output: {}", e))?;

        Ok((mapped_reads, total_unique))
    }

    fn vote_weighted(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
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

    fn vote_majority(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
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
            b.1.cmp(&a.1)
                .then_with(|| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
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

    fn vote_best_score(&self, k_results: &[KResult]) -> Vec<ConsensusResult> {
        let best = k_results
            .iter()
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
            .unwrap();

        vec![ConsensusResult {
            genome_id: best.genome_id,
            genome_name: best.genome_name.clone(),
            vote_count: 1,
            k_count: k_results.len(),
            consensus_score: best.score,
            k_results: vec![best.clone()],
        }]
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

    fn write_sam_line(
        &self,
        writer: &mut BufWriter<File>,
        read_name: &str,
        read_seq: &str,
        cr: &ConsensusResult,
        is_primary: bool,
    ) -> Result<(), String> {
        let best = cr.k_results.first().unwrap();

        let mut flag: u16 = if best.is_reverse { 0x10 } else { 0 };
        if !is_primary {
            flag |= 0x800;
        }

        let mapq = (cr.consensus_score * 60.0) as u16;
        let pos = (best.position + 1) as u32;
        let cigar = &best.cigar;
        let nm = self.compute_nm_from_cigar(&best.cigar);

        let rk_tags: String = cr
            .k_results
            .iter()
            .map(|kr| format!("\tRK{}:Z:{}", kr.k, kr.genome_name))
            .collect();

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
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tNM:i:{}{}{}{}{}{}{}",
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

    fn write_unmapped_line(
        &self,
        writer: &mut BufWriter<File>,
        read_name: &str,
        read_seq: &str,
    ) -> Result<(), String> {
        // Write unmapped SAM line: flag=4, genome="*", pos=0, mapq=0, cigar="*"
        writeln!(
            writer,
            "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t*\tAS:f:0.0000",
            read_name, read_seq,
        )
        .map_err(|e| format!("Failed to write unmapped SAM line: {}", e))?;
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

    /// Main entry: run map for each index, combine, write output.
    pub fn run(
        &self,
        reads_path: &Path,
        output_path: &PathBuf,
        threads: usize,
    ) -> Result<(usize, usize), String> {
        println!("Phase 1: Mapping each index with `bit-pop map`");
        println!("{}", "=".repeat(50));
        let sam_paths = self.phase1_map(reads_path, threads)?;

        println!();
        println!("Phase 2: Parsing temp SAMs");
        println!("{}", "=".repeat(50));
        let (read_results, name_to_id) = self.phase2_parse(&sam_paths)?;

        println!();
        println!("Phase 3: Writing consensus SAM");
        println!("{}", "=".repeat(50));
        let (mapped, total) =
            self.phase3_combine(&read_results, &name_to_id, output_path, &sam_paths)?;

        println!();
        println!("Temp SAM files (kept for further processing):");
        for p in &sam_paths {
            println!("  {}", p.display());
        }
        println!();
        println!("Mapped: {} / {} reads", mapped, total);
        println!("Output: {}", output_path.display());

        Ok((mapped, total))
    }
}
