use jni::objects::{JClass, JString};
use jni::sys::{jint, jstring};
use jni::JNIEnv;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::{fastq, BitPop};

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_mapReads(
    mut env: JNIEnv,
    _: JClass,
    index_path: JString,
    reads_path: JString,
    output_path: JString,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let idx: String = env.get_string(&index_path)?.into();
        let rds: String = env.get_string(&reads_path)?.into();
        let out: String = env.get_string(&output_path)?.into();

        let bp = BitPop::deserialize_from_file(&idx)?;
        let reads = fastq::parse_fastq(&rds)?;

        let mut results = Vec::new();
        for (name, seq, _) in reads {
            let mapped = bp.map_read(&seq, 0);
            if let Some(best) = mapped.first() {
                let genome_name = bp.genome_name(best.genome_id).unwrap_or("unknown");
                results.push(format!("{}\t{}\t{:.2}", name, genome_name, best.score));
            }
        }

        std::fs::write(&out, results.join("\n") + "\n")?;
        Ok(format!("Mapped {} reads to {}", results.len(), out))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_buildIndex(
    mut env: JNIEnv,
    _: JClass,
    fasta_path: JString,
    output_path: JString,
    k_mer_size: jint,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let fasta: String = env.get_string(&fasta_path)?.into();
        let output: String = env.get_string(&output_path)?.into();

        let mut bp = BitPop::new(k_mer_size as usize);
        let genomes = fastq::parse_fasta(&fasta)?;

        for (name, seq) in genomes {
            bp.add_genome(&name, &seq);
        }

        bp.build();
        bp.serialize_to_file(&output)?;

        Ok(format!("Index saved to {}", output))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_getGenomeNames(
    mut env: JNIEnv,
    _: JClass,
    index_path: JString,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let idx: String = env.get_string(&index_path)?.into();
        let bp = BitPop::deserialize_from_file(&idx)?;
        let names: Vec<String> = bp.genome_names_ordered();
        Ok(names.join(","))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_fastConMap(
    mut env: JNIEnv,
    _: JClass,
    index_paths: JString,
    reads_path: JString,
    output_path: JString,
    chunk_pct: jint,
    chunk_min: jint,
    chunk_max: jint,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let idx_str: String = env.get_string(&index_paths)?.into();
        let rds: String = env.get_string(&reads_path)?.into();
        let out: String = env.get_string(&output_path)?.into();
        let chunk_pct_val = chunk_pct as f64;
        let chunk_min_val = chunk_min as usize;
        let chunk_max_val = chunk_max as usize;

        let paths: Vec<String> = idx_str.split(',').map(|s| s.trim().to_string()).collect();
        if paths.len() < 2 {
            return Err("Need at least 2 index paths".into());
        }

        let reads = fastq::parse_fastq(&rds)?;
        let total = reads.len();

        // Phase 1: Map each index, save SAM
        let mut sam_results: Vec<Vec<(String, String, f64)>> = Vec::new();

        for (i, p) in paths.iter().enumerate() {
            let bp = BitPop::deserialize_from_file(p)?;
            let mut results: Vec<(String, String, f64)> = Vec::new();

            for (name, seq, _) in &reads {
                let chunks = if chunk_pct_val > 0.0 {
                    make_chunks(seq, chunk_pct_val, chunk_min_val, chunk_max_val)
                } else {
                    vec![seq.as_str()]
                };

                let mut best_score = 0.0f64;
                let mut best_genome_id: Option<u32> = None;

                for chunk in &chunks {
                    let mapped = bp.map_read(chunk, 4);
                    if let Some(r) = mapped.first() {
                        if best_genome_id.is_none() || r.score > best_score {
                            best_score = r.score;
                            best_genome_id = Some(r.genome_id);
                        }
                    }
                }

                if let Some(gid) = best_genome_id {
                    let genome = bp.genome_name(gid).unwrap_or("?");
                    results.push((name.clone(), genome.to_string(), best_score));
                }
            }

            sam_results.push(results);
        }

        // Phase 2: Consensus
        let mut read_votes: HashMap<String, Vec<(String, f64)>> = HashMap::new();
        for results in &sam_results {
            for (name, genome, score) in results {
                read_votes
                    .entry(name.clone())
                    .or_default()
                    .push((genome.clone(), *score));
            }
        }

        // Phase 3: Write output - all reads, consensus where both indexes agree
        let file = File::create(&out)?;
        let mut writer = BufWriter::new(file);
        let mut mapped = 0usize;

        for (read_name, votes) in &read_votes {
            if votes.is_empty() {
                continue;
            }

            let mut genome_scores: HashMap<String, (usize, f64)> = HashMap::new();
            for (genome, score) in votes {
                let entry = genome_scores.entry(genome.clone()).or_insert((0, 0.0));
                entry.0 += 1;
                entry.1 += score;
            }

            let mut candidates: Vec<_> = genome_scores.into_iter().collect();
            candidates.sort_by(|a, b| {
                b.1 .1
                    .partial_cmp(&a.1 .1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            let limit = std::cmp::min(2, candidates.len());
            for (genome, (vote_count, total_score)) in candidates.iter().take(limit) {
                let avg_score = total_score / votes.len() as f64;
                writeln!(
                    writer,
                    "{}\t{}\t{:.4}\t{}",
                    read_name, genome, avg_score, vote_count
                )?;
                if *vote_count == candidates[0].1 .0 {
                    mapped += 1;
                }
            }
        }

        writer.flush()?;
        Ok(format!("FastCon: {} / {} reads mapped", mapped, total))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

fn make_chunks(seq: &str, chunk_pct: f64, chunk_min: usize, chunk_max: usize) -> Vec<&str> {
    let len = seq.len();
    if len <= chunk_max {
        return vec![seq];
    }
    // Same logic as desktop: chunk_pct controls chunk size, not step
    // chunk_pct is in percent (1 = 1%), convert to fraction
    let chunk_size =
        (len as f64 * chunk_pct / 100.0).clamp(chunk_min as f64, chunk_max as f64) as usize;
    let overlap = (chunk_size as f64 * 0.2) as usize;
    let step = chunk_size - overlap;
    let mut chunks = Vec::new();
    let mut start = 0usize;
    while start < len {
        let end = std::cmp::min(start + chunk_size, len);
        chunks.push(&seq[start..end]);
        start += step;
        if start >= len {
            break;
        }
    }
    chunks
}
