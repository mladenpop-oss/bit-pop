use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult};

/// Parse a FASTQ file and return read names, sequences, and quality scores.
pub fn parse_fastq(path: &str) -> IoResult<Vec<(String, String, Vec<u8>)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();
    let mut lines = reader.lines();

    while let Some(header_line) = lines.next() {
        let header = match header_line {
            Ok(h) => h,
            Err(_) => break,
        };

        if !header.starts_with('@') {
            break;
        }

        let read_name = header[1..].trim().to_string();

        let seq_line = match lines.next() {
            Some(Ok(s)) => s.trim().to_string(),
            _ => break,
        };

        let _plus_line = match lines.next() {
            Some(Ok(p)) => p,
            _ => break,
        };

        let qual_line = match lines.next() {
            Some(Ok(q)) => q,
            _ => break,
        };

        // Phred+33 encoding: quality byte = ascii_code - 33
        let quality: Vec<u8> = qual_line
            .bytes()
            .map(|b| b.saturating_sub(33))
            .collect();

        reads.push((read_name, seq_line, quality));
    }

    Ok(reads)
}

/// Parse a FASTA file and return read names and sequences (no quality scores).
pub fn parse_fasta(path: &str) -> IoResult<Vec<(String, String)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();
    let mut header = String::new();
    let mut sequence = String::new();

    for line in reader.lines() {
        match line {
            Ok(l) => {
                if l.starts_with('>') {
                    if !header.is_empty() {
                        reads.push((header.clone(), sequence.clone()));
                    }
                    header = l[1..].trim().to_string();
                    sequence = String::new();
                } else if !l.trim().is_empty() {
                    sequence.push_str(l.trim());
                }
            }
            Err(_) => break,
        }
    }

    if !header.is_empty() {
        reads.push((header, sequence));
    }

    Ok(reads)
}

/// Auto-detect file format and parse accordingly.
pub fn parse_reads(path: &str) -> IoResult<ReadsFormat> {
    let ext = path.to_lowercase();
    
    if ext.ends_with(".fastq") || ext.ends_with(".fq") {
        let reads = parse_fastq(path)?;
        Ok(ReadsFormat::Fastq(reads))
    } else {
        let reads = parse_fasta(path)?;
        Ok(ReadsFormat::Fasta(reads))
    }
}

/// Parsed reads in either FASTA or FASTQ format.
pub enum ReadsFormat {
    Fasta(Vec<(String, String)>),
    Fastq(Vec<(String, String, Vec<u8>)>),
}

impl ReadsFormat {
    pub fn count(&self) -> usize {
        match self {
            ReadsFormat::Fasta(reads) => reads.len(),
            ReadsFormat::Fastq(reads) => reads.len(),
        }
    }

    pub fn iter_fasta(&self) -> Box<dyn Iterator<Item = (&str, &str)> + '_> {
        match self {
            ReadsFormat::Fasta(reads) => Box::new(reads.iter().map(|(name, seq)| (name.as_str(), seq.as_str()))),
            ReadsFormat::Fastq(reads) => Box::new(reads.iter().map(|(name, seq, _)| (name.as_str(), seq.as_str()))),
        }
    }

    pub fn iter_fastq(&self) -> Option<Box<dyn Iterator<Item = (&str, &str, &[u8])> + '_>> {
        match self {
            ReadsFormat::Fastq(reads) => Some(Box::new(reads.iter().map(|(name, seq, qual)| (name.as_str(), seq.as_str(), qual.as_slice())))),
            ReadsFormat::Fasta(_) => None,
        }
    }

    pub fn has_quality(&self) -> bool {
        matches!(self, ReadsFormat::Fastq(_))
    }
}

/// Filter reads by minimum average quality score.
/// Returns indices of reads that pass the filter.
pub fn filter_by_quality(reads: &[(String, String, Vec<u8>)], min_avg_qual: u8) -> Vec<usize> {
    reads
        .iter()
        .enumerate()
        .filter_map(|(i, (_, _, qual))| {
            if qual.is_empty() {
                return Some(i);
            }
            let avg: u32 = qual.iter().map(|&q| q as u32).sum();
            let avg = avg / qual.len() as u32;
            if avg >= min_avg_qual as u32 {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

/// Normalize a read name by stripping /1, /2, :1, :2 suffixes and .fq/.fastq extensions.
fn normalize_read_name(name: &str) -> String {
    let name = name.trim();
    // Strip common paired-end suffixes: "read/1", "read/2", "read:1", "read:2"
    let name = if name.ends_with("/1") || name.ends_with("/2") {
        &name[..name.len() - 2]
    } else if name.ends_with(":1") || name.ends_with(":2") {
        &name[..name.len() - 2]
    } else {
        name
    };
    name.to_string()
}

/// Check if a read name looks like it's from the first or second read in a pair.
fn is_read_1(name: &str) -> bool {
    name.ends_with("/1") || name.ends_with(":1")
}

fn is_read_2(name: &str) -> bool {
    name.ends_with("/2") || name.ends_with(":2")
}

/// Parse a FASTQ file for paired-end reads and return paired read pairs.
/// Matches R1/R2 by normalized name (strips /1, /2, :1, :2 suffixes).
/// Returns Vec of (normalized_name, read1_seq, read1_qual, read2_seq, read2_qual).
pub fn parse_paired_fastq(path1: &str, path2: &str) -> IoResult<Vec<(String, String, Vec<u8>, String, Vec<u8>)>> {
    let reads1 = parse_fastq(path1)?;
    let reads2 = parse_fastq(path2)?;

    // Build lookup from normalized name -> read data
    let mut map1: HashMap<String, (String, Vec<u8>)> = HashMap::new();
    for (name, seq, qual) in &reads1 {
        let norm = normalize_read_name(name);
        if is_read_2(name) {
            continue; // skip /2 entries from file1
        }
        map1.insert(norm, (seq.clone(), qual.clone()));
    }

    let mut pairs: Vec<(String, String, Vec<u8>, String, Vec<u8>)> = Vec::new();
    for (name, seq, qual) in &reads2 {
        if !is_read_2(name) {
            continue; // skip /1 entries from file2
        }
        let norm = normalize_read_name(name);
        if let Some((seq1, qual1)) = map1.get(&norm) {
            pairs.push((norm.clone(), seq1.clone(), qual1.clone(), seq.clone(), qual.clone()));
        }
    }

    // Sort by name for deterministic output
    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(pairs)
}

/// Parse single FASTQ file and return paired-end reads assuming interleaved format.
/// In interleaved FASTQ, read1 and read2 alternate: R1, R2, R1, R2, ...
pub fn parse_interleaved_paired_fastq(path: &str) -> IoResult<Vec<(String, String, Vec<u8>, String, Vec<u8>)>> {
    let reads = parse_fastq(path)?;

    if reads.len() % 2 != 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Interleaved paired-end FASTQ must have even number of reads",
        ));
    }

    let mut pairs: Vec<(String, String, Vec<u8>, String, Vec<u8>)> = Vec::new();
    let mut i = 0;
    while i + 1 < reads.len() {
        let (name1, seq1, qual1) = &reads[i];
        let (name2, seq2, qual2) = &reads[i + 1];
        let norm1 = normalize_read_name(name1);
        let norm2 = normalize_read_name(name2);

        if norm1 == norm2 {
            pairs.push((norm1.clone(), seq1.clone(), qual1.clone(), seq2.clone(), qual2.clone()));
        } else {
            // Different names — treat as unmatched, skip or pair by position
            pairs.push((format!("unmatched_{}", i), seq1.clone(), qual1.clone(), seq2.clone(), qual2.clone()));
        }
        i += 2;
    }

    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(pairs)
}

use std::collections::HashMap;

/// Compute Phred-scaled mismatch penalty for a single base pair.
/// Higher quality mismatches are penalized more.
pub fn phred_mismatch_penalty(read_qual: u8, genome_qual: u8) -> f64 {
    // Average of read and genome quality (genome qual assumed from reference)
    let avg_qual = ((read_qual as f64) + (genome_qual as f64)) / 2.0;
    
    // Scale: Q10 = -1, Q20 = -2, Q30 = -3, Q40 = -4
    (avg_qual as f64 / 10.0).min(5.0) * -1.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_test_fastq() -> (String, Vec<(String, String, Vec<u8>)>) {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_reads_{}_{}.fastq", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "@read1").unwrap();
        writeln!(file, "ACGTACGTACGT").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIIIIIIIIIII").unwrap();
        writeln!(file, "@read2").unwrap();
        writeln!(file, "TTTTGGGGAAAA").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "!!!!!!!!!!!!").unwrap();
        
        let reads = parse_fastq(path.to_str().unwrap()).unwrap();
        (path.to_str().unwrap().to_string(), reads)
    }

    #[test]
    fn test_parse_fastq() {
        let (_path, reads) = make_test_fastq();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].0, "read1");
        assert_eq!(reads[0].1, "ACGTACGTACGT");
        // I = ascii 73, quality = 73 - 33 = 40
        assert_eq!(reads[0].2.len(), 12);
        assert_eq!(reads[0].2[0], 40);

        // ! = ascii 33, quality = 33 - 33 = 0
        assert_eq!(reads[1].2[0], 0);
    }

    #[test]
    fn test_parse_fasta() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_fasta_{}_{}.fasta", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, ">read1 description").unwrap();
        writeln!(file, "ACGTACGT").unwrap();
        writeln!(file, ">read2").unwrap();
        writeln!(file, "TTTTGGGG").unwrap();
        
        let reads = parse_fasta(path.to_str().unwrap()).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].0, "read1 description");
        assert_eq!(reads[0].1, "ACGTACGT");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_filter_by_quality() {
        let reads = vec![
            ("r1".to_string(), "ACGT".to_string(), vec![30, 30, 30, 30]), // avg 30
            ("r2".to_string(), "ACGT".to_string(), vec![10, 10, 10, 10]), // avg 10
            ("r3".to_string(), "ACGT".to_string(), vec![25, 25, 25, 25]), // avg 25
        ];

        let passed = filter_by_quality(&reads, 20);
        assert_eq!(passed.len(), 2);
        assert!(passed.contains(&0)); // r1 passes
        assert!(!passed.contains(&1)); // r2 fails
        assert!(passed.contains(&2)); // r3 passes
    }

    #[test]
    fn test_phred_mismatch_penalty() {
        // Q20 read, Q20 genome → penalty ≈ -2.0
        let penalty = phred_mismatch_penalty(20, 20);
        assert!((penalty - (-2.0)).abs() < 0.01);

        // Q30 read, Q30 genome → penalty ≈ -3.0
        let penalty = phred_mismatch_penalty(30, 30);
        assert!((penalty - (-3.0)).abs() < 0.01);

        // Q10 read, Q10 genome → penalty ≈ -1.0
        let penalty = phred_mismatch_penalty(10, 10);
        assert!((penalty - (-1.0)).abs() < 0.01);
    }

    #[test]
    fn test_reads_format_iterators() {
        let (_path, fastq_reads) = make_test_fastq();
        let format = ReadsFormat::Fastq(fastq_reads);

        assert!(format.has_quality());
        assert_eq!(format.count(), 2);

        // Test FASTQ iterator
        let fastq_iter = format.iter_fastq().unwrap();
        for (name, seq, qual) in fastq_iter {
            assert!(!name.is_empty());
            assert!(!seq.is_empty());
            assert!(!qual.is_empty());
        }

        // Test FASTA iterator (also works on Fastq data, ignores quality)
        let fasta_iter = format.iter_fasta();
        let count: usize = fasta_iter.count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_normalize_read_name() {
        assert_eq!(normalize_read_name("read1/1"), "read1");
        assert_eq!(normalize_read_name("read1/2"), "read1");
        assert_eq!(normalize_read_name("read1:1"), "read1");
        assert_eq!(normalize_read_name("read1:2"), "read1");
        assert_eq!(normalize_read_name("read1"), "read1");
    }

    #[test]
    fn test_is_read_1_and_2() {
        assert!(is_read_1("read/1"));
        assert!(!is_read_1("read/2"));
        assert!(is_read_2("read/2"));
        assert!(!is_read_2("read/1"));
    }

    #[test]
    fn test_parse_paired_fastq() {
        let dir = std::env::temp_dir();
        let path1 = dir.join(format!("paired_r1_{}_{}.fastq", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        let path2 = dir.join(format!("paired_r2_{}_{}.fastq", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));

        {
            let mut f = std::fs::File::create(&path1).unwrap();
            writeln!(f, "@read1/1").unwrap();
            writeln!(f, "ACGTACGT").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "IIIIIIII").unwrap();
            writeln!(f, "@read2/1").unwrap();
            writeln!(f, "TTTTGGGG").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "!!!!!!!!").unwrap();
        }

        {
            let mut f = std::fs::File::create(&path2).unwrap();
            writeln!(f, "@read1/2").unwrap();
            writeln!(f, "ACGTACGT").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "IIIIIIII").unwrap();
            writeln!(f, "@read2/2").unwrap();
            writeln!(f, "TTTTGGGG").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "!!!!!!!!").unwrap();
        }

        let pairs = parse_paired_fastq(path1.to_str().unwrap(), path2.to_str().unwrap()).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "read1");
        assert_eq!(pairs[0].1, "ACGTACGT"); // R1 seq
        assert_eq!(pairs[0].3, "ACGTACGT"); // R2 seq

        let _ = std::fs::remove_file(&path1);
        let _ = std::fs::remove_file(&path2);
    }
}
