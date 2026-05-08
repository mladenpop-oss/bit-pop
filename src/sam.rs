use std::fs::File;
use std::io::{self, Write};

use crate::{MappingResult, PairedMappingResult, InsertSizeStats, QualityMappingResult};

/// SAM FLAG bits
pub mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MAT_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MAT_REVERSE: u16 = 0x20;
    pub const FIRST: u16 = 0x40;
    pub const LAST: u16 = 0x80;
    pub const SUPPLEMENTARY: u16 = 0x800;
    pub const SECONDARY: u16 = 0x100;
}

/// Writes SAM format output to a file.
pub struct SamWriter {
    file: File,
}

impl SamWriter {
    /// Create a new SamWriter that writes to the given file path.
    pub fn new(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self { file })
    }

    /// Write SAM header @SQ lines for each reference genome.
    /// Each entry is (genome_name, sequence_length).
    pub fn write_header(&mut self, genomes: &[(&str, usize)]) -> io::Result<()> {
        for (name, len) in genomes {
            writeln!(self.file, "@SQ\tSN:{}\tLN:{}", name, len)?;
        }
        Ok(())
    }

    /// Write a single SAM line for one mapping result.
    /// is_supplementary: set true for secondary mappings (not the best match).
    fn write_sam_line(
        &mut self,
        read_name: &str,
        read_seq: &str,
        result: &MappingResult,
        genome_name: &str,
        is_supplementary: bool,
    ) -> io::Result<()> {
        let mut sam_flag: u16 = 0;

        if result.cigar.is_empty() {
            sam_flag |= flag::UNMAPPED;
        }

        if is_supplementary {
            sam_flag |= flag::SUPPLEMENTARY;
        }

        if result.is_reverse {
            sam_flag |= flag::REVERSE;
        }

        // MAPQ: convert score (0.0-1.0) to Phred-scaled quality (0-60)
        let mapq = (result.score * 60.0) as u16;

        // POS is 1-based in SAM
        let pos = result.position + 1;

        writeln!(
            self.file,
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*",
            read_name,        // QNAME
            sam_flag,         // FLAG
            genome_name,      // RNAME
            pos,              // POS (1-based)
            mapq,             // MAPQ
            result.cigar,     // CIGAR
            read_seq,         // SEQ
        )
    }

    /// Write SAM lines for all mapping results of a single read.
    /// The first (best) result gets no supplementary flag.
    /// Additional results are marked as supplementary.
    pub fn write_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[MappingResult],
        genome_names: &[&str],
    ) -> io::Result<()> {
        if results.is_empty() {
            // Write unmapped read
            writeln!(
                self.file,
                "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t*",
                read_name,
                flag::UNMAPPED,
                read_seq,
            )?;
            return Ok(());
        }

        for (i, result) in results.iter().enumerate() {
            let is_supplementary = i > 0;
            let gname = if result.genome_id < genome_names.len() as u32 {
                genome_names[result.genome_id as usize]
            } else {
                "*"
            };

            self.write_sam_line(read_name, read_seq, result, gname, is_supplementary)?;
        }

        Ok(())
    }

    /// Write quality-aware mapping results with QUAL field.
    pub fn write_quality_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[QualityMappingResult],
        genome_names: &[&str],
    ) -> io::Result<()> {
        if results.is_empty() {
            writeln!(
                self.file,
                "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
                read_name,
                flag::UNMAPPED,
                read_seq,
                String::from_utf8_lossy(&results.first().map(|r| &r.quality_scores).cloned().unwrap_or_default()),
            )?;
            return Ok(());
        }

        for (i, result) in results.iter().enumerate() {
            let is_supplementary = i > 0;
            let gname = if result.genome_id < genome_names.len() as u32 {
                genome_names[result.genome_id as usize]
            } else {
                "*"
            };

            let mut sam_flag: u16 = 0;
            if is_supplementary {
                sam_flag |= flag::SUPPLEMENTARY;
            }

            if result.is_reverse {
                sam_flag |= flag::REVERSE;
            }

            let mapq = ((result.combined_score * 60.0) as u16).min(60);
            let pos = result.position + 1;

            // QUAL field: original quality scores from FASTQ
            let qual_str = String::from_utf8_lossy(&result.quality_scores);

            writeln!(
                self.file,
                "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tMQ:f:{}",
                read_name,
                sam_flag,
                gname,
                pos,
                mapq,
                result.cigar,
                read_seq,
                qual_str,
                result.quality_penalty,
            )?;
        }

        Ok(())
    }

    /// Write SAM lines for a paired-end read mapping.
    pub fn write_paired_mappings(
        &mut self,
        read_name: &str,
        pair_result: &PairedMappingResult,
        genome_names: &[&str],
        insert_stats: &InsertSizeStats,
    ) -> io::Result<()> {
        let is_paired = true;
        let mut flag1: u16 = flag::FIRST;
        let mut flag2: u16 = flag::LAST;

        if is_paired {
            flag1 |= flag::PAIRED;
            flag2 |= flag::PAIRED;
        }

        // Determine if each read is mapped
        let m1 = &pair_result.map1;
        let m2 = &pair_result.map2;

        match (m1, m2) {
            (Some(map1), Some(map2)) => {
                // Both reads mapped
                let proper_pair = insert_stats.is_proper_pair(pair_result.tlen);
                if proper_pair {
                    flag1 |= flag::PROPER_PAIR;
                    flag2 |= flag::PROPER_PAIR;
                }

                // Determine strand orientation
                if map1.is_reverse {
                    flag1 |= flag::REVERSE;
                }
                if map2.is_reverse {
                    flag2 |= flag::REVERSE;
                    flag1 |= flag::MAT_REVERSE;
                }

                // MRNM: mate reference name (same genome for proper pairs)
                let mate_name = if map2.genome_id < genome_names.len() as u32 {
                    genome_names[map2.genome_id as usize]
                } else {
                    "*"
                };
                let rname1 = if map1.genome_id < genome_names.len() as u32 {
                    genome_names[map1.genome_id as usize]
                } else {
                    "*"
                };

                let pos1 = (map1.position + 1) as i32;
                let pos2 = (map2.position + 1) as i32;
                let mapq1 = ((map1.score * 60.0) as u16).min(60);
                let mapq2 = ((map2.score * 60.0) as u16).min(60);
                let tlen = pair_result.tlen as i32;

                writeln!(
                    self.file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t*",
                    read_name, flag1, rname1, pos1, mapq1, map1.cigar, mate_name, pos2, tlen, map2.cigar
                )?;

                writeln!(
                    self.file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t*\t{}\t{}\t\t*",
                    read_name, flag2, mate_name, pos2, mapq2, map2.cigar, 0, 0
                )?;
            }
            (Some(map1), None) => {
                // R1 mapped, R2 unmapped
                flag1 |= flag::PROPER_PAIR;
                flag2 |= flag::MAT_UNMAPPED;

                let rname1 = if map1.genome_id < genome_names.len() as u32 {
                    genome_names[map1.genome_id as usize]
                } else {
                    "*"
                };

                let pos1 = (map1.position + 1) as i32;
                let mapq1 = ((map1.score * 60.0) as u16).min(60);

                writeln!(
                    self.file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*",
                    read_name, flag1, rname1, pos1, mapq1, map1.cigar
                )?;

                writeln!(
                    self.file,
                    "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t*\t*",
                    read_name, flag2
                )?;
            }
            (None, Some(map2)) => {
                // R1 unmapped, R2 mapped
                flag1 |= flag::UNMAPPED;
                flag1 |= flag::MAT_UNMAPPED;
                flag2 |= flag::PROPER_PAIR;

                let rname2 = if map2.genome_id < genome_names.len() as u32 {
                    genome_names[map2.genome_id as usize]
                } else {
                    "*"
                };

                let pos2 = (map2.position + 1) as i32;
                let mapq2 = ((map2.score * 60.0) as u16).min(60);

                writeln!(
                    self.file,
                    "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t*\t*",
                    read_name, flag1
                )?;

                writeln!(
                    self.file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*",
                    read_name, flag2, rname2, pos2, mapq2, map2.cigar
                )?;
            }
            (None, None) => {
                // Both unmapped
                flag1 |= flag::UNMAPPED;
                flag2 |= flag::UNMAPPED;

                writeln!(
                    self.file,
                    "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t*\t*",
                    read_name, flag1
                )?;

                writeln!(
                    self.file,
                    "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t*\t*",
                    read_name, flag2
                )?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn create_temp_path() -> String {
        let dir = std::env::temp_dir();
        let ns = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        dir.join(format!("bitpop_sam_test_{}_{}.sam", std::process::id(), ns))
            .to_str()
            .unwrap()
            .to_string()
    }

    #[test]
    fn test_sam_header() {
        let path = create_temp_path();
        let mut writer = SamWriter::new(&path).unwrap();

        let genomes = vec![("chr1", 1000usize), ("chr2", 2000usize)];
        writer.write_header(&genomes).unwrap();
        drop(writer);

        let mut content = String::new();
        std::fs::File::open(&path).unwrap().read_to_string(&mut content).unwrap();
        let _ = std::fs::remove_file(&path);

        let lines: Vec<&str> = content.trim().lines().collect();
        assert_eq!(lines[0], "@SQ\tSN:chr1\tLN:1000");
        assert_eq!(lines[1], "@SQ\tSN:chr2\tLN:2000");
    }

    #[test]
    fn test_sam_single_mapping() {
        let path = create_temp_path();
        let mut writer = SamWriter::new(&path).unwrap();

        let result = MappingResult {
            genome_id: 0,
            position: 100,
            score: 0.95,
            cigar: "50M".to_string(),
            context: String::new(),
            is_reverse: false,
        };

        writer
            .write_mappings("read1", "ACGTACGT", &[result], &["chr1"])
            .unwrap();
        drop(writer);

        let mut content = String::new();
        std::fs::File::open(&path).unwrap().read_to_string(&mut content).unwrap();
        let _ = std::fs::remove_file(&path);

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "read1");   // QNAME
        assert_eq!(fields[1], "0");       // FLAG (no flags)
        assert_eq!(fields[2], "chr1");    // RNAME
        assert_eq!(fields[3], "101");     // POS (1-based)
        assert_eq!(fields[5], "50M");     // CIGAR
        assert_eq!(fields[9], "ACGTACGT"); // SEQ
    }

    #[test]
    fn test_sam_multiple_mappings() {
        let path = create_temp_path();
        let mut writer = SamWriter::new(&path).unwrap();

        let results = vec![
            MappingResult {
                genome_id: 0,
                position: 100,
                score: 0.95,
                cigar: "50M".to_string(),
                context: String::new(),
                is_reverse: false,
            },
            MappingResult {
                genome_id: 1,
                position: 200,
                score: 0.80,
                cigar: "50M".to_string(),
                context: String::new(),
                is_reverse: false,
            },
        ];

        writer
            .write_mappings("read1", "ACGT", &results, &["chr1", "chr2"])
            .unwrap();
        drop(writer);

        let mut content = String::new();
        std::fs::File::open(&path).unwrap().read_to_string(&mut content).unwrap();
        let _ = std::fs::remove_file(&path);

        let lines: Vec<&str> = content.trim().lines().collect();
        assert_eq!(lines.len(), 2);

        // First line: no supplementary flag
        let flags1: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(flags1[1], "0");

        // Second line: supplementary flag (0x800 = 2048)
        let flags2: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(flags2[1], "2048");
    }

    #[test]
    fn test_sam_unmapped_read() {
        let path = create_temp_path();
        let mut writer = SamWriter::new(&path).unwrap();

        writer
            .write_mappings("read1", "ACGT", &[], &["chr1"])
            .unwrap();
        drop(writer);

        let mut content = String::new();
        std::fs::File::open(&path).unwrap().read_to_string(&mut content).unwrap();
        let _ = std::fs::remove_file(&path);

        let line = content.trim();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "read1");
        assert_eq!(fields[1], "4");        // FLAG: UNMAPPED
        assert_eq!(fields[2], "*");        // RNAME: unmapped
        assert_eq!(fields[5], "*");        // CIGAR: unmapped
    }

    #[test]
    fn test_sam_mapq_calculation() {
        let path = create_temp_path();
        let mut writer = SamWriter::new(&path).unwrap();

        let result = MappingResult {
            genome_id: 0,
            position: 0,
            score: 1.0,
            cigar: "100M".to_string(),
            context: String::new(),
            is_reverse: false,
        };

        writer
            .write_mappings("read1", "ACGT", &[result], &["chr1"])
            .unwrap();
        drop(writer);

        let mut content = String::new();
        std::fs::File::open(&path).unwrap().read_to_string(&mut content).unwrap();
        let _ = std::fs::remove_file(&path);

        let fields: Vec<&str> = content.trim().split('\t').collect();
        assert_eq!(fields[4], "60"); // Perfect score → MAPQ 60
    }

    #[test]
    fn test_sam_empty_path() {
        let result = SamWriter::new("/nonexistent/dir/file.sam");
        assert!(result.is_err());
    }
}
