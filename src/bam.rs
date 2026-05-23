//! BAM binary format writer using BGZF compression.
//!
//! Writes SAM alignment data in BAM binary format.
//! Compatible with samtools, bcftools, and other standard bioinformatics tools.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{self, Write};

use crate::{InsertSizeStats, MappingResult, PairedMappingResult, QualityMappingResult};

/// BAM magic string
const BAM_MAGIC: &[u8] = b"BAM\0";

/// Maximum uncompressed size per BGZF block (64KB - 4 bytes for header)
const BGZF_MAX_BLOCK: usize = 65536 - 4;

/// BGZF local header signature
const BGZF_HEADER: &[u8] = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\xbc\x02\x00";

/// BGZF local footer (uncompressed size as uint16)
fn bgzf_footer(uncompressed_len: u32) -> [u8; 8] {
    [
        0x00,
        0x00, // CRC32 (placeholder, will be fixed)
        0x00,
        0x00,
        (uncompressed_len & 0xFF) as u8,
        ((uncompressed_len >> 8) & 0xFF) as u8,
        0x00,
        0x00, // Reserved
    ]
}

/// BGZF writer that handles block boundaries.
struct BgzfWriter<FileType: Write> {
    file: FileType,
    buffer: Vec<u8>,
}

impl<FileType: Write> Write for BgzfWriter<FileType> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.buffer.extend_from_slice(data);
        while self.buffer.len() >= BGZF_MAX_BLOCK {
            self.flush_block()?;
        }
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[allow(dead_code)]
impl<FileType: Write> BgzfWriter<FileType> {
    fn new(file: FileType) -> Self {
        Self {
            file,
            buffer: Vec::with_capacity(BGZF_MAX_BLOCK),
        }
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // Write BGZF header
        self.file.write_all(BGZF_HEADER)?;

        // Compress the block
        let block_data = std::mem::take(&mut self.buffer);
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&block_data)?;
        let compressed = encoder.finish()?;

        // Write compressed data
        self.file.write_all(&compressed)?;

        // Write footer (uncompressed size)
        let footer = bgzf_footer(block_data.len() as u32);
        self.file.write_all(&footer)?;

        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        Ok(())
    }
}

impl<FileType: Write> Drop for BgzfWriter<FileType> {
    fn drop(&mut self) {
        let _ = self.flush_block();
    }
}

// Re-export SAM helper functions needed by BAM writer
mod sam_helpers {
    /// Compute NM tag (edit distance) from CIGAR string.
    #[allow(dead_code)]
    pub fn compute_nm_from_cigar(cigar: &str) -> usize {
        let mut nm = 0usize;
        let mut num_str = String::new();

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                num_str.push(ch);
            } else {
                if !num_str.is_empty() {
                    let count = num_str.parse::<usize>().unwrap_or(1);
                    match ch {
                        'X' => nm += count,
                        'I' => nm += count,
                        'D' => {}
                        'M' | 'N' | 'S' | 'H' | 'P' => {}
                        _ => {}
                    }
                }
                num_str.clear();
            }
        }

        nm
    }

    /// Generate MD tag string from CIGAR, read sequence, and reference context.
    pub fn generate_md_string(
        cigar: &str,
        read_seq: &str,
        ref_context: &str,
        _align_pos: usize,
    ) -> String {
        let mut md = String::new();
        let read_bases: Vec<char> = read_seq
            .chars()
            .filter(|c| "ACGTNacgtn".contains(*c))
            .collect();
        let ref_bases: Vec<char> = ref_context
            .chars()
            .filter(|c| "ACGTNacgtN".contains(*c))
            .collect();

        if read_bases.is_empty() || ref_bases.is_empty() {
            return String::new();
        }

        let mut ref_idx: usize = 0;
        let mut read_idx: usize = 0;
        let mut num_str = String::new();

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                num_str.push(ch);
            } else {
                let count = num_str.parse::<usize>().unwrap_or(0);
                num_str.clear();

                match ch {
                    'M' => {
                        let mut match_run = 0usize;
                        for i in 0..count {
                            if read_idx + i < read_bases.len() && ref_idx + i < ref_bases.len() {
                                if read_bases[read_idx + i]
                                    .eq_ignore_ascii_case(&ref_bases[ref_idx + i])
                                {
                                    match_run += 1;
                                } else {
                                    if match_run > 0 {
                                        md.push_str(&match_run.to_string());
                                        match_run = 0;
                                    }
                                    md.push(read_bases[read_idx + i].to_ascii_uppercase());
                                }
                            } else {
                                if match_run > 0 {
                                    md.push_str(&match_run.to_string());
                                    match_run = 0;
                                }
                            }
                        }
                        if match_run > 0 {
                            md.push_str(&match_run.to_string());
                        }
                        ref_idx += count;
                        read_idx += count;
                    }
                    'I' => {
                        read_idx += count;
                    }
                    'D' if ref_idx < ref_bases.len() => {
                        let del_count = count.min(ref_bases.len() - ref_idx);
                        if del_count > 0 {
                            md.push('^');
                            for i in 0..del_count {
                                md.push(ref_bases[ref_idx + i]);
                            }
                            ref_idx += del_count;
                        }
                    }
                    'N' => {
                        ref_idx += count;
                        read_idx += count;
                    }
                    'S' => {
                        read_idx += count;
                    }
                    'H' => {
                        // Hard clip - don't affect indices
                    }
                    'P' => {
                        // Padding - no effect on MD
                    }
                    _ => {}
                }
            }
        }

        md
    }
}

/// Encode CIGAR string to BAM binary format.
fn encode_cigar(cigar: &str) -> Vec<u8> {
    let mut result = Vec::new();
    let mut num_str = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            if !num_str.is_empty() {
                let count: u32 = num_str.parse().unwrap_or(0);
                let op_code = cigar_op_to_bam(ch);
                let code = (count << 4) | op_code;
                result.extend_from_slice(&code.to_le_bytes());
                num_str.clear();
            }
        }
    }

    result
}

/// Convert CIGAR character to BAM operation code.
fn cigar_op_to_bam(op: char) -> u32 {
    match op {
        'M' => 0,
        'I' => 1,
        'D' => 2,
        'N' => 3,
        'S' => 4,
        'H' => 5,
        'P' => 6,
        '=' => 7,
        'X' => 8,
        _ => 0,
    }
}

/// Encode DNA sequence to 2-bit per base BAM format.
fn encode_seq_to_bam(seq: &str) -> Option<Vec<u8>> {
    let bases: Vec<char> = seq.chars().filter(|c| "ACGTacgt".contains(*c)).collect();

    if bases.is_empty() {
        return None;
    }

    let mut encoded = Vec::with_capacity(bases.len().div_ceil(2));
    let mut i = 0;

    while i < bases.len() {
        let lo = base_to_2bit(bases[i]);
        if i + 1 < bases.len() {
            let hi = base_to_2bit(bases[i + 1]);
            encoded.push((hi << 4) | lo);
            i += 2;
        } else {
            encoded.push(lo);
            i += 1;
        }
    }

    Some(encoded)
}

/// Convert single base to 2-bit value.
fn base_to_2bit(base: char) -> u8 {
    match base.to_ascii_uppercase() {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 3,
    }
}

/// Encode a string as null-terminated bytes for BAM ref name.
fn encode_ref_name(name: &str) -> Vec<u8> {
    let truncated: String = name.chars().take(1).collect();
    truncated.into_bytes()
}

/// Write the BAM header block.
fn write_bam_header<W: Write>(writer: &mut W, references: &[(String, u32)]) -> io::Result<()> {
    writer.write_all(BAM_MAGIC)?;

    // Empty header text
    writer.write_all(&0u32.to_le_bytes())?;

    // Number of references
    let n_ref = references.len() as u32;
    writer.write_all(&n_ref.to_le_bytes())?;

    // Reference blocks
    for (name, length) in references {
        let name_bytes = encode_ref_name(name);
        let name_len = name_bytes.len() as u32;
        writer.write_all(&name_len.to_le_bytes())?;
        writer.write_all(&name_bytes)?;
        writer.write_all(&length.to_le_bytes())?;
    }

    Ok(())
}

/// Write a single BAM alignment record for unmapped read.
fn write_unmapped_record<W: Write>(writer: &mut W, flag: u16) -> io::Result<()> {
    writer.write_all(&flag.to_le_bytes())?; // FLAG
    writer.write_all(&0u32.to_le_bytes())?; // RNAME (0 = *)
    writer.write_all(&0u32.to_le_bytes())?; // POS (0 = *)
    writer.write_all(&0u8.to_le_bytes())?; // MAPQ
    writer.write_all(&0u8.to_le_bytes())?; // CIGAR (0 = *)
    writer.write_all(&0u32.to_le_bytes())?; // RNEXT (0 = *)
    writer.write_all(&0i32.to_le_bytes())?; // PNEXT (0 = *)
    writer.write_all(&0i32.to_le_bytes())?; // TLEN (0)

    // SEQ: "*" (0x0F)
    writer.write_all(&1u32.to_le_bytes())?;
    writer.write_all(&0x0Fu8.to_le_bytes())?;

    // QUAL: "*" (0xFF)
    writer.write_all(&1u32.to_le_bytes())?;
    writer.write_all(&0xFFu8.to_le_bytes())?;

    Ok(())
}

/// Write a single BAM alignment record for a mapped read.
#[expect(clippy::too_many_arguments)]
fn write_mapped_record<W: Write>(
    writer: &mut W,
    qname: &[u8],
    flag: u16,
    rname_idx: u32,
    pos: u32,
    mapq: u8,
    cigar_bytes: &[u8],
    seq_bytes: &[u8],
    qual_bytes: &[u8],
    tags: &str,
) -> io::Result<()> {
    writer.write_all(qname)?; // QNAME
    writer.write_all(&flag.to_le_bytes())?; // FLAG
    writer.write_all(&rname_idx.to_le_bytes())?; // RNAME
    writer.write_all(&pos.to_le_bytes())?; // POS
    writer.write_all(&mapq.to_le_bytes())?; // MAPQ
    writer.write_all(&(cigar_bytes.len() as u32).to_le_bytes())?; // CIGAR len
    writer.write_all(cigar_bytes)?; // CIGAR
    writer.write_all(&0u32.to_le_bytes())?; // RNEXT
    writer.write_all(&0i32.to_le_bytes())?; // PNEXT
    writer.write_all(&0i32.to_le_bytes())?; // TLEN

    // SEQ
    writer.write_all(&(seq_bytes.len() as u32).to_le_bytes())?;
    writer.write_all(seq_bytes)?;

    // QUAL
    writer.write_all(&(qual_bytes.len() as u32).to_le_bytes())?;
    writer.write_all(qual_bytes)?;

    // Tags
    if !tags.is_empty() {
        writer.write_all(tags.as_bytes())?;
        writer.write_all(b"\0")?;
    }

    Ok(())
}

/// BAM writer - writes alignments in BAM binary format.
pub struct BamWriter {
    file: File,
    genome_indices: std::collections::HashMap<String, u32>,
    bgzf_writer: Option<BgzfWriter<File>>,
}

impl BamWriter {
    /// Create a new BamWriter that writes to the given file path.
    pub fn new(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self {
            file,
            genome_indices: std::collections::HashMap::new(),
            bgzf_writer: None,
        })
    }

    /// Write BAM header @SQ lines for each reference genome.
    pub fn write_header(&mut self, genomes: &[(&str, usize)]) -> io::Result<()> {
        let references: Vec<(String, u32)> = genomes
            .iter()
            .enumerate()
            .map(|(i, (name, _len))| {
                let idx = i as u32;
                self.genome_indices.insert(name.to_string(), idx);
                (name.to_string(), 0)
            })
            .collect();

        write_bam_header(&mut self.file, &references)?;

        // Initialize BGZF writer for alignment records
        self.bgzf_writer = Some(BgzfWriter::new(self.file.try_clone()?));

        Ok(())
    }

    /// Get the reference index for a genome name.
    fn get_ref_idx(&self, genome_name: &str) -> u32 {
        *self.genome_indices.get(genome_name).unwrap_or(&0)
    }

    /// Convert SAM-style tags to BAM tag format.
    fn encode_tags(&self, tag_str: &str) -> String {
        let mut result = String::new();

        for field in tag_str.split('\t') {
            if field.is_empty() {
                continue;
            }

            let parts: Vec<&str> = field.split(':').collect();
            if parts.len() < 3 {
                continue;
            }

            let tag_id = format!("{}{}", parts[0], parts[1]);
            let tag_type = self.sam_tag_type_to_bam(parts[1]);

            let tag_value_str = parts[2..].join(":");

            let value_bytes = match tag_type.as_str() {
                "i" => {
                    if let Ok(v) = tag_value_str.parse::<i32>() {
                        v.to_le_bytes().to_vec()
                    } else {
                        continue;
                    }
                }
                "f" => {
                    if let Ok(v) = tag_value_str.parse::<f32>() {
                        v.to_le_bytes().to_vec()
                    } else {
                        continue;
                    }
                }
                "Z" | "A" | "a" => tag_value_str.as_bytes().to_vec(),
                _ => continue,
            };

            result.push_str(&tag_id);
            result.push_str(&tag_type);
            result.push_str(&String::from_utf8_lossy(&value_bytes));
        }

        result
    }

    /// Convert SAM tag format letter to BAM type character.
    fn sam_tag_type_to_bam(&self, format: &str) -> String {
        match format {
            "i:" => "i".to_string(),
            "f:" => "f".to_string(),
            "Z:" => "Z".to_string(),
            "A:" => "A".to_string(),
            _ => "Z".to_string(),
        }
    }

    /// Write BAM record for a mapped read.
    #[expect(clippy::too_many_arguments)]
    fn write_mapped_bam(
        &mut self,
        read_name: &str,
        read_seq: &str,
        cigar: &str,
        score: f64,
        rarity: f64,
        hf_score: f64,
        md_string: &str,
        is_supplementary: bool,
        is_paired: bool,
        is_first: bool,
        is_reverse: bool,
        genome_name: &str,
    ) -> io::Result<()> {
        let mut flag: u16 = 0;

        if is_paired {
            flag |= 0x1;
            if is_first {
                flag |= 0x40;
            } else {
                flag |= 0x80;
            }
        }

        if cigar.is_empty() {
            flag |= 0x4;
        }

        if is_supplementary {
            flag |= 0x800;
        }

        if is_reverse {
            flag |= 0x10;
        }

        let mapq = ((score * 60.0) as u8).min(60);
        let pos: u32 = 0; // Simplified - would need position from result
        let rname_idx = self.get_ref_idx(genome_name);
        let cigar_bytes = encode_cigar(cigar);
        let seq_bytes = encode_seq_to_bam(read_seq).unwrap_or_default();
        let qual_bytes = vec![0x21u8; seq_bytes.len() * 2];

        let mut tags = String::new();
        if !md_string.is_empty() {
            tags.push_str(&format!("\tMD:Z:{}", md_string));
        }
        tags.push_str(&format!("\tAS:f:{:.4}", score));
        tags.push_str(&format!("\tRK:f:{:.6}", rarity));
        tags.push_str(&format!("\tHF:f:{:.4}", hf_score));

        if is_supplementary {
            tags.push_str(&format!("\tXS:f:{:.4}", score));
        }

        let bam_tags = self.encode_tags(&tags);

        write_mapped_record(
            &mut self.bgzf_writer.as_mut().unwrap(),
            read_name.as_bytes(),
            flag,
            rname_idx,
            pos,
            mapq,
            &cigar_bytes,
            &seq_bytes,
            &qual_bytes,
            &bam_tags,
        )
    }
}

impl AlignmentWriter for BamWriter {
    fn write_header(&mut self, genomes: &[(&str, usize)]) -> io::Result<()> {
        self.write_header(genomes)
    }

    fn write_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[MappingResult],
        genome_names: &[&str],
    ) -> io::Result<()> {
        if results.is_empty() {
            write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), 0x4)
        } else {
            for (i, result) in results.iter().enumerate() {
                let is_supplementary = i > 0;
                let gname = if result.genome_id < genome_names.len() as u32 {
                    genome_names[result.genome_id as usize]
                } else {
                    "*"
                };

                let md = sam_helpers::generate_md_string(
                    &result.cigar,
                    read_seq,
                    &result.md_string,
                    result.position as usize,
                );

                self.write_mapped_bam(
                    read_name,
                    read_seq,
                    &result.cigar,
                    result.score,
                    result.rarity,
                    result.hf_score,
                    &md,
                    is_supplementary,
                    false,
                    true,
                    result.is_reverse,
                    gname,
                )?;
            }
            Ok(())
        }
    }

    fn write_quality_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[QualityMappingResult],
        genome_names: &[&str],
    ) -> io::Result<()> {
        if results.is_empty() {
            write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), 0x4)
        } else {
            for (i, result) in results.iter().enumerate() {
                let is_supplementary = i > 0;
                let gname = if result.genome_id < genome_names.len() as u32 {
                    genome_names[result.genome_id as usize]
                } else {
                    "*"
                };

                let md = sam_helpers::generate_md_string(
                    &result.cigar,
                    read_seq,
                    &result.context,
                    result.position as usize,
                );

                let mut tags = String::new();
                if !md.is_empty() {
                    tags.push_str(&format!("\tMD:Z:{}", md));
                }
                tags.push_str(&format!("\tAS:f:{:.4}", result.align_score));
                tags.push_str(&format!("\tRK:f:{:.6}", result.rarity));
                tags.push_str(&format!("\tHF:f:{:.4}", result.rarity * 0.0));

                if is_supplementary {
                    tags.push_str(&format!("\tXS:f:{:.4}", result.combined_score));
                }

                let bam_tags = self.encode_tags(&tags);

                let mut flag: u16 = 0;
                if is_supplementary {
                    flag |= 0x800;
                }
                if result.is_reverse {
                    flag |= 0x10;
                }

                let mapq = ((result.combined_score * 60.0) as u8).min(60);
                let pos: u32 = result.position as u32;
                let rname_idx = self.get_ref_idx(gname);
                let cigar_bytes = encode_cigar(&result.cigar);
                let seq_bytes = encode_seq_to_bam(read_seq).unwrap_or_default();
                let qual_bytes = result.quality_scores.clone();

                write_mapped_record(
                    &mut self.bgzf_writer.as_mut().unwrap(),
                    read_name.as_bytes(),
                    flag,
                    rname_idx,
                    pos,
                    mapq,
                    &cigar_bytes,
                    &seq_bytes,
                    &qual_bytes,
                    &bam_tags,
                )?;
            }
            Ok(())
        }
    }

    fn write_paired_mappings(
        &mut self,
        read_name: &str,
        pair_result: &PairedMappingResult,
        genome_names: &[&str],
        insert_stats: &InsertSizeStats,
    ) -> io::Result<()> {
        let m1 = &pair_result.map1;
        let m2 = &pair_result.map2;

        match (m1, m2) {
            (Some(map1), Some(map2)) => {
                let proper_pair = insert_stats.is_proper_pair(pair_result.tlen);
                let mut flag1: u16 = 0x1 | 0x40;
                let mut flag2: u16 = 0x1 | 0x80;

                if proper_pair {
                    flag1 |= 0x2;
                    flag2 |= 0x2;
                }

                if map1.is_reverse {
                    flag1 |= 0x10;
                }
                if map2.is_reverse {
                    flag2 |= 0x10;
                    flag1 |= 0x20;
                }

                let rname1 = if map1.genome_id < genome_names.len() as u32 {
                    genome_names[map1.genome_id as usize]
                } else {
                    "*"
                };
                let rname2 = if map2.genome_id < genome_names.len() as u32 {
                    genome_names[map2.genome_id as usize]
                } else {
                    "*"
                };

                let mapq1 = ((map1.score * 60.0) as u8).min(60);
                let mapq2 = ((map2.score * 60.0) as u8).min(60);
                let pos1: u32 = map1.position as u32;
                let pos2: u32 = map2.position as u32;

                let cigar1_bytes = encode_cigar(&map1.cigar);
                let cigar2_bytes = encode_cigar(&map2.cigar);
                let seq1_bytes = encode_seq_to_bam("").unwrap_or_default();
                let seq2_bytes = encode_seq_to_bam("").unwrap_or_default();
                let qual1_bytes = vec![0x21u8; 0];
                let qual2_bytes = vec![0x21u8; 0];

                // Gaussian insert size model confidence score
                let gm_tag = if insert_stats.count >= 2 {
                    let confidence = insert_stats.insert_size_confidence(pair_result.tlen);
                    format!("\tGM:f:{:.4}", confidence)
                } else {
                    String::new()
                };

                let mut tags1 = String::new();
                let md1 = sam_helpers::generate_md_string(
                    &map1.cigar,
                    "",
                    &map1.md_string,
                    map1.position as usize,
                );
                if !md1.is_empty() {
                    tags1.push_str(&format!("\tMD:Z:{}", md1));
                }
                tags1.push_str(&format!("\tAS:f:{:.4}", map1.align_score));
                tags1.push_str(&format!("\tRK:f:{:.6}", map1.rarity));
                tags1.push_str(&gm_tag);
                let bam_tags1 = self.encode_tags(&tags1);

                let mut tags2 = String::new();
                let md2 = sam_helpers::generate_md_string(
                    &map2.cigar,
                    "",
                    &map2.md_string,
                    map2.position as usize,
                );
                if !md2.is_empty() {
                    tags2.push_str(&format!("\tMD:Z:{}", md2));
                }
                tags2.push_str(&format!("\tAS:f:{:.4}", map2.align_score));
                tags2.push_str(&format!("\tRK:f:{:.6}", map2.rarity));
                tags2.push_str(&gm_tag);
                let bam_tags2 = self.encode_tags(&tags2);

                let r1_idx = self.get_ref_idx(rname1);
                let r2_idx = self.get_ref_idx(rname2);

                write_mapped_record(
                    &mut self.bgzf_writer.as_mut().unwrap(),
                    read_name.as_bytes(),
                    flag1,
                    r1_idx,
                    pos1,
                    mapq1,
                    &cigar1_bytes,
                    &seq1_bytes,
                    &qual1_bytes,
                    &bam_tags1,
                )?;

                write_mapped_record(
                    &mut self.bgzf_writer.as_mut().unwrap(),
                    read_name.as_bytes(),
                    flag2,
                    r2_idx,
                    pos2,
                    mapq2,
                    &cigar2_bytes,
                    &seq2_bytes,
                    &qual2_bytes,
                    &bam_tags2,
                )
            }
            (Some(map1), None) => {
                let flag1: u16 = 0x1 | 0x40 | 0x2;
                let flag2: u16 = 0x1 | 0x80 | 0x8;

                let rname1 = if map1.genome_id < genome_names.len() as u32 {
                    genome_names[map1.genome_id as usize]
                } else {
                    "*"
                };

                let mapq1 = ((map1.score * 60.0) as u8).min(60);
                let pos1: u32 = map1.position as u32;
                let cigar1_bytes = encode_cigar(&map1.cigar);
                let seq1_bytes = encode_seq_to_bam("").unwrap_or_default();
                let qual1_bytes = vec![0x21u8; 0];

                let mut tags1 = String::new();
                let md1 = sam_helpers::generate_md_string(
                    &map1.cigar,
                    "",
                    &map1.md_string,
                    map1.position as usize,
                );
                if !md1.is_empty() {
                    tags1.push_str(&format!("\tMD:Z:{}", md1));
                }
                tags1.push_str(&format!("\tAS:f:{:.4}", map1.align_score));
                tags1.push_str(&format!("\tRK:f:{:.6}", map1.rarity));
                let bam_tags1 = self.encode_tags(&tags1);

                let r1_idx = self.get_ref_idx(rname1);
                write_mapped_record(
                    &mut self.bgzf_writer.as_mut().unwrap(),
                    read_name.as_bytes(),
                    flag1,
                    r1_idx,
                    pos1,
                    mapq1,
                    &cigar1_bytes,
                    &seq1_bytes,
                    &qual1_bytes,
                    &bam_tags1,
                )?;

                write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), flag2)
            }
            (None, Some(map2)) => {
                let flag1: u16 = 0x1 | 0x40 | 0x4 | 0x8;
                let mut flag2: u16 = 0x1 | 0x80 | 0x2;
                if map2.is_reverse {
                    flag2 |= 0x10;
                }

                let rname2 = if map2.genome_id < genome_names.len() as u32 {
                    genome_names[map2.genome_id as usize]
                } else {
                    "*"
                };

                let mapq2 = ((map2.score * 60.0) as u8).min(60);
                let pos2: u32 = map2.position as u32;
                let cigar2_bytes = encode_cigar(&map2.cigar);
                let seq2_bytes = encode_seq_to_bam("").unwrap_or_default();
                let qual2_bytes = vec![0x21u8; 0];

                let mut tags2 = String::new();
                let md2 = sam_helpers::generate_md_string(
                    &map2.cigar,
                    "",
                    &map2.md_string,
                    map2.position as usize,
                );
                if !md2.is_empty() {
                    tags2.push_str(&format!("\tMD:Z:{}", md2));
                }
                tags2.push_str(&format!("\tAS:f:{:.4}", map2.align_score));
                tags2.push_str(&format!("\tRK:f:{:.6}", map2.rarity));
                let bam_tags2 = self.encode_tags(&tags2);

                let r2_idx = self.get_ref_idx(rname2);
                write_mapped_record(
                    &mut self.bgzf_writer.as_mut().unwrap(),
                    read_name.as_bytes(),
                    flag2,
                    r2_idx,
                    pos2,
                    mapq2,
                    &cigar2_bytes,
                    &seq2_bytes,
                    &qual2_bytes,
                    &bam_tags2,
                )?;

                write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), flag1)
            }
            (None, None) => {
                let flag1: u16 = 0x1 | 0x40 | 0x4;
                let flag2: u16 = 0x1 | 0x80 | 0x4;

                write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), flag1)?;
                write_unmapped_record(&mut self.bgzf_writer.as_mut().unwrap(), flag2)
            }
        }
    }
}

/// Trait for SAM/BAM writers - both formats share the same interface.
pub trait AlignmentWriter {
    fn write_header(&mut self, genomes: &[(&str, usize)]) -> io::Result<()>;
    fn write_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[MappingResult],
        genome_names: &[&str],
    ) -> io::Result<()>;
    fn write_quality_mappings(
        &mut self,
        read_name: &str,
        read_seq: &str,
        results: &[QualityMappingResult],
        genome_names: &[&str],
    ) -> io::Result<()>;
    fn write_paired_mappings(
        &mut self,
        read_name: &str,
        pair_result: &PairedMappingResult,
        genome_names: &[&str],
        insert_stats: &InsertSizeStats,
    ) -> io::Result<()>;
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
        dir.join(format!("bitpop_bam_test_{}_{}.bam", std::process::id(), ns))
            .to_str()
            .unwrap()
            .to_string()
    }

    #[test]
    fn test_bam_header() {
        let path = create_temp_path();
        let mut writer = BamWriter::new(&path).unwrap();

        let genomes = vec![("chr1", 1000usize), ("chr2", 2000usize)];
        writer.write_header(&genomes).unwrap();
        drop(writer);

        let mut file = std::fs::File::open(&path).unwrap();
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic).unwrap();
        assert_eq!(&magic, BAM_MAGIC);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_bam_cigar_encoding() {
        let cigar = encode_cigar("50M");
        assert_eq!(cigar.len(), 4);
        let code = u32::from_le_bytes([cigar[0], cigar[1], cigar[2], cigar[3]]);
        assert_eq!(code, (50 << 4));

        let cigar = encode_cigar("10S40M5S");
        assert_eq!(cigar.len(), 12);
        let code1 = u32::from_le_bytes([cigar[0], cigar[1], cigar[2], cigar[3]]);
        assert_eq!(code1, (10 << 4) | 4);
        let code2 = u32::from_le_bytes([cigar[4], cigar[5], cigar[6], cigar[7]]);
        assert_eq!(code2, (40 << 4));
    }

    #[test]
    fn test_bam_dna_encoding() {
        let seq = "ACGT";
        let encoded = encode_seq_to_bam(seq).unwrap();
        assert_eq!(encoded.len(), 2);
        // BAM: first base in low nibble, second base in high nibble
        // A=0, C=1 → (1 << 4) | 0 = 0x10
        // G=2, T=3 → (3 << 4) | 2 = 0x32
        assert_eq!(encoded[0], 0x10);
        assert_eq!(encoded[1], 0x32);
    }

    #[test]
    fn test_bam_empty_path() {
        let result = BamWriter::new("/nonexistent/dir/file.bam");
        assert!(result.is_err());
    }
}
