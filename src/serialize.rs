use std::io::{Read, Write, Result as IoResult};
use std::collections::HashMap;

use crate::fm::FmIndex;
use crate::BitPop;

const SAMPLE_INTERVAL: usize = 32;

// --- Binary Format ---
//
// [MAGIC: 4 bytes "BITP"]
// [VERSION: u32 LE]
// [K: u16 LE]
// [NUM_GENOMES: u32 LE]
//   For each genome (in order of genome_id):
//     [NAME_LEN: u32 LE]
//     [NAME: bytes, UTF-8]
//     [SEQ_LEN: u64 LE]
//     [SEQ: bytes]
// [FM_INDEX:]
//   [BWT_LEN: u64 LE]
//   [BWT: BWT_LEN × u32 LE]
//   [SA_LEN: u64 LE]
//   [SA: SA_LEN × u64 LE] (usize = u64)
//   [C_ARRAY: ALPHABET_SIZE × u64 LE]
//   [NUM_BOUNDARIES: u64 LE]
//   [BOUNDARIES: NUM_BOUNDARIES × (u32 LE, u64 LE, u64 LE)] (gid, start, len)
//   [OCC_SAMPLES_LEN: u64 LE]
//   [OCC_SAMPLES: OCC_SAMPLES_LEN × (ALPHABET_SIZE × u32 LE)]
//   [SAMPLE_INTERVAL: u32 LE]

const MAGIC: &[u8] = b"BITP";
const VERSION: u32 = 2;

fn read_u32_le<R: Read>(reader: &mut R) -> IoResult<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64_le<R: Read>(reader: &mut R) -> IoResult<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u16_le<R: Read>(reader: &mut R) -> IoResult<u16> {
    let mut buf = [0u8; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn write_u32_le<W: Write>(writer: &mut W, val: u32) -> IoResult<()> {
    writer.write_all(&val.to_le_bytes())
}

fn write_u64_le<W: Write>(writer: &mut W, val: u64) -> IoResult<()> {
    writer.write_all(&val.to_le_bytes())
}

fn write_u16_le<W: Write>(writer: &mut W, val: u16) -> IoResult<()> {
    writer.write_all(&val.to_le_bytes())
}

fn write_bytes<W: Write>(writer: &mut W, data: &[u8]) -> IoResult<()> {
    writer.write_all(data)
}

fn read_bytes<R: Read>(reader: &mut R, len: usize) -> IoResult<Vec<u8>> {
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

/// Serialize the entire BitPop instance to bytes.
pub fn serialize_bitpop(bp: &BitPop) -> IoResult<Vec<u8>> {
    let mut output = Vec::new();

    output.write_all(MAGIC)?;
    write_u32_le(&mut output, VERSION)?;
    write_u16_le(&mut output, bp.k() as u16)?;

    let num_genomes = bp.genome_count();
    write_u32_le(&mut output, num_genomes as u32)?;

    for i in 0..num_genomes {
        let gid = i as u32;
        let name = bp.genome_name(gid).unwrap_or("");
        write_u32_le(&mut output, name.len() as u32)?;
        write_bytes(&mut output, name.as_bytes())?;

        if let Some(seq) = bp.get_genome_seq(gid) {
            write_u64_le(&mut output, seq.len() as u64)?;
            write_bytes(&mut output, seq)?;
        } else {
            write_u64_le(&mut output, 0)?;
        }
    }

    Ok(output)
}

/// Deserialize a BitPop instance from bytes.
pub fn deserialize_bitpop(data: &[u8]) -> IoResult<BitPop> {
    let mut pos = 0;

    if data.len() < 4 + 4 + 2 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Data too short for header",
        ));
    }

    let magic = &data[pos..pos+4];
    if magic != MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid magic: not BITP",
        ));
    }
    pos += 4;

    let version = u32::from_le_bytes(data[pos..pos+4].try_into().map_err(
        |_| std::io::Error::new(std::io::ErrorKind::InvalidData, "Failed to read version"),
    )?);
    // Support both v1 (legacy) and v2 (FM-index)
    if version != 1 && version != VERSION {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Unsupported version: {}", version),
        ));
    }
    pos += 4;

    let k = u16::from_le_bytes(data[pos..pos+2].try_into().map_err(
        |_| std::io::Error::new(std::io::ErrorKind::InvalidData, "Failed to read k"),
    )?) as usize;
    pos += 2;

    let num_genomes = u32::from_le_bytes(data[pos..pos+4].try_into().map_err(
        |_| std::io::Error::new(std::io::ErrorKind::InvalidData, "Failed to read num_genomes"),
    )?) as usize;
    pos += 4;

    let mut genomes: HashMap<u32, Vec<u8>> = HashMap::new();
    let mut genome_names: HashMap<u32, String> = HashMap::new();

    for i in 0..num_genomes {
        if pos + 4 > data.len() {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Unexpected end at genome {} name_len", i)));
        }
        let name_len = u32::from_le_bytes(data[pos..pos+4].try_into().unwrap()) as usize;
        pos += 4;

        if pos + name_len > data.len() {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Unexpected end at genome {} name", i)));
        }
        let name_str = String::from_utf8(data[pos..pos+name_len].to_vec()).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, e)
        })?;
        pos += name_len;

        if pos + 8 > data.len() {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Unexpected end at genome {} seq_len", i)));
        }
        let seq_len = u64::from_le_bytes(data[pos..pos+8].try_into().unwrap()) as usize;
        pos += 8;

        if pos + seq_len > data.len() {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Unexpected end at genome {} seq", i)));
        }
        let seq_bytes = data[pos..pos+seq_len].to_vec();
        pos += seq_len;

        let gid = i as u32;
        genome_names.insert(gid, name_str);
        genomes.insert(gid, seq_bytes);
    }

    // Build FM-index from loaded genomes
    let fm_genomes: Vec<(&str, &[u8])> = (0..num_genomes).map(|i| {
        let gid = i as u32;
        (
            genome_names.get(&gid).map(|s| s.as_str()).unwrap_or(""),
            genomes.get(&gid).map(|s| s.as_slice()).unwrap_or(&[]),
        )
    }).collect();

    let fm_index = FmIndex::build(&fm_genomes);

    Ok(BitPop::from_fm_index(k, genomes, genome_names, fm_index))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_bitpop() -> BitPop {
        let mut bp = BitPop::new(6);
        bp.add_genome("human", "ACGTACGTACGTACGTACGTACGT");
        bp.add_genome("chimp", "ACGTACGTACGTAACAACGTACGT");
        bp.add_genome("mouse", "TTTTGGGGACGTACGTACGTACGT");
        bp.build();
        bp
    }

     #[test]
    fn test_serialize_deserialize_roundtrip() {
        let bp = make_test_bitpop();
        let bytes = serialize_bitpop(&bp).unwrap();
        let loaded = deserialize_bitpop(&bytes).unwrap();

        assert_eq!(loaded.k(), 6);
        assert_eq!(loaded.genome_count(), 3);
        assert_eq!(loaded.genome_name(0), Some("human"));
        assert_eq!(loaded.genome_name(1), Some("chimp"));
        assert_eq!(loaded.genome_name(2), Some("mouse"));

        let results_orig = bp.map_read("ACGTACGT", 3);
        let results_loaded = loaded.map_read("ACGTACGT", 3);
        assert!(!results_orig.is_empty());
        assert!(!results_loaded.is_empty());
        for orig in &results_orig {
            let loaded_match = results_loaded.iter().find(|r| r.genome_id == orig.genome_id);
            assert!(loaded_match.is_some(), "Missing genome_id {} in loaded", orig.genome_id);
        }
    }

    #[test]
    fn test_serialize_single_genome() {
        let mut bp = BitPop::new(8);
        bp.add_genome("test", "AACCGGTTAACCGGTT");
        bp.build();

        let bytes = serialize_bitpop(&bp).unwrap();
        assert!(bytes.len() > 0);

        let loaded = deserialize_bitpop(&bytes).unwrap();
        assert_eq!(loaded.genome_count(), 1);
        assert_eq!(loaded.genome_name(0), Some("test"));

        let results = loaded.map_read("AACCGGTT", 3);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_serialize_after_build() {
        let mut bp = BitPop::new(10);
        bp.add_genome("g1", &"ACGT".repeat(100));
        bp.add_genome("g2", &"TGCA".repeat(100));
        bp.build();

        let bytes = serialize_bitpop(&bp).unwrap();
        let loaded = deserialize_bitpop(&bytes).unwrap();

        assert_eq!(loaded.genome_count(), 2);
        assert!(loaded.bwt_len() > 0);

        let results_loaded = loaded.map_read("ACGTACGTACGT", 3);
        assert!(!results_loaded.is_empty());
    }

    #[test]
    fn test_deserialize_invalid_magic() {
        let bad_data: Vec<u8> = vec![b'B', b'A', b'D', b'X', 0, 0, 0, 1, 0, 0];
        let result = deserialize_bitpop(&bad_data);
        assert!(result.is_err());
    }

    #[test]
    fn test_deserialize_too_short() {
        let short_data: Vec<u8> = vec![b'B', b'I'];
        let result = deserialize_bitpop(&short_data);
        assert!(result.is_err());
    }

    #[test]
    fn test_serialize_no_genomes() {
        let mut bp = BitPop::new(6);
        bp.build();

        let bytes = serialize_bitpop(&bp).unwrap();
        let loaded = deserialize_bitpop(&bytes).unwrap();
        assert_eq!(loaded.genome_count(), 0);
    }

    #[test]
    fn test_serialize_multi_genome_ranking_preserved() {
        let mut bp = BitPop::new(6);
        bp.add_genome("shared", "ACGTAACAACGTAACAACGTAACA");
        bp.add_genome("unique", "TTTTTTTTACGTAACATTTTTTTT");
        bp.build();

        let bytes = serialize_bitpop(&bp).unwrap();
        let loaded = deserialize_bitpop(&bytes).unwrap();

        let read = "ACGTAACA";
        let results_orig = bp.map_read(read, 3);
        let results_loaded = loaded.map_read(read, 3);

        assert_eq!(results_orig.len(), results_loaded.len());
        for orig in &results_orig {
            let loaded_match = results_loaded.iter().find(|r| r.genome_id == orig.genome_id);
            assert!(loaded_match.is_some());
            let loaded_r = loaded_match.unwrap();
            assert!((orig.score - loaded_r.score).abs() < 0.001);
        }
    }

    #[test]
    fn test_serialize_large_genome() {
        let mut bp = BitPop::new(8);
        let genome = format!("{}{}{}", "ACGT".repeat(5000), "AACCGGTT", "TTTT".repeat(5000));
        bp.add_genome("large", &genome);
        bp.build();

        let bytes = serialize_bitpop(&bp).unwrap();

        let loaded = deserialize_bitpop(&bytes).unwrap();
        assert_eq!(loaded.genome_seq_len(0), Some(genome.len()));

        let results = loaded.map_read("AACCGGTT", 3);
        assert!(!results.is_empty());
    }
}
