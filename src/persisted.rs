use std::collections::HashMap;
use std::io::{Read, Result as IoResult};

use memmap2::Mmap;
use sha2::{Digest, Sha256};

use crate::fm::{FmIndex, OccCounter};
use crate::BitPop;

type GenomeData = (HashMap<u32, Vec<u8>>, HashMap<u32, String>);

// --- File Format Constants ---

const MAGIC: [u8; 4] = *b"BITP";
const VERSION: u32 = 8;
const HEADER_SIZE: usize = 64;
const SECTION_NAME_LEN: usize = 16;
const SECTION_HEADER_SIZE: usize = 48; // name[16] + offset(8) + comp_size(8) + decomp_size(8) + flags(8)

// Section names (padded to SECTION_NAME_LEN=16)
const SECTION_BWT_UNCOMP: [u8; 16] = *b"BWT_UNCOMP\0\0\0\0\0\0";
const SECTION_SA_UNCOMP: [u8; 16] = *b"SA_UNCOMP\0\0\0\0\0\0\0";
const SECTION_FM_INDEX: [u8; 16] = *b"FM_INDEX\0\0\0\0\0\0\0\0";
const SECTION_GENOMES: [u8; 16] = *b"GENOMES\0\0\0\0\0\0\0\0\0";
const SECTION_SPACED_SEED: [u8; 16] = *b"SPACED_SEED\0\0\0\0\0";
const SECTION_HF: [u8; 16] = *b"HF_PROFILES\0\0\0\0\0";

/// Number of sections in v5+ format (BWT_UNCOMP + SA_UNCOMP + FM_INDEX + GENOMES + SPACED_SEED + HF)
const NUM_SECTIONS_V5: usize = 6;

/// Represents a section in the persisted file.
#[allow(dead_code)]
struct SectionInfo {
    name: [u8; SECTION_NAME_LEN],
    offset: u64,
    compressed_size: u64,
    decompressed_size: u64,
    flags: u64,
}

/// Header at the start of every persisted file.
#[repr(C)]
struct FileHeader {
    magic: [u8; 4],      // "BITP"
    version: u32,        // format version
    k: u16,              // k-mer size
    num_genomes: u32,    // number of genomes
    _reserved: [u8; 46], // padding to 64 bytes
}

impl FileHeader {
    fn new(k: usize, num_genomes: usize) -> Self {
        Self {
            magic: MAGIC,
            version: VERSION,
            k: k as u16,
            num_genomes: num_genomes as u32,
            _reserved: [0u8; 46],
        }
    }
}

// --- Serialization (save) ---

fn serialize_spaced_seed(bp: &BitPop) -> Vec<u8> {
    let mut data = Vec::new();

    // Spaced seed pattern as binary string
    let pattern_str: String = bp
        .spaced_seed_pattern
        .iter()
        .map(|&b| if b { '1' } else { '0' })
        .collect();
    data.extend_from_slice(&(pattern_str.len() as u32).to_le_bytes());
    data.extend_from_slice(pattern_str.as_bytes());

    // Spaced seed hash table
    if let Some(hash_table) = &bp.spaced_seed_hash {
        data.extend_from_slice(&1u32.to_le_bytes()); // has_hash = true
        data.extend_from_slice(&(hash_table.len() as u64).to_le_bytes());

        for (key, positions) in hash_table {
            data.extend_from_slice(&key.to_le_bytes());
            data.extend_from_slice(&(positions.len() as u64).to_le_bytes());
            for &(gid, pos) in positions {
                data.extend_from_slice(&gid.to_le_bytes());
                data.extend_from_slice(&pos.to_le_bytes());
            }
        }
    } else {
        data.extend_from_slice(&0u32.to_le_bytes()); // has_hash = false
    }

    data
}

fn serialize_hf(bp: &BitPop) -> Vec<u8> {
    let mut data = Vec::new();
    data.extend_from_slice(&(bp.hf_profiles.len() as u64).to_le_bytes());
    for (gid, profile) in &bp.hf_profiles {
        data.extend_from_slice(&gid.to_le_bytes());
        let profile_bytes = crate::hf::serialize_hf_profile(profile);
        data.extend_from_slice(&(profile_bytes.len() as u64).to_le_bytes());
        data.extend_from_slice(&profile_bytes);
    }
    data
}

fn deserialize_spaced_seed(data: &[u8], bp: &mut BitPop) {
    let mut pos = 0;

    // Parse pattern
    if pos + 4 > data.len() {
        return;
    }
    let pattern_len = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
    pos += 4;

    if pos + pattern_len > data.len() {
        return;
    }
    let pattern_str = String::from_utf8(data[pos..pos + pattern_len].to_vec()).unwrap_or_default();
    pos += pattern_len;

    bp.spaced_seed_pattern = pattern_str.chars().map(|c| c == '1').collect();

    // Parse hash table
    if pos + 4 > data.len() {
        return;
    }
    let has_hash = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) != 0;
    pos += 4;

    if has_hash {
        if pos + 8 > data.len() {
            return;
        }
        let num_entries = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()) as usize;
        pos += 8;

        let mut hash_table = std::collections::HashMap::new();

        for _ in 0..num_entries {
            if pos + 8 + 8 + 8 > data.len() {
                break;
            }
            let key = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
            pos += 8;

            let num_positions = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()) as usize;
            pos += 8;

            let mut positions = Vec::with_capacity(num_positions);
            for _ in 0..num_positions {
                if pos + 16 > data.len() {
                    break;
                }
                let gid = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
                pos += 4;
                let _reserved = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
                pos += 8;
                positions.push((gid, _reserved));
            }

            hash_table.insert(key, positions);
        }

        bp.spaced_seed_hash = Some(hash_table);
    }
}

fn deserialize_hf_profiles(data: &[u8], bp: &mut BitPop) {
    let mut pos = 0;
    if pos + 8 > data.len() {
        return;
    }
    let num_profiles: usize = {
        let bytes: [u8; 8] = data[pos..pos + 8].try_into().unwrap_or_default();
        u64::from_le_bytes(bytes) as usize
    };
    pos += 8;

    for _ in 0..num_profiles {
        if pos + 8 + 8 > data.len() {
            break;
        }
        let gid: u32 = {
            let bytes: [u8; 4] = data[pos..pos + 4].try_into().unwrap_or_default();
            u32::from_le_bytes(bytes)
        };
        pos += 4;
        let _reserved: u32 = {
            let bytes: [u8; 4] = data[pos..pos + 4].try_into().unwrap_or_default();
            u32::from_le_bytes(bytes)
        };
        pos += 4;
        let profile_len: usize = {
            let bytes: [u8; 8] = data[pos..pos + 8].try_into().unwrap_or_default();
            u64::from_le_bytes(bytes) as usize
        };
        pos += 8;
        if pos + profile_len > data.len() {
            break;
        }
        let profile_bytes = &data[pos..pos + profile_len];
        pos += profile_len;
        if let Some(profile) = crate::hf::deserialize_hf_profile(profile_bytes) {
            bp.hf_profiles.insert(gid, profile);
        }
    }
}

/// Save a BitPop instance to a file using the memmap-friendly format.
/// Format v7: [header][section_table][BWT_COMPRESSED][SA_COMPRESSED][FM_INDEX][GENOMES][checksum]
/// Both BWT and SA use zstd compression for smaller file size.
pub fn save_bitpop(bp: &BitPop, path: &str) -> IoResult<()> {
    // 1. Serialize FM-Index (compressed fallback)
    let fm_data = serialize_fm_index(bp)?;
    let fm_compressed = zstd::encode_all(fm_data.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd FM compress failed: {}", e)))?;

    // 2. Serialize genomes (compressed)
    let genomes_data = serialize_genomes(bp)?;
    let genomes_compressed = zstd::encode_all(genomes_data.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd genomes compress failed: {}", e)))?;

    // 3. Serialize BWT with zstd compression
    let bwt_compressed = serialize_bwt_compressed(bp)?;
    let bwt_decompressed_len = bp.get_fm_index().map(|f| f.len()).unwrap_or(0);

    // 4. Serialize SA with delta + VLI compression
    let sa_compressed = serialize_sa_compressed(bp)?;

    // 5. Build section table (4 sections: BWT_UNCOMP, SA_UNCOMP, FM_INDEX, GENOMES)
    let mut section_table = Vec::new();

    let base_offset: u64 = (HEADER_SIZE + (NUM_SECTIONS_V5 * SECTION_HEADER_SIZE)) as u64;

    let mut offset = base_offset;
    write_section_header(
        &mut section_table,
        &SECTION_BWT_UNCOMP,
        offset,
        bwt_compressed.len() as u64,
        bwt_decompressed_len as u64,
        0,
    );
    offset += bwt_compressed.len() as u64;

    write_section_header(
        &mut section_table,
        &SECTION_SA_UNCOMP,
        offset,
        sa_compressed.len() as u64,
        sa_compressed.len() as u64,
        0,
    );
    offset += sa_compressed.len() as u64;

    write_section_header(
        &mut section_table,
        &SECTION_FM_INDEX,
        offset,
        fm_compressed.len() as u64,
        fm_data.len() as u64,
        0,
    );
    offset += fm_compressed.len() as u64;

    write_section_header(
        &mut section_table,
        &SECTION_GENOMES,
        offset,
        genomes_compressed.len() as u64,
        genomes_data.len() as u64,
        0,
    );
    offset += genomes_compressed.len() as u64;

    // 5b. Spaced seed section (optional - only if spaced_seed_hash exists)
    let spaced_seed_data = serialize_spaced_seed(bp);
    let spaced_seed_compressed = zstd::encode_all(spaced_seed_data.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd spaced seed compress failed: {}", e)))?;

    write_section_header(
        &mut section_table,
        &SECTION_SPACED_SEED,
        offset,
        spaced_seed_compressed.len() as u64,
        spaced_seed_data.len() as u64,
        0,
    );
    offset += spaced_seed_compressed.len() as u64;

    // 5c. HF profiles section (optional - only if hf_profiles exist)
    let hf_data = serialize_hf(bp);
    let hf_compressed = zstd::encode_all(hf_data.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd hf compress failed: {}", e)))?;

    write_section_header(
        &mut section_table,
        &SECTION_HF,
        offset,
        hf_compressed.len() as u64,
        hf_data.len() as u64,
        0,
    );
    offset += hf_compressed.len() as u64;

    // 6. Assemble file: header + section_table + sections
    let mut all_data = Vec::with_capacity((offset + hf_compressed.len() as u64 + 32) as usize);

    let header_placeholder = vec![0u8; HEADER_SIZE];
    all_data.extend_from_slice(&header_placeholder);
    all_data.extend_from_slice(&section_table);
    all_data.extend_from_slice(&bwt_compressed);
    all_data.extend_from_slice(&sa_compressed);
    all_data.extend_from_slice(&fm_compressed);
    all_data.extend_from_slice(&genomes_compressed);
    all_data.extend_from_slice(&spaced_seed_compressed);
    all_data.extend_from_slice(&hf_compressed);

    // 7. Fill in header
    let header = FileHeader::new(bp.k(), bp.genome_count());
    let mut header_bytes = Vec::with_capacity(HEADER_SIZE);
    header_bytes.extend_from_slice(&header.magic);
    header_bytes.extend_from_slice(&header.version.to_le_bytes());
    header_bytes.extend_from_slice(&header.k.to_le_bytes());
    header_bytes.extend_from_slice(&header.num_genomes.to_le_bytes());
    header_bytes.resize(HEADER_SIZE, 0u8);

    all_data[..HEADER_SIZE].copy_from_slice(&header_bytes);

    // 8. Compute checksum
    let mut hasher = Sha256::new();
    hasher.update(&all_data);
    let hash_bytes = hasher.finalize();
    let mut checksum = [0u8; 32];
    checksum.copy_from_slice(hash_bytes.as_ref());

    all_data.extend_from_slice(&checksum);

    std::fs::write(path, &all_data)?;
    Ok(())
}

fn write_section_header(
    buf: &mut Vec<u8>,
    name: &[u8],
    offset: u64,
    comp_size: u64,
    decomp_size: u64,
    _flags: u64,
) {
    let mut name_bytes = [0u8; SECTION_NAME_LEN];
    let len = name.len().min(SECTION_NAME_LEN);
    name_bytes[..len].copy_from_slice(&name[..len]);
    buf.extend_from_slice(&name_bytes);
    buf.extend_from_slice(&offset.to_le_bytes());
    buf.extend_from_slice(&comp_size.to_le_bytes());
    buf.extend_from_slice(&decomp_size.to_le_bytes());
    buf.extend_from_slice(&_flags.to_le_bytes());
}

fn serialize_fm_index(bp: &BitPop) -> IoResult<Vec<u8>> {
    let fm = match bp.get_fm_index() {
        Some(fm) => fm,
        None => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FM-Index not built",
            ))
        }
    };

    let mut data = Vec::new();
    let bwt_len = fm.len();

    // BWT length
    data.extend_from_slice(&(bwt_len as u64).to_le_bytes());

    // BWT: 2 bits per entry (4 values per byte)
    let bwt_packed_len = bwt_len.div_ceil(4);
    let mut bwt_packed = vec![0u8; bwt_packed_len];
    for i in 0..bwt_len {
        let v = fm.bwt_at(i) & 3;
        let byte_idx = i / 4;
        let bit_offset = 6 - (i % 4) * 2;
        bwt_packed[byte_idx] |= v << bit_offset;
    }
    data.extend_from_slice(&bwt_packed);

    // SA: u32 per entry
    let sa_len = fm.sa_len();
    data.extend_from_slice(&(sa_len as u64).to_le_bytes());
    for rank in 0..sa_len {
        data.extend_from_slice(&(fm.sa_at(rank) as u32).to_le_bytes());
    }

    // C-array: u32 x 5
    for j in 0..5 {
        data.extend_from_slice(&(fm.c_array(j) as u32).to_le_bytes());
    }

    // Genome boundaries
    let boundaries = fm.genome_boundaries();
    data.extend_from_slice(&(boundaries.len() as u64).to_le_bytes());
    for &(start, len, gid) in boundaries {
        data.extend_from_slice(&(start as u32).to_le_bytes());
        data.extend_from_slice(&(len as u32).to_le_bytes());
        data.extend_from_slice(&gid.to_le_bytes());
    }

    // Sample interval
    data.extend_from_slice(&32u32.to_le_bytes());

    // Sentinel mask
    let sentinel_mask_len = bwt_len.div_ceil(8);
    let mut sentinel_mask = vec![0u8; sentinel_mask_len];
    for i in 0..bwt_len {
        if fm.bwt_at(i) == 0 {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            sentinel_mask[byte_idx] |= 1u8 << bit_idx;
        }
    }
    data.extend_from_slice(&(sentinel_mask_len as u64).to_le_bytes());
    data.extend_from_slice(&sentinel_mask);

    Ok(data)
}

fn serialize_genomes(bp: &BitPop) -> IoResult<Vec<u8>> {
    let mut data = Vec::new();

    for i in 0..bp.genome_count() {
        let gid = i as u32;
        let name = bp.genome_name(gid).unwrap_or("");
        data.extend_from_slice(&(name.len() as u32).to_le_bytes());
        data.extend_from_slice(name.as_bytes());
        if let Some(seq) = bp.get_genome_seq(gid) {
            data.extend_from_slice(&(seq.len() as u64).to_le_bytes());
            data.extend_from_slice(seq);
        } else {
            data.extend_from_slice(&0u64.to_le_bytes());
        }
    }

    Ok(data)
}

/// Serialize BWT using zstd compression.
/// Format: [bwt_len: u64][comp_size: u64][zstd_data]
fn serialize_bwt_compressed(bp: &BitPop) -> IoResult<Vec<u8>> {
    let fm = match bp.get_fm_index() {
        Some(fm) => fm,
        None => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FM-Index not built",
            ))
        }
    };

    let bwt_len = fm.len();
    let mut raw = Vec::with_capacity(bwt_len);

    // BWT as raw bytes (one value per byte, values 0-4)
    for i in 0..bwt_len {
        raw.push(fm.bwt_at(i));
    }

    // zstd compress
    let compressed = zstd::encode_all(raw.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd BWT compress failed: {}", e)))?;

    let mut data = Vec::new();
    data.extend_from_slice(&(bwt_len as u64).to_le_bytes());
    data.extend_from_slice(&(compressed.len() as u64).to_le_bytes());
    data.extend_from_slice(&compressed);

    Ok(data)
}

/// Serialize SA using zstd compression.
/// Format: [sa_len: u64][zstd_compressed_sa_data]
fn serialize_sa_compressed(bp: &BitPop) -> IoResult<Vec<u8>> {
    let fm = match bp.get_fm_index() {
        Some(fm) => fm,
        None => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FM-Index not built",
            ))
        }
    };

    let sa_len = fm.sa_len();
    let mut raw = Vec::new();

    // SA entries as u32
    for rank in 0..sa_len {
        raw.extend_from_slice(&(fm.sa_at(rank) as u32).to_le_bytes());
    }

    // zstd compress
    let compressed = zstd::encode_all(raw.as_slice(), 3)
        .map_err(|e| std::io::Error::other(format!("zstd SA compress failed: {}", e)))?;

    let mut data = Vec::new();
    data.extend_from_slice(&(sa_len as u64).to_le_bytes());
    data.extend_from_slice(&(compressed.len() as u64).to_le_bytes());
    data.extend_from_slice(&compressed);

    Ok(data)
}

// --- Deserialization (load with memmap2) ---

/// Fast metadata-only load using memmap2. Returns just the header info without loading sections.
pub fn load_header(path: &str) -> IoResult<(usize, usize)> {
    let file = std::fs::File::open(path)?;
    let mmap = unsafe { Mmap::map(&file) }?;

    if mmap.len() < HEADER_SIZE + 32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "File too short",
        ));
    }

    // Validate magic
    let magic = &mmap[0..4];
    if magic != MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid magic: not BITP",
        ));
    }

    let version = u32::from_le_bytes(mmap[4..8].try_into().unwrap());
    if !(5..=VERSION).contains(&version) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Unsupported version: {} (expected 5-{})", version, VERSION),
        ));
    }

    let k = u16::from_le_bytes(mmap[8..10].try_into().unwrap()) as usize;
    let num_genomes = u32::from_le_bytes(mmap[10..14].try_into().unwrap()) as usize;

    Ok((k, num_genomes))
}

/// Load a BitPop instance from a persisted file using memmap2.
/// Format v5: BWT/SA stored uncompressed for direct memmap (<10ms load).
/// Format v4: BWT/SA stored compressed in FM_INDEX section (decompress on load).
pub fn load_bitpop(path: &str) -> IoResult<BitPop> {
    let file = std::fs::File::open(path)?;
    let mmap = unsafe { Mmap::map(&file) }?;

    if mmap.len() < HEADER_SIZE + 32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "File too short",
        ));
    }

    // Validate magic and version
    let magic = &mmap[0..4];
    if magic != MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid magic: not BITP",
        ));
    }

    let version = u32::from_le_bytes(mmap[4..8].try_into().unwrap());
    if version < 4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Unsupported version: {} (minimum 4)", version),
        ));
    }

    let k = u16::from_le_bytes(mmap[8..10].try_into().unwrap()) as usize;
    let num_genomes = u32::from_le_bytes(mmap[10..14].try_into().unwrap()) as usize;

    // Verify checksum
    let checksum_offset = mmap.len() - 32;
    let stored_checksum = &mmap[checksum_offset..];
    let mut hasher = Sha256::new();
    hasher.update(&mmap[..checksum_offset]);
    let hash_bytes = hasher.finalize();
    let computed_checksum: &[u8; 32] = hash_bytes.as_ref();

    if stored_checksum != computed_checksum.as_slice() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Checksum mismatch: file may be corrupted",
        ));
    }

    // Parse all section headers
    let num_sections = if version >= 5 { NUM_SECTIONS_V5 } else { 2 };

    let mut sections: [Option<SectionInfo>; 6] = [None, None, None, None, None, None];

    for i in 0..num_sections {
        let offset = HEADER_SIZE + (i * SECTION_HEADER_SIZE);
        let section = parse_section_header(&mmap, offset)?;

        if section.name == SECTION_BWT_UNCOMP {
            sections[0] = Some(section);
        } else if section.name == SECTION_SA_UNCOMP {
            sections[1] = Some(section);
        } else if section.name == SECTION_FM_INDEX {
            sections[2] = Some(section);
        } else if section.name == SECTION_GENOMES {
            sections[3] = Some(section);
        } else if section.name == SECTION_SPACED_SEED {
            sections[4] = Some(section);
        } else if section.name == SECTION_HF {
            sections[5] = Some(section);
        }
    }

    let genomes_section = sections[3].as_ref().ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing GENOMES section")
    })?;

    // Parse genomes from compressed GENOMES section
    let genomes_start = genomes_section.offset as usize;
    let genomes_end = genomes_start + genomes_section.compressed_size as usize;
    let genomes_compressed = &mmap[genomes_start..genomes_end];
    let genomes_decompressed = zstd::decode_all(genomes_compressed)
        .map_err(|e| std::io::Error::other(format!("Genomes decompression failed: {}", e)))?;

    let (genomes_map, genome_names_map) =
        parse_genomes_from_bytes(&genomes_decompressed, num_genomes)?;

    // Build FM-Index using v5/v6 (uncompressed memmap) or v4 (decompressed) approach
    if version >= 5 {
        // V5+ format: use uncompressed BWT/SA from memmap directly
        let bwt_section = sections[0].as_ref().ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Missing BWT_UNCOMP section",
            )
        })?;
        let sa_section = sections[1].as_ref().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing SA_UNCOMP section")
        })?;

        let fm_section = sections[2].as_ref().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing FM_INDEX section")
        })?;

        // Load BWT directly from memmap (no decompression!)
        let mut bwt = load_bwt_from_mmap(&mmap, bwt_section)?;
        let bwt_len = bwt.len();

        // Load SA directly from memmap (no decompression!)
        let sa = load_sa_from_mmap(&mmap, sa_section)?;

        // Parse C-array and boundaries from compressed FM_INDEX fallback
        let fm_start = fm_section.offset as usize;
        let fm_end = fm_start + fm_section.compressed_size as usize;
        let fm_compressed = &mmap[fm_start..fm_end];
        let fm_data = zstd::decode_all(fm_compressed)
            .map_err(|e| std::io::Error::other(format!("FM-index decompression failed: {}", e)))?;

        // Skip bwt_len (8 bytes) + bwt_packed data in FM_INDEX section
        // We already have BWT from the uncompressed section
        let bwt_len_u64 = u64::from_le_bytes(fm_data[0..8].try_into().unwrap());
        let bwt_packed_len = bwt_len_u64.div_ceil(4) as usize;
        let mut fm_pos: usize = 8 + bwt_packed_len;

        // SA_LEN + SA entries in FM_INDEX section — skip (we already have SA from uncompressed section)
        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA_LEN",
            ));
        }
        let sa_len_from_fm =
            u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8 + (sa_len_from_fm * 4); // skip SA_LEN + SA entries

        // C-array
        let mut c_array = [0usize; 5];
        for j in 0..5 {
            if fm_pos + 4 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index C_ARRAY",
                ));
            }
            c_array[j] =
                u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
            fm_pos += 4;
        }

        // Genome boundaries
        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARIES",
            ));
        }
        let num_boundaries =
            u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut genome_boundaries = Vec::with_capacity(num_boundaries);
        for _ in 0..num_boundaries {
            if fm_pos + 12 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index BOUNDARY",
                ));
            }
            let start =
                u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
            let len =
                u32::from_le_bytes(fm_data[fm_pos + 4..fm_pos + 8].try_into().unwrap()) as usize;
            let gid = u32::from_le_bytes(fm_data[fm_pos + 8..fm_pos + 12].try_into().unwrap());
            fm_pos += 12;
            genome_boundaries.push((start, len, gid));
        }

        // Sample interval (skip)
        if fm_pos + 4 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SAMPLE_INTERVAL",
            ));
        }
        fm_pos += 4;

        // Sentinel mask — needed for OccCounter to identify BWT terminators
        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at sentinel_mask_len",
            ));
        }
        let sentinel_mask_len =
            u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8;

        if fm_pos + sentinel_mask_len > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at sentinel_mask",
            ));
        }
        let sentinel_mask = &fm_data[fm_pos..fm_pos + sentinel_mask_len];

        // Mark terminators in BWT (value 4)
        for i in 0..bwt_len {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            if byte_idx < sentinel_mask.len() && (sentinel_mask[byte_idx] & (1u8 << bit_idx)) != 0 {
                bwt[i] = 0;
            }
        }

        // Build FM-Index from memmapped components
        let occ = OccCounter::new(&bwt, 32);
        let fm_index =
            FmIndex::from_components(bwt, sa, c_array, occ, genome_boundaries, num_genomes);

        let mut bp = BitPop::from_fm_index(k, genomes_map, genome_names_map, fm_index);

        // Load spaced seed section if present
        if let Some(spaced_section) = &sections[4] {
            let sp_start = spaced_section.offset as usize;
            let sp_end = sp_start + spaced_section.compressed_size as usize;
            let sp_compressed = &mmap[sp_start..sp_end];
            let sp_data = zstd::decode_all(sp_compressed).map_err(|e| {
                std::io::Error::other(format!("Spaced seed decompression failed: {}", e))
            })?;
            deserialize_spaced_seed(&sp_data, &mut bp);
        }

        // Load HF profiles section if present
        if let Some(hf_section) = &sections[5] {
            let hf_start = hf_section.offset as usize;
            let hf_end = hf_start + hf_section.compressed_size as usize;
            let hf_compressed = &mmap[hf_start..hf_end];
            let hf_data = zstd::decode_all(hf_compressed).map_err(|e| {
                std::io::Error::other(format!("HF profiles decompression failed: {}", e))
            })?;
            deserialize_hf_profiles(&hf_data, &mut bp);
        }

        Ok(bp)
    } else {
        // V4 format: decompress FM_INDEX section (original behavior)
        let fm_section = sections[2].as_ref().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing FM_INDEX section")
        })?;

        let fm_start = fm_section.offset as usize;
        let fm_end = fm_start + fm_section.compressed_size as usize;
        let fm_compressed = &mmap[fm_start..fm_end];
        let fm_data = zstd::decode_all(fm_compressed)
            .map_err(|e| std::io::Error::other(format!("FM-index decompression failed: {}", e)))?;

        // Parse FM-Index components from fm_data (same as original load_bitpop)
        let bwt_len = u64::from_le_bytes(fm_data[0..8].try_into().unwrap()) as usize;
        let mut fm_pos = 8;

        let bwt_packed_len = bwt_len.div_ceil(4);
        let bwt_data_start = fm_pos;
        fm_pos += bwt_packed_len;

        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA_LEN",
            ));
        }
        let sa_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut sa = Vec::with_capacity(sa_len);
        for _ in 0..sa_len {
            if fm_pos + 4 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index SA",
                ));
            }
            sa.push(u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize);
            fm_pos += 4;
        }

        let mut c_array = [0usize; 5];
        for j in 0..5 {
            if fm_pos + 4 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index C_ARRAY",
                ));
            }
            c_array[j] =
                u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
            fm_pos += 4;
        }

        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARIES",
            ));
        }
        let num_boundaries =
            u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut genome_boundaries = Vec::with_capacity(num_boundaries);
        for _ in 0..num_boundaries {
            if fm_pos + 12 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index BOUNDARY",
                ));
            }
            let start =
                u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
            let len =
                u32::from_le_bytes(fm_data[fm_pos + 4..fm_pos + 8].try_into().unwrap()) as usize;
            let gid = u32::from_le_bytes(fm_data[fm_pos + 8..fm_pos + 12].try_into().unwrap());
            fm_pos += 12;
            genome_boundaries.push((start, len, gid));
        }

        if fm_pos + 4 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SAMPLE_INTERVAL",
            ));
        }
        fm_pos += 4;

        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at sentinel_mask_len",
            ));
        }
        let sentinel_mask_len =
            u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
        fm_pos += 8;

        if fm_pos + sentinel_mask_len > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at sentinel_mask",
            ));
        }
        let sentinel_mask = &fm_data[fm_pos..fm_pos + sentinel_mask_len];

        let mut bwt = vec![0u8; bwt_len];
        for i in 0..bwt_len {
            let byte_idx = bwt_data_start + i / 4;
            let bit_offset = 6 - (i % 4) * 2;
            bwt[i] = ((fm_data[byte_idx] >> bit_offset) & 0x03) as u8;
        }
        for i in 0..bwt_len {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            if byte_idx < sentinel_mask.len() && (sentinel_mask[byte_idx] & (1u8 << bit_idx)) != 0 {
                bwt[i] = 0;
            }
        }

        let occ = OccCounter::new(&bwt, 32);
        let fm_index =
            FmIndex::from_components(bwt, sa, c_array, occ, genome_boundaries, num_genomes);

        Ok(BitPop::from_fm_index(
            k,
            genomes_map,
            genome_names_map,
            fm_index,
        ))
    }
}

fn parse_section_header(mmap: &Mmap, offset: usize) -> IoResult<SectionInfo> {
    let end = offset + SECTION_HEADER_SIZE;
    if end > mmap.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at section header",
        ));
    }

    let mut name = [0u8; SECTION_NAME_LEN];
    name.copy_from_slice(&mmap[offset..offset + SECTION_NAME_LEN]);

    let section_offset = u64::from_le_bytes(
        mmap[offset + SECTION_NAME_LEN..offset + SECTION_NAME_LEN + 8]
            .try_into()
            .unwrap(),
    );
    let compressed_size = u64::from_le_bytes(
        mmap[offset + SECTION_NAME_LEN + 8..offset + SECTION_NAME_LEN + 16]
            .try_into()
            .unwrap(),
    );
    let decompressed_size = u64::from_le_bytes(
        mmap[offset + SECTION_NAME_LEN + 16..offset + SECTION_NAME_LEN + 24]
            .try_into()
            .unwrap(),
    );
    let flags = u64::from_le_bytes(
        mmap[offset + SECTION_NAME_LEN + 24..offset + SECTION_NAME_LEN + 32]
            .try_into()
            .unwrap(),
    );

    Ok(SectionInfo {
        name,
        offset: section_offset,
        compressed_size,
        decompressed_size,
        flags,
    })
}

/// Load BWT directly from memmapped v5 format (no decompression).
fn load_bwt_from_mmap(mmap: &Mmap, section: &SectionInfo) -> IoResult<Vec<u8>> {
    let start = section.offset as usize;
    let end = start + section.compressed_size as usize;
    let data = &mmap[start..end];

    if data.len() < 8 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "BWT section too short",
        ));
    }

    let bwt_len = u64::from_le_bytes(data[0..8].try_into().unwrap()) as usize;

    // Detect v5 (raw) vs v7 (compressed) format using section header
    // v5: decompressed_size == compressed_size == 8 + bwt_len (raw data, no compression)
    // v7: decompressed_size == bwt_len (actual BWT length, NOT including 8-byte header)
    let is_v7 = section.decompressed_size == bwt_len as u64;

    if is_v7 {
        load_bwt_zstd(data, bwt_len)
    } else {
        // v5 backward compatibility - raw bytes
        let expected_size = 8 + bwt_len;
        if data.len() < expected_size {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "BWT data truncated",
            ));
        }
        let bwt_bytes = &data[8..expected_size];
        Ok(bwt_bytes.to_vec())
    }
}

/// Load zstd-compressed BWT from v7 format.
fn load_bwt_zstd(data: &[u8], bwt_len: usize) -> IoResult<Vec<u8>> {
    if data.len() < 16 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "BWT data too short for header",
        ));
    }

    let comp_size = u64::from_le_bytes(data[8..16].try_into().unwrap()) as usize;
    if 16 + comp_size > data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "BWT compressed data truncated",
        ));
    }
    let compressed = &data[16..16 + comp_size];

    let raw = zstd::decode_all(compressed)
        .map_err(|e| std::io::Error::other(format!("zstd BWT decompress failed: {}", e)))?;

    if raw.len() != bwt_len {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BWT decompressed size {} != expected {}",
                raw.len(),
                bwt_len
            ),
        ));
    }

    Ok(raw)
}

/// Load SA directly from memmapped v5 format (no decompression).
fn load_sa_from_mmap(mmap: &Mmap, section: &SectionInfo) -> IoResult<Vec<usize>> {
    let start = section.offset as usize;
    let end = start + section.compressed_size as usize;

    if end > mmap.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "SA section extends beyond file: offset={}, size={}, file_len={}",
                start,
                section.compressed_size,
                mmap.len()
            ),
        ));
    }

    let data = &mmap[start..end];

    if data.len() < 8 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "SA section too short",
        ));
    }

    let sa_len_u64 = u64::from_le_bytes(data[0..8].try_into().unwrap());

    // Sanity check
    if sa_len_u64 > (mmap.len() as u64) / 4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "SA length {} is impossibly large for file of {} bytes",
                sa_len_u64,
                mmap.len()
            ),
        ));
    }

    let sa_len = sa_len_u64 as usize;
    if sa_len == 0 {
        return Ok(Vec::new());
    }

    // Check if this is v6 compressed format or v5 raw format
    // v6: [sa_len: u64][comp_size: u64][zstd_data]
    // v5: [sa_len: u64][sa_entries: u32 x sa_len]
    let expected_v5_size = 8 + (sa_len * 4);
    let is_v6 = data.len() < expected_v5_size;

    if is_v6 {
        // Compressed format (v6)
        load_sa_zstd(data, sa_len)
    } else {
        // Uncompressed format (v5 backward compatibility)
        let mut sa = Vec::with_capacity(sa_len);
        for i in 0..sa_len {
            let byte_offset = 8 + (i * 4);
            sa.push(
                u32::from_le_bytes(data[byte_offset..byte_offset + 4].try_into().unwrap()) as usize,
            );
        }
        Ok(sa)
    }
}

/// Load zstd-compressed SA from v6 format.
fn load_sa_zstd(data: &[u8], sa_len: usize) -> IoResult<Vec<usize>> {
    if data.len() < 16 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "SA data too short for header",
        ));
    }

    let comp_size = u64::from_le_bytes(data[8..16].try_into().unwrap()) as usize;
    let compressed = &data[16..16 + comp_size];

    let raw = zstd::decode_all(compressed)
        .map_err(|e| std::io::Error::other(format!("zstd SA decompress failed: {}", e)))?;

    if raw.len() != sa_len * 4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "SA decompressed size {} != expected {}",
                raw.len(),
                sa_len * 4
            ),
        ));
    }

    let mut sa = Vec::with_capacity(sa_len);
    for i in 0..sa_len {
        let off = i * 4;
        sa.push(u32::from_le_bytes(raw[off..off + 4].try_into().unwrap()) as usize);
    }

    Ok(sa)
}

/// Parse genomes from decompressed bytes, returning (genomes_map, names_map).
fn parse_genomes_from_bytes(data: &[u8], num_genomes: usize) -> IoResult<GenomeData> {
    let mut genomes: HashMap<u32, Vec<u8>> = HashMap::new();
    let mut genome_names: HashMap<u32, String> = HashMap::new();

    let mut gpos = 0;
    for i in 0..num_genomes {
        if gpos + 4 > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name_len", i),
            ));
        }
        let name_len = u32::from_le_bytes(data[gpos..gpos + 4].try_into().unwrap()) as usize;
        gpos += 4;

        if gpos + name_len > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name", i),
            ));
        }
        let name_str = String::from_utf8(data[gpos..gpos + name_len].to_vec())
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        gpos += name_len;

        if gpos + 8 > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq_len", i),
            ));
        }
        let seq_len = u64::from_le_bytes(data[gpos..gpos + 8].try_into().unwrap()) as usize;
        gpos += 8;

        if gpos + seq_len > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq data", i),
            ));
        }
        let seq_bytes = data[gpos..gpos + seq_len].to_vec();
        gpos += seq_len;

        let gid = i as u32;
        genome_names.insert(gid, name_str);
        genomes.insert(gid, seq_bytes);
    }

    Ok((genomes, genome_names))
}

// --- Legacy format support (PLAN2 format) ---

/// Load from the old PLAN2 format (single zstd block, no memmap).
/// Kept for backward compatibility.
pub fn load_legacy_bitpop(path: &str) -> IoResult<BitPop> {
    let mut file = std::fs::File::open(path)?;
    let mut raw_data = Vec::new();
    file.read_to_end(&mut raw_data)?;

    if raw_data.len() < HEADER_SIZE + 32 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "File too short",
        ));
    }

    let magic = &raw_data[0..4];
    if magic != MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid magic: not BITP",
        ));
    }

    let version = u32::from_le_bytes(raw_data[4..8].try_into().unwrap());
    if version != 3 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Unsupported legacy version: {}", version),
        ));
    }

    let k = u16::from_le_bytes(raw_data[8..10].try_into().unwrap()) as usize;
    let num_genomes = u32::from_le_bytes(raw_data[10..14].try_into().unwrap()) as usize;

    let compressed_start = HEADER_SIZE;
    let checksum_start = raw_data.len() - 32;
    let compressed = &raw_data[compressed_start..checksum_start];

    let decompressed = zstd::decode_all(compressed)
        .map_err(|e| std::io::Error::other(format!("zstd decompression failed: {}", e)))?;

    // Parse outer all_data to get fm_compressed and genomes_compressed
    if decompressed.len() < 16 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "All-data too short",
        ));
    }
    let fm_compressed_len = u64::from_le_bytes(decompressed[0..8].try_into().unwrap()) as usize;
    let fm_compressed = &decompressed[8..8 + fm_compressed_len];
    let genomes_compressed_offset = 8 + fm_compressed_len;

    let fm_data = zstd::decode_all(fm_compressed)
        .map_err(|e| std::io::Error::other(format!("FM-index decompression failed: {}", e)))?;

    // Parse FM-Index (same as in load_bitpop)
    let bwt_len = u64::from_le_bytes(fm_data[0..8].try_into().unwrap()) as usize;
    let mut fm_pos = 8;
    let bwt_packed_len = bwt_len.div_ceil(4);
    let bwt_data_start = fm_pos;
    fm_pos += bwt_packed_len;

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index SA_LEN",
        ));
    }
    let sa_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
    fm_pos += 8;

    let mut sa = Vec::with_capacity(sa_len);
    for _ in 0..sa_len {
        if fm_pos + 4 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA",
            ));
        }
        sa.push(u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize);
        fm_pos += 4;
    }

    let mut c_array = [0usize; 5];
    for j in 0..5 {
        if fm_pos + 4 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index C_ARRAY",
            ));
        }
        c_array[j] = u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
        fm_pos += 4;
    }

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index BOUNDARIES",
        ));
    }
    let num_boundaries =
        u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
    fm_pos += 8;

    let mut genome_boundaries = Vec::with_capacity(num_boundaries);
    for _ in 0..num_boundaries {
        if fm_pos + 12 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARY",
            ));
        }
        let start = u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap()) as usize;
        let len = u32::from_le_bytes(fm_data[fm_pos + 4..fm_pos + 8].try_into().unwrap()) as usize;
        let gid = u32::from_le_bytes(fm_data[fm_pos + 8..fm_pos + 12].try_into().unwrap());
        fm_pos += 12;
        genome_boundaries.push((start, len, gid));
    }

    if fm_pos + 4 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index SAMPLE_INTERVAL",
        ));
    }
    let _sample_interval = u32::from_le_bytes(fm_data[fm_pos..fm_pos + 4].try_into().unwrap());
    fm_pos += 4;

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at sentinel_mask_len",
        ));
    }
    let sentinel_mask_len =
        u64::from_le_bytes(fm_data[fm_pos..fm_pos + 8].try_into().unwrap()) as usize;
    fm_pos += 8;

    if fm_pos + sentinel_mask_len > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at sentinel_mask",
        ));
    }
    let sentinel_mask = &fm_data[fm_pos..fm_pos + sentinel_mask_len];

    let mut bwt = vec![0u8; bwt_len];
    for i in 0..bwt_len {
        let byte_idx = bwt_data_start + i / 4;
        let bit_offset = 6 - (i % 4) * 2;
        bwt[i] = ((fm_data[byte_idx] >> bit_offset) & 0x03) as u8;
    }
    for i in 0..bwt_len {
        let byte_idx = i / 8;
        let bit_idx = i % 8;
        if byte_idx < sentinel_mask.len() && (sentinel_mask[byte_idx] & (1u8 << bit_idx)) != 0 {
            bwt[i] = 0;
        }
    }

    let occ = OccCounter::new(&bwt, 32);
    let fm_index = FmIndex::from_components(bwt, sa, c_array, occ, genome_boundaries, num_genomes);

    let mut genomes: std::collections::HashMap<u32, Vec<u8>> = std::collections::HashMap::new();
    let mut genome_names: std::collections::HashMap<u32, String> = std::collections::HashMap::new();

    if genomes_compressed_offset + 16 > decompressed.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at genomes compressed length",
        ));
    }
    let genomes_compressed_len = u64::from_le_bytes(
        decompressed[genomes_compressed_offset..genomes_compressed_offset + 8]
            .try_into()
            .unwrap(),
    ) as usize;
    let genomes_compressed = &decompressed
        [genomes_compressed_offset + 8..genomes_compressed_offset + 8 + genomes_compressed_len];

    let genomes_decompressed = zstd::decode_all(genomes_compressed)
        .map_err(|e| std::io::Error::other(format!("Genomes decompression failed: {}", e)))?;

    let mut gpos = 0;
    for i in 0..num_genomes {
        if gpos + 4 > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name_len", i),
            ));
        }
        let name_len =
            u32::from_le_bytes(genomes_decompressed[gpos..gpos + 4].try_into().unwrap()) as usize;
        gpos += 4;

        if gpos + name_len > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name", i),
            ));
        }
        let name_str = String::from_utf8(genomes_decompressed[gpos..gpos + name_len].to_vec())
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        gpos += name_len;

        if gpos + 8 > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq_len", i),
            ));
        }
        let seq_len =
            u64::from_le_bytes(genomes_decompressed[gpos..gpos + 8].try_into().unwrap()) as usize;
        gpos += 8;

        if gpos + seq_len > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq data", i),
            ));
        }
        let seq_bytes = genomes_decompressed[gpos..gpos + seq_len].to_vec();
        gpos += seq_len;

        let gid = i as u32;
        genome_names.insert(gid, name_str);
        genomes.insert(gid, seq_bytes);
    }

    Ok(BitPop::from_fm_index(k, genomes, genome_names, fm_index))
}

// --- Auto-detect format and load ---

/// Load a BitPop instance, auto-detecting the format (new memmap or legacy).
pub fn load_bitpop_auto(path: &str) -> IoResult<BitPop> {
    // Try new format first
    if let Ok(meta) = std::fs::metadata(path) {
        if meta.len() >= (HEADER_SIZE + 32) as u64 {
            // Peek at the version field to detect format
            let mut file = std::fs::File::open(path)?;
            let mut header_buf = [0u8; 14];
            if file.read_exact(&mut header_buf).is_ok() && header_buf[0..4] == MAGIC[..] {
                let version = u32::from_le_bytes(header_buf[4..8].try_into().unwrap());
                if version == VERSION {
                    return load_bitpop(path); // new format with memmap2
                } else if version == 3 {
                    return load_legacy_bitpop(path); // legacy format
                }
            }
        }
    }

    // Fallback: try legacy
    load_legacy_bitpop(path)
}

/// Compute SHA256 hash of data.
#[allow(dead_code)]
fn sha256(data: &[u8]) -> [u8; 32] {
    use sha2::Digest;
    let mut hasher = sha2::Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}

// --- Tests ---

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
    fn test_persisted_roundtrip() {
        let bp = make_test_bitpop();
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_v4_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap();

        save_bitpop(&bp, path_str).unwrap();

        // Test header-only load (memmap2)
        let (k, num_genomes) = load_header(path_str).unwrap();
        assert_eq!(k, 6);
        assert_eq!(num_genomes, 3);

        // Test full load
        let loaded = load_bitpop(path_str).unwrap();
        assert_eq!(loaded.k(), 6);
        assert_eq!(loaded.genome_count(), 3);
        assert_eq!(loaded.genome_name(0), Some("human"));
        assert_eq!(loaded.genome_name(1), Some("chimp"));
        assert_eq!(loaded.genome_name(2), Some("mouse"));

        let results_orig = bp.map_read("ACGTACGT", 3);
        let results_loaded = loaded.map_read("ACGTACGT", 3);
        assert!(!results_orig.is_empty());
        assert!(!results_loaded.is_empty());
        assert_eq!(results_loaded[0].genome_id, results_orig[0].genome_id);

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_persisted_single_genome() {
        let mut bp = BitPop::new(8);
        bp.add_genome("test", "AACCGGTTAACCGGTT");
        bp.build();

        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_v4_single_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap();

        save_bitpop(&bp, path_str).unwrap();
        let loaded = load_bitpop(path_str).unwrap();

        assert_eq!(loaded.genome_count(), 1);
        assert_eq!(loaded.genome_name(0), Some("test"));

        let results = loaded.map_read("AACCGGTT", 3);
        assert!(!results.is_empty());

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_persisted_invalid_magic() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_bad_magic_v4_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap();

        std::fs::write(path_str, b"BADX").unwrap();

        // Header load should fail
        let result = load_header(path_str);
        assert!(result.is_err());

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_persisted_compression_ratio() {
        let mut bp = BitPop::new(8);
        let genome = format!(
            "{}{}{}",
            "ACGT".repeat(5000),
            "AACCGGTT",
            "TTTT".repeat(5000)
        );
        bp.add_genome("large", &genome);
        bp.build();

        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_v4_comp_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap();

        save_bitpop(&bp, path_str).unwrap();
        let file_size = std::fs::metadata(path_str).unwrap().len();

        // v5 format: stores both uncompressed BWT/SA (for memmap) + compressed FM_INDEX/GENOMES
        // 40K bases -> ~170KB SA+BWT uncompressed + ~20KB compressed sections = ~200KB total
        assert!(
            file_size < 300000,
            "Compressed file {} bytes is too large",
            file_size
        );

        let loaded = load_bitpop(path_str).unwrap();
        assert_eq!(loaded.genome_seq_len(0), Some(genome.len()));

        let results = loaded.map_read("AACCGGTT", 3);
        assert!(!results.is_empty());

        let _ = std::fs::remove_file(path_str);
    }

    #[test]
    fn test_legacy_format_compatibility() {
        // Load a legacy (v3) file using the auto-detect function
        // First create a v3 format file manually
        let bp = make_test_bitpop();
        let dir = std::env::temp_dir();

        // Create a v3 format file by writing data in the old nested format
        let legacy_path = dir.join(format!(
            "bitpop_v3_compat_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let legacy_path_str = legacy_path.to_str().unwrap();

        // Build v3 format data manually
        let mut header = Vec::with_capacity(64);
        header.extend_from_slice(b"BITP");
        header.extend_from_slice(&3u32.to_le_bytes()); // version 3
        header.extend_from_slice(&6u16.to_le_bytes()); // k=6
        header.extend_from_slice(&(3u32).to_le_bytes()); // 3 genomes
        header.resize(64, 0u8);

        // Build nested compressed data (v3 format)
        let mut genomes_data = Vec::new();
        for name in &["human", "chimp", "mouse"] {
            let seq = match *name {
                "human" => b"ACGTACGTACGTACGTACGTACGT".to_vec(),
                "chimp" => b"ACGTACGTACGTAACAACGTACGT".to_vec(),
                "mouse" => b"TTTTGGGGACGTACGTACGTACGT".to_vec(),
                _ => vec![],
            };
            genomes_data.extend_from_slice(&(name.len() as u32).to_le_bytes());
            genomes_data.extend_from_slice(name.as_bytes());
            genomes_data.extend_from_slice(&(seq.len() as u64).to_le_bytes());
            genomes_data.extend_from_slice(&seq);
        }

        let fm = bp.get_fm_index().unwrap();
        let bwt_len = fm.len();
        let mut fm_data = Vec::new();
        fm_data.extend_from_slice(&(bwt_len as u64).to_le_bytes());

        let bwt_packed_len = bwt_len.div_ceil(4);
        let mut bwt_packed = vec![0u8; bwt_packed_len];
        for i in 0..bwt_len {
            let v = fm.bwt_at(i) & 3;
            let byte_idx = i / 4;
            let bit_offset = 6 - (i % 4) * 2;
            bwt_packed[byte_idx] |= v << bit_offset;
        }
        fm_data.extend_from_slice(&bwt_packed);

        let sa_len = fm.sa_len();
        fm_data.extend_from_slice(&(sa_len as u64).to_le_bytes());
        for rank in 0..sa_len {
            fm_data.extend_from_slice(&(fm.sa_at(rank) as u32).to_le_bytes());
        }
        for j in 0..5 {
            fm_data.extend_from_slice(&(fm.c_array(j) as u32).to_le_bytes());
        }
        let boundaries = fm.genome_boundaries();
        fm_data.extend_from_slice(&(boundaries.len() as u64).to_le_bytes());
        for &(start, len, gid) in boundaries {
            fm_data.extend_from_slice(&(start as u32).to_le_bytes());
            fm_data.extend_from_slice(&(len as u32).to_le_bytes());
            fm_data.extend_from_slice(&gid.to_le_bytes());
        }
        fm_data.extend_from_slice(&32u32.to_le_bytes());

        let sentinel_mask_len = bwt_len.div_ceil(8);
        let mut sentinel_mask = vec![0u8; sentinel_mask_len];
        for i in 0..bwt_len {
            if fm.bwt_at(i) == 0 {
                let byte_idx = i / 8;
                let bit_idx = i % 8;
                sentinel_mask[byte_idx] |= 1u8 << bit_idx;
            }
        }
        fm_data.extend_from_slice(&(sentinel_mask_len as u64).to_le_bytes());
        fm_data.extend_from_slice(&sentinel_mask);

        let fm_compressed = zstd::encode_all(fm_data.as_slice(), 3).unwrap();
        let genomes_compressed = zstd::encode_all(genomes_data.as_slice(), 3).unwrap();

        let mut all_data = Vec::new();
        all_data.extend_from_slice(&(fm_compressed.len() as u64).to_le_bytes());
        all_data.extend_from_slice(&fm_compressed);
        all_data.extend_from_slice(&(genomes_compressed.len() as u64).to_le_bytes());
        all_data.extend_from_slice(&genomes_compressed);

        let nested = zstd::encode_all(all_data.as_slice(), 3).unwrap();

        let mut file_data = Vec::new();
        file_data.extend_from_slice(&header);
        file_data.extend_from_slice(&nested);
        let checksum = sha256(&file_data);
        file_data.extend_from_slice(&checksum);

        std::fs::write(legacy_path_str, &file_data).unwrap();

        // Now load with auto-detect — should use legacy loader
        let loaded = load_bitpop_auto(legacy_path_str).unwrap();
        assert_eq!(loaded.k(), 6);
        assert_eq!(loaded.genome_count(), 3);
        assert_eq!(loaded.genome_name(0), Some("human"));

        let _ = std::fs::remove_file(legacy_path_str);
    }

    #[test]
    fn test_new_format_after_legacy() {
        // Create a v4 file, verify auto-detect uses new loader
        let bp = make_test_bitpop();
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "bitpop_v4_auto_{}_{}.bitpop",
            std::process::id(),
            std::time::SystemTime::now().elapsed().unwrap().as_nanos()
        ));
        let path_str = path.to_str().unwrap();

        save_bitpop(&bp, path_str).unwrap();

        // Auto-detect should use new memmap format
        let loaded = load_bitpop_auto(path_str).unwrap();
        assert_eq!(loaded.genome_count(), 3);

        let results = loaded.map_read("ACGTACGT", 3);
        assert!(!results.is_empty());

        let _ = std::fs::remove_file(path_str);
    }
}
