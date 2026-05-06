use std::io::{Read, Result as IoResult};

use memmap2::Mmap;
use sha2::{Digest, Sha256};

use crate::fm::{FmIndex, OccCounter};
use crate::BitPop;

// --- File Format Constants ---

const MAGIC: [u8; 4] = *b"BITP";
const VERSION: u32 = 5;
const HEADER_SIZE: usize = 64;
const SECTION_NAME_LEN: usize = 16;
const SECTION_HEADER_SIZE: usize = 48; // name[16] + offset(8) + comp_size(8) + decomp_size(8) + flags(8)

// Section names (padded to SECTION_NAME_LEN=16)
const SECTION_BWT_UNCOMP: [u8; 16] = *b"BWT_UNCOMP\0\0\0\0\0\0";
const SECTION_SA_UNCOMP: [u8; 16] = *b"SA_UNCOMP\0\0\0\0\0\0\0";
const SECTION_FM_INDEX: [u8; 16] = *b"FM_INDEX\0\0\0\0\0\0\0\0";
const SECTION_GENOMES: [u8; 16] = *b"GENOMES\0\0\0\0\0\0\0\0\0";

/// Number of sections in v5+ format (BWT_UNCOMP + SA_UNCOMP + FM_INDEX + GENOMES)
const NUM_SECTIONS_V5: usize = 4;

/// Represents a section in the persisted file.
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
    magic: [u8; 4],       // "BITP"
    version: u32,         // format version
    k: u16,               // k-mer size
    num_genomes: u32,     // number of genomes
    _reserved: [u8; 46],  // padding to 64 bytes
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

/// Save a BitPop instance to a file using the memmap-friendly format.
/// Format v5: [header][section_table][BWT_UNCOMP][SA_UNCOMP][FM_INDEX][GENOMES][checksum]
/// BWT_UNCOMP and SA_UNCOMP are stored uncompressed for direct memmap access (<10ms load).
pub fn save_bitpop(bp: &BitPop, path: &str) -> IoResult<()> {
    // 1. Serialize FM-Index (compressed fallback)
    let fm_data = serialize_fm_index(bp)?;
    let fm_compressed = zstd::encode_all(fm_data.as_slice(), 3).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("zstd FM compress failed: {}", e))
    })?;

    // 2. Serialize genomes (compressed)
    let genomes_data = serialize_genomes(bp)?;
    let genomes_compressed = zstd::encode_all(genomes_data.as_slice(), 3).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("zstd genomes compress failed: {}", e))
    })?;

    // 3. Serialize BWT uncompressed (4 values per byte, 2 bits each)
    let bwt_uncomp = serialize_bwt_uncompressed(bp)?;

    // 4. Serialize SA uncompressed (u32 per entry)
    let sa_uncomp = serialize_sa_uncompressed(bp)?;

    // 5. Build section table (4 sections: BWT_UNCOMP, SA_UNCOMP, FM_INDEX, GENOMES)
    let mut section_table = Vec::new();
    
    let base_offset: u64 = (HEADER_SIZE + (NUM_SECTIONS_V5 * SECTION_HEADER_SIZE)) as u64;
    
    let mut offset = base_offset;
    write_section_header(&mut section_table, &SECTION_BWT_UNCOMP, offset, bwt_uncomp.len() as u64, bwt_uncomp.len() as u64, 0);
    offset += bwt_uncomp.len() as u64;
    
    write_section_header(&mut section_table, &SECTION_SA_UNCOMP, offset, sa_uncomp.len() as u64, sa_uncomp.len() as u64, 0);
    offset += sa_uncomp.len() as u64;
    
    write_section_header(&mut section_table, &SECTION_FM_INDEX, offset, fm_compressed.len() as u64, fm_data.len() as u64, 0);
    offset += fm_compressed.len() as u64;
    
    write_section_header(&mut section_table, &SECTION_GENOMES, offset, genomes_compressed.len() as u64, genomes_data.len() as u64, 0);

    // 6. Assemble file: header + section_table + sections
    let mut all_data = Vec::new();
    all_data.reserve((offset + genomes_compressed.len() as u64 + 32) as usize);
    
    let header_placeholder = vec![0u8; HEADER_SIZE];
    all_data.extend_from_slice(&header_placeholder);
    all_data.extend_from_slice(&section_table);
    all_data.extend_from_slice(&bwt_uncomp);
    all_data.extend_from_slice(&sa_uncomp);
    all_data.extend_from_slice(&fm_compressed);
    all_data.extend_from_slice(&genomes_compressed);

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

fn write_section_header(buf: &mut Vec<u8>, name: &[u8], offset: u64, comp_size: u64, decomp_size: u64, _flags: u64) {
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
        None => return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "FM-Index not built",
        )),
    };

    let mut data = Vec::new();
    let bwt_len = fm.len();

    // BWT length
    data.extend_from_slice(&(bwt_len as u64).to_le_bytes());

    // BWT: 2 bits per entry (4 values per byte)
    let bwt_packed_len = (bwt_len + 3) / 4;
    let mut bwt_packed = vec![0u8; bwt_packed_len];
    for i in 0..bwt_len {
        let v = (fm.bwt_at(i) & 3) as u8;
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
        data.extend_from_slice(&(gid as u32).to_le_bytes());
    }

    // Sample interval
    data.extend_from_slice(&32u32.to_le_bytes());

    // Sentinel mask
    let sentinel_mask_len = (bwt_len + 7) / 8;
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

/// Serialize BWT in uncompressed raw format for direct memmap.
/// Format: [bwt_len: u64][bwt_bytes: u8 x bwt_len, one 2-bit value per byte]
fn serialize_bwt_uncompressed(bp: &BitPop) -> IoResult<Vec<u8>> {
    let fm = match bp.get_fm_index() {
        Some(fm) => fm,
        None => return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "FM-Index not built",
        )),
    };

    let bwt_len = fm.len();
    let mut data = Vec::new();
    
    // Length header
    data.extend_from_slice(&(bwt_len as u64).to_le_bytes());
    
    // BWT as raw bytes (one value per byte, values 0-4)
    for i in 0..bwt_len {
        data.push(fm.bwt_at(i) as u8);
    }

    Ok(data)
}

/// Serialize SA in uncompressed u32 format for direct memmap.
/// Format: [sa_len: u64][sa_entries: u32 x sa_len]
fn serialize_sa_uncompressed(bp: &BitPop) -> IoResult<Vec<u8>> {
    let fm = match bp.get_fm_index() {
        Some(fm) => fm,
        None => return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "FM-Index not built",
        )),
    };

    let sa_len = fm.sa_len();
    let mut data = Vec::new();
    
    // Length header
    data.extend_from_slice(&(sa_len as u64).to_le_bytes());
    
    // SA entries as u32
    for rank in 0..sa_len {
        data.extend_from_slice(&(fm.sa_at(rank) as u32).to_le_bytes());
    }

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
    if version != VERSION {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Unsupported version: {} (expected {})", version, VERSION),
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
    let num_sections = if version >= VERSION { NUM_SECTIONS_V5 } else { 2 };
    
    let mut sections: [Option<SectionInfo>; 4] = [None, None, None, None];

    for i in 0..num_sections {
        let offset = HEADER_SIZE + (i * SECTION_HEADER_SIZE);
        let section = parse_section_header(&mmap, offset)?;
        
        if section.name == SECTION_BWT_UNCOMP { sections[0] = Some(section); }
        else if section.name == SECTION_SA_UNCOMP { sections[1] = Some(section); }
        else if section.name == SECTION_FM_INDEX { sections[2] = Some(section); }
        else if section.name == SECTION_GENOMES { sections[3] = Some(section); }
    }

    let genomes_section = sections[3].as_ref().ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing GENOMES section")
    })?;

    // Parse genomes from compressed GENOMES section
    let genomes_start = genomes_section.offset as usize;
    let genomes_end = genomes_start + genomes_section.compressed_size as usize;
    let genomes_compressed = &mmap[genomes_start..genomes_end];
    let genomes_decompressed = zstd::decode_all(genomes_compressed).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("Genomes decompression failed: {}", e))
    })?;

    let (genomes_map, genome_names_map) = parse_genomes_from_bytes(&genomes_decompressed, num_genomes)?;

    // Build FM-Index using v5 (uncompressed memmap) or v4 (decompressed) approach
    if version >= VERSION {
        // V5+ format: use uncompressed BWT/SA from memmap directly
        let bwt_section = sections[0].as_ref().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing BWT_UNCOMP section")
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
        let fm_data = zstd::decode_all(fm_compressed).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::Other, format!("FM-index decompression failed: {}", e))
        })?;

        // Skip bwt_len (8 bytes) + bwt_packed data in FM_INDEX section
        // We already have BWT from the uncompressed section
        let bwt_len_u64 = u64::from_le_bytes(fm_data[0..8].try_into().unwrap());
        let bwt_packed_len = ((bwt_len_u64 + 3) / 4) as usize;
        let mut fm_pos: usize = 8 + bwt_packed_len;

        // SA_LEN + SA entries in FM_INDEX section — skip (we already have SA from uncompressed section)
        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA_LEN",
            ));
        }
        let sa_len_from_fm = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
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
            c_array[j] = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
            fm_pos += 4;
        }

        // Genome boundaries
        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARIES",
            ));
        }
        let num_boundaries = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut genome_boundaries = Vec::with_capacity(num_boundaries);
        for _ in 0..num_boundaries {
            if fm_pos + 12 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index BOUNDARY",
                ));
            }
            let start = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
            let len = u32::from_le_bytes(fm_data[fm_pos+4..fm_pos+8].try_into().unwrap()) as usize;
            let gid = u32::from_le_bytes(fm_data[fm_pos+8..fm_pos+12].try_into().unwrap());
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
        let sentinel_mask_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
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
        let fm_index = FmIndex::from_components(
            bwt, sa, c_array, occ, genome_boundaries, num_genomes,
        );

        Ok(BitPop::from_fm_index(k, genomes_map, genome_names_map, fm_index))
    } else {
        // V4 format: decompress FM_INDEX section (original behavior)
        let fm_section = sections[2].as_ref().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing FM_INDEX section")
        })?;

        let fm_start = fm_section.offset as usize;
        let fm_end = fm_start + fm_section.compressed_size as usize;
        let fm_compressed = &mmap[fm_start..fm_end];
        let fm_data = zstd::decode_all(fm_compressed).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::Other, format!("FM-index decompression failed: {}", e))
        })?;

        // Parse FM-Index components from fm_data (same as original load_bitpop)
        let bwt_len = u64::from_le_bytes(fm_data[0..8].try_into().unwrap()) as usize;
        let mut fm_pos = 8;

        let bwt_packed_len = (bwt_len + 3) / 4;
        let bwt_data_start = fm_pos;
        fm_pos += bwt_packed_len;

        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA_LEN",
            ));
        }
        let sa_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut sa = Vec::with_capacity(sa_len);
        for _ in 0..sa_len {
            if fm_pos + 4 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index SA",
                ));
            }
            sa.push(u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize);
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
            c_array[j] = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
            fm_pos += 4;
        }

        if fm_pos + 8 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARIES",
            ));
        }
        let num_boundaries = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
        fm_pos += 8;

        let mut genome_boundaries = Vec::with_capacity(num_boundaries);
        for _ in 0..num_boundaries {
            if fm_pos + 12 > fm_data.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Unexpected end at FM-Index BOUNDARY",
                ));
            }
            let start = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
            let len = u32::from_le_bytes(fm_data[fm_pos+4..fm_pos+8].try_into().unwrap()) as usize;
            let gid = u32::from_le_bytes(fm_data[fm_pos+8..fm_pos+12].try_into().unwrap());
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
        let sentinel_mask_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
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
        let fm_index = FmIndex::from_components(
            bwt, sa, c_array, occ, genome_boundaries, num_genomes,
        );

        Ok(BitPop::from_fm_index(k, genomes_map, genome_names_map, fm_index))
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

    let section_offset = u64::from_le_bytes(mmap[offset + SECTION_NAME_LEN..offset + SECTION_NAME_LEN + 8].try_into().unwrap());
    let compressed_size = u64::from_le_bytes(mmap[offset + SECTION_NAME_LEN + 8..offset + SECTION_NAME_LEN + 16].try_into().unwrap());
    let decompressed_size = u64::from_le_bytes(mmap[offset + SECTION_NAME_LEN + 16..offset + SECTION_NAME_LEN + 24].try_into().unwrap());
    let flags = u64::from_le_bytes(mmap[offset + SECTION_NAME_LEN + 24..offset + SECTION_NAME_LEN + 32].try_into().unwrap());

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
    let end = start + section.compressed_size as usize; // compressed_size == uncompressed size for v5 BWT
    let data = &mmap[start..end];

    if data.len() < 8 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "BWT section too short",
        ));
    }

    let bwt_len = u64::from_le_bytes(data[0..8].try_into().unwrap()) as usize;

    if 8 + bwt_len > data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "BWT data truncated",
        ));
    }

    let bwt_bytes = &data[8..8 + bwt_len];
    Ok(bwt_bytes.to_vec())
}

/// Load SA directly from memmapped v5 format (no decompression).
fn load_sa_from_mmap(mmap: &Mmap, section: &SectionInfo) -> IoResult<Vec<usize>> {
    let start = section.offset as usize;
    let end = start + section.compressed_size as usize; // compressed_size == uncompressed size for v5 SA
    
    if end > mmap.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("SA section extends beyond file: offset={}, size={}, file_len={}", start, section.compressed_size, mmap.len()),
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
    
    // Sanity check: SA length should not exceed file size / 4
    if sa_len_u64 > (mmap.len() as u64) / 4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("SA length {} is impossibly large for file of {} bytes", sa_len_u64, mmap.len()),
        ));
    }
    
    let sa_len = sa_len_u64 as usize;
    let data_start = 8;
    let expected_size = data_start + (sa_len * 4);

    if expected_size > data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("SA data truncated: need {} bytes, have {}", expected_size, data.len()),
        ));
    }

    let mut sa = Vec::with_capacity(sa_len);
    for i in 0..sa_len {
        let byte_offset = data_start + (i * 4);
        sa.push(u32::from_le_bytes(
            data[byte_offset..byte_offset + 4].try_into().unwrap()
        ) as usize);
    }

    Ok(sa)
}

/// Parse genomes from decompressed bytes, returning (genomes_map, names_map).
fn parse_genomes_from_bytes(
    data: &[u8],
    num_genomes: usize,
) -> IoResult<(std::collections::HashMap<u32, Vec<u8>>, std::collections::HashMap<u32, String>)> {
    let mut genomes: std::collections::HashMap<u32, Vec<u8>> = std::collections::HashMap::new();
    let mut genome_names: std::collections::HashMap<u32, String> = std::collections::HashMap::new();

    let mut gpos = 0;
    for i in 0..num_genomes {
        if gpos + 4 > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name_len", i),
            ));
        }
        let name_len = u32::from_le_bytes(data[gpos..gpos+4].try_into().unwrap()) as usize;
        gpos += 4;

        if gpos + name_len > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name", i),
            ));
        }
        let name_str = String::from_utf8(data[gpos..gpos+name_len].to_vec()).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, e)
        })?;
        gpos += name_len;

        if gpos + 8 > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq_len", i),
            ));
        }
        let seq_len = u64::from_le_bytes(data[gpos..gpos+8].try_into().unwrap()) as usize;
        gpos += 8;

        if gpos + seq_len > data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq data", i),
            ));
        }
        let seq_bytes = data[gpos..gpos+seq_len].to_vec();
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

    let decompressed = zstd::decode_all(compressed).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("zstd decompression failed: {}", e))
    })?;

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

    let fm_data = zstd::decode_all(fm_compressed).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("FM-index decompression failed: {}", e))
    })?;

    // Parse FM-Index (same as in load_bitpop)
    let bwt_len = u64::from_le_bytes(fm_data[0..8].try_into().unwrap()) as usize;
    let mut fm_pos = 8;
    let bwt_packed_len = (bwt_len + 3) / 4;
    let bwt_data_start = fm_pos;
    fm_pos += bwt_packed_len;

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index SA_LEN",
        ));
    }
    let sa_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
    fm_pos += 8;

    let mut sa = Vec::with_capacity(sa_len);
    for _ in 0..sa_len {
        if fm_pos + 4 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index SA",
            ));
        }
        sa.push(u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize);
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
        c_array[j] = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
        fm_pos += 4;
    }

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index BOUNDARIES",
        ));
    }
    let num_boundaries = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
    fm_pos += 8;

    let mut genome_boundaries = Vec::with_capacity(num_boundaries);
    for _ in 0..num_boundaries {
        if fm_pos + 12 > fm_data.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected end at FM-Index BOUNDARY",
            ));
        }
        let start = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap()) as usize;
        let len = u32::from_le_bytes(fm_data[fm_pos+4..fm_pos+8].try_into().unwrap()) as usize;
        let gid = u32::from_le_bytes(fm_data[fm_pos+8..fm_pos+12].try_into().unwrap());
        fm_pos += 12;
        genome_boundaries.push((start, len, gid));
    }

    if fm_pos + 4 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at FM-Index SAMPLE_INTERVAL",
        ));
    }
    let _sample_interval = u32::from_le_bytes(fm_data[fm_pos..fm_pos+4].try_into().unwrap());
    fm_pos += 4;

    if fm_pos + 8 > fm_data.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at sentinel_mask_len",
        ));
    }
    let sentinel_mask_len = u64::from_le_bytes(fm_data[fm_pos..fm_pos+8].try_into().unwrap()) as usize;
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
    let fm_index = FmIndex::from_components(
        bwt, sa, c_array, occ, genome_boundaries, num_genomes,
    );

    let mut genomes: std::collections::HashMap<u32, Vec<u8>> = std::collections::HashMap::new();
    let mut genome_names: std::collections::HashMap<u32, String> = std::collections::HashMap::new();

    if genomes_compressed_offset + 16 > decompressed.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Unexpected end at genomes compressed length",
        ));
    }
    let genomes_compressed_len = u64::from_le_bytes(decompressed[genomes_compressed_offset..genomes_compressed_offset+8].try_into().unwrap()) as usize;
    let genomes_compressed = &decompressed[genomes_compressed_offset+8..genomes_compressed_offset+8+genomes_compressed_len];

    let genomes_decompressed = zstd::decode_all(genomes_compressed).map_err(|e| {
        std::io::Error::new(std::io::ErrorKind::Other, format!("Genomes decompression failed: {}", e))
    })?;

    let mut gpos = 0;
    for i in 0..num_genomes {
        if gpos + 4 > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name_len", i),
            ));
        }
        let name_len = u32::from_le_bytes(genomes_decompressed[gpos..gpos+4].try_into().unwrap()) as usize;
        gpos += 4;

        if gpos + name_len > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} name", i),
            ));
        }
        let name_str = String::from_utf8(genomes_decompressed[gpos..gpos+name_len].to_vec()).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, e)
        })?;
        gpos += name_len;

        if gpos + 8 > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq_len", i),
            ));
        }
        let seq_len = u64::from_le_bytes(genomes_decompressed[gpos..gpos+8].try_into().unwrap()) as usize;
        gpos += 8;

        if gpos + seq_len > genomes_decompressed.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unexpected end at genome {} seq data", i),
            ));
        }
        let seq_bytes = genomes_decompressed[gpos..gpos+seq_len].to_vec();
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
    match std::fs::metadata(path) {
        Ok(meta) => {
            if meta.len() as u64 >= (HEADER_SIZE + 32) as u64 {
                // Peek at the version field to detect format
                let mut file = std::fs::File::open(path)?;
                let mut header_buf = [0u8; 14];
                if file.read_exact(&mut header_buf).is_ok() {
                    if &header_buf[0..4] == &MAGIC[..] {
                        let version = u32::from_le_bytes(header_buf[4..8].try_into().unwrap());
                        if version == VERSION {
                            return load_bitpop(path); // new format with memmap2
                        } else if version == 3 {
                            return load_legacy_bitpop(path); // legacy format
                        }
                    }
                }
            }
        }
        Err(_) => {}
    }
    
    // Fallback: try legacy
    load_legacy_bitpop(path)
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
        let path = dir.join(format!("bitpop_v4_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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
        let path = dir.join(format!("bitpop_v4_single_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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
        let path = dir.join(format!("bitpop_bad_magic_v4_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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
        let genome = format!("{}{}{}", "ACGT".repeat(5000), "AACCGGTT", "TTTT".repeat(5000));
        bp.add_genome("large", &genome);
        bp.build();

        let dir = std::env::temp_dir();
        let path = dir.join(format!("bitpop_v4_comp_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        let path_str = path.to_str().unwrap();

        save_bitpop(&bp, path_str).unwrap();
        let file_size = std::fs::metadata(path_str).unwrap().len();

        // v5 format: stores both uncompressed BWT/SA (for memmap) + compressed FM_INDEX/GENOMES
        // 40K bases -> ~170KB SA+BWT uncompressed + ~20KB compressed sections = ~200KB total
        assert!(file_size < 300000, "Compressed file {} bytes is too large", file_size);

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
        let legacy_path = dir.join(format!("bitpop_v3_compat_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
        let legacy_path_str = legacy_path.to_str().unwrap();

        // Build v3 format data manually
        let mut header = Vec::with_capacity(64);
        header.extend_from_slice(b"BITP");
        header.extend_from_slice(&3u32.to_le_bytes()); // version 3
        header.extend_from_slice(&6u16.to_le_bytes());  // k=6
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
        
        let bwt_packed_len = (bwt_len + 3) / 4;
        let mut bwt_packed = vec![0u8; bwt_packed_len];
        for i in 0..bwt_len {
            let v = (fm.bwt_at(i) & 3) as u8;
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
            fm_data.extend_from_slice(&(gid as u32).to_le_bytes());
        }
        fm_data.extend_from_slice(&32u32.to_le_bytes());

        let sentinel_mask_len = (bwt_len + 7) / 8;
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
        let path = dir.join(format!("bitpop_v4_auto_{}_{}.bitpop", std::process::id(), std::time::SystemTime::now().elapsed().unwrap().as_nanos()));
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

/// Compute SHA256 hash of data.
fn sha256(data: &[u8]) -> [u8; 32] {
    use sha2::Digest;
    let mut hasher = sha2::Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}
