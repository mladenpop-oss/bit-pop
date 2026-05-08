/// Delta encoding with Variable-Length Integer (VLI) compression.
///
/// Instead of storing absolute positions [57, 58, 59, 60], we store [57, +1, +1, +1].
/// Small deltas are extremely common in text (consecutive tokens), so VLI encoding
/// gives us massive compression: most deltas fit in 1-2 bytes instead of 4 bytes (u32).
///
/// VLI format (MSB-first, 7 data bits + 1 continuation bit per byte):
///   Byte: [continuation][7 bits of data]
///   If continuation bit is 1, more bytes follow.
///   Example: 300 = 0b100101101 → split into 7-bit groups → 0b00010 0b0101101
///   Encoded as: 0x32 0x2D (continuation=1 on first byte)
///
/// Encode a u64 into VLI format. Returns bytes in MSB-first order.
/// Uses 7-bit groups, max 10 bytes for u64.
pub fn vli_encode_u64(value: u64) -> Vec<u8> {
    if value == 0 {
        return vec![0];
    }

    let mut result = Vec::with_capacity(10);
    let mut val = value;

    let mut groups = Vec::new();
    while val > 0 {
        groups.push((val & 0x7F) as u8);
        val >>= 7;
    }

    for (i, &group) in groups.iter().rev().enumerate() {
        if i < groups.len() - 1 {
            result.push(group | 0x80);
        } else {
            result.push(group);
        }
    }

    result
}

/// Decode a VLI-encoded u64 value from a byte slice.
/// Returns (value, bytes_consumed).
pub fn vli_decode_u64(bytes: &[u8]) -> Option<(u64, usize)> {
    if bytes.is_empty() {
        return None;
    }

    let mut result: u64 = 0;
    let mut i = 0;

    for &byte in bytes {
        if i >= 10 {
            return None;
        }
        result = (result << 7) | (byte & 0x7F) as u64;
        i += 1;
        if byte & 0x80 == 0 {
            break;
        }
    }

    Some((result, i))
}

/// Block-based delta + VLI compression for u64 positions with skip index.
pub struct BlockCompressedPositionsU64 {
    pub data: Vec<u8>,
    pub skip_offsets: Vec<usize>,
    pub block_size: usize,
    /// Total number of positions, stored at compression time for O(1) lookup.
    pub num_positions: usize,
}

impl BlockCompressedPositionsU64 {
    pub fn num_blocks(&self) -> usize {
        self.skip_offsets.len()
    }

    /// Total number of positions — O(1), no decompression needed.
    pub fn num_positions(&self) -> usize {
        self.num_positions
    }

    pub fn block_positions(&self, block_idx: usize) -> usize {
        let start = self.skip_offsets[block_idx];
        let end = if block_idx + 1 < self.skip_offsets.len() {
            self.skip_offsets[block_idx + 1]
        } else {
            self.data.len()
        };
        let block_data = &self.data[start..end];
        let mut count = 0usize;
        let mut offset = 0;
        while offset < block_data.len() {
            if let Some((_, consumed)) = vli_decode_u64(&block_data[offset..]) {
                count += 1;
                offset += consumed;
            } else {
                break;
            }
        }
        count
    }

    pub fn decompress_block(&self, block_idx: usize) -> Vec<u64> {
        if block_idx >= self.num_blocks() {
            return Vec::new();
        }
        let start = self.skip_offsets[block_idx];
        let end = if block_idx + 1 < self.skip_offsets.len() {
            self.skip_offsets[block_idx + 1]
        } else {
            self.data.len()
        };
        let block_data = &self.data[start..end];

        let mut values = Vec::new();
        let mut offset = 0;
        while offset < block_data.len() {
            if let Some((val, consumed)) = vli_decode_u64(&block_data[offset..]) {
                values.push(val);
                offset += consumed;
            } else {
                break;
            }
        }

        if values.is_empty() {
            return Vec::new();
        }

        let mut positions = Vec::with_capacity(values.len());
        positions.push(values[0]);
        let mut current = values[0];
        for &val in &values[1..] {
            current += val;
            positions.push(current);
        }
        positions
    }

    pub fn decompress_all(&self) -> Vec<u64> {
        let mut result = Vec::new();
        for block_idx in 0..self.num_blocks() {
            result.extend(self.decompress_block(block_idx));
        }
        result
    }
}

/// Compress u64 positions into independent blocks with skip index.
pub fn compress_positions_blocked_u64(
    positions: &[u64],
    block_size: usize,
) -> BlockCompressedPositionsU64 {
    if positions.is_empty() {
        return BlockCompressedPositionsU64 {
            data: Vec::new(),
            skip_offsets: Vec::new(),
            block_size,
            num_positions: 0,
        };
    }

    let total = positions.len();
    let num_blocks = total.div_ceil(block_size);
    let mut skip_offsets = Vec::with_capacity(num_blocks);
    let mut data = Vec::new();

    for block_idx in 0..num_blocks {
        skip_offsets.push(data.len());
        let start = block_idx * block_size;
        let end = (start + block_size).min(total);
        let block = &positions[start..end];

        if block.is_empty() {
            continue;
        }

        data.extend_from_slice(&vli_encode_u64(block[0]));

        for w in block.windows(2) {
            let delta = w[1] - w[0];
            data.extend_from_slice(&vli_encode_u64(delta));
        }
    }

    BlockCompressedPositionsU64 {
        data,
        skip_offsets,
        block_size,
        num_positions: total,
    }
}

/// Encode a u32 into VLI format. Returns bytes in MSB-first order.
pub fn vli_encode(value: u32) -> Vec<u8> {
    if value == 0 {
        return vec![0];
    }

    let mut result = Vec::with_capacity(5); // max 5 bytes for u32
    let mut val = value;

    // Extract 7-bit groups from LSB
    let mut groups = Vec::new();
    while val > 0 {
        groups.push((val & 0x7F) as u8);
        val >>= 7;
    }

    // Reverse to get MSB-first, set continuation bits
    for (i, &group) in groups.iter().rev().enumerate() {
        if i < groups.len() - 1 {
            result.push(group | 0x80); // continuation bit set
        } else {
            result.push(group); // last byte, no continuation
        }
    }

    result
}

/// Decode a VLI-encoded value from a byte slice.
/// Returns (value, bytes_consumed).
pub fn vli_decode(bytes: &[u8]) -> Option<(u32, usize)> {
    if bytes.is_empty() {
        return None;
    }

    let mut result: u32 = 0;
    let mut i = 0;

    for &byte in bytes {
        if i >= 5 {
            return None; // overflow protection
        }
        result = (result << 7) | (byte & 0x7F) as u32;
        i += 1;
        if byte & 0x80 == 0 {
            break; // no more continuation
        }
    }

    Some((result, i))
}

/// Delta-encode a sequence of positions.
/// Input: [10, 12, 13, 15, 20]
/// Output: [10, 2, 1, 2, 5] (first value is absolute, rest are deltas)
pub fn delta_encode(positions: &[u32]) -> Vec<u32> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::with_capacity(positions.len());
    encoded.push(positions[0]);

    let mut prev = positions[0];
    for &pos in &positions[1..] {
        let delta = pos - prev;
        encoded.push(delta);
        prev = pos;
    }

    encoded
}

/// Delta-decode a sequence back to absolute positions.
pub fn delta_decode(encoded: &[u32]) -> Vec<u32> {
    if encoded.is_empty() {
        return Vec::new();
    }

    let mut decoded = Vec::with_capacity(encoded.len());
    let mut pos = encoded[0];
    decoded.push(pos);

    for &delta in &encoded[1..] {
        pos += delta;
        decoded.push(pos);
    }

    decoded
}

/// Full compression pipeline: positions → delta → VLI bytes.
pub fn compress_positions(positions: &[u32]) -> Vec<u8> {
    if positions.is_empty() {
        return Vec::new();
    }

    let deltas = delta_encode(positions);
    let mut compressed = Vec::new();

    for &val in &deltas {
        let encoded = vli_encode(val);
        compressed.extend_from_slice(&encoded);
    }

    compressed
}

/// Full decompression pipeline: VLI bytes → delta decode → absolute positions.
pub fn decompress_positions(compressed: &[u8]) -> Vec<u32> {
    if compressed.is_empty() {
        return Vec::new();
    }

    let mut deltas = Vec::new();
    let mut i = 0;

    while i < compressed.len() {
        if let Some((val, consumed)) = vli_decode(&compressed[i..]) {
            deltas.push(val);
            i += consumed;
        } else {
            break;
        }
    }

    delta_decode(&deltas)
}

/// Iterator over VLI-encoded values without full decompression.
/// This allows early-exit during search without materializing all positions.
pub struct VliIterator<'a> {
    data: &'a [u8],
    offset: usize,
}

impl<'a> VliIterator<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, offset: 0 }
    }
}

impl<'a> Iterator for VliIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        if self.offset >= self.data.len() {
            return None;
        }

        let (val, consumed) = vli_decode(&self.data[self.offset..])?;
        self.offset += consumed;
        Some(val)
    }
}

/// Delta-aware iterator that yields absolute positions from VLI-encoded delta stream.
pub struct DeltaVliIterator<'a> {
    vli: VliIterator<'a>,
    current_pos: Option<u32>,
    first: bool,
}

impl<'a> DeltaVliIterator<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            vli: VliIterator::new(data),
            current_pos: None,
            first: true,
        }
    }
}

impl<'a> Iterator for DeltaVliIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        let raw = self.vli.next()?;

        if self.first {
            self.current_pos = Some(raw);
            self.first = false;
        } else {
            self.current_pos = self.current_pos.map(|p| p + raw);
        }

        self.current_pos
    }
}

/// Block-based delta + VLI compression with skip index for random access.
///
/// Instead of one giant VLI stream that must be decoded sequentially,
/// positions are divided into independent blocks. Each block stores:
///   [absolute_first_position: VLI][delta_1: VLI][delta_2: VLI]...
///
/// A skip index maps block_index → byte_offset, allowing O(1) jumps
/// to any block without decoding preceding blocks.
///
/// Example (block_size=4):
///   Positions: [0, 1, 2, 5, 100, 101, 103, 200]
///   Block 0: [0, +1, +1, +3]        → absolute positions [0, 1, 2, 5]
///   Block 1: [100, +1, +2]          → absolute positions [100, 101, 103]
///   Block 2: [200]                  → absolute positions [200]
///   Skip offsets: [0, byte_offset_of_block1, byte_offset_of_block2]
///
/// Benefits:
/// - Random access: decompress any block without touching others
/// - Anchor search jumps to relevant block instead of scanning all deltas
/// - Still gets delta + VLI compression benefits within each block
pub struct BlockCompressedPositions {
    /// Concatenated compressed blocks
    pub data: Vec<u8>,
    /// Byte offset of each block in `data`. skip_offsets[i] = byte start of block i.
    pub skip_offsets: Vec<usize>,
    /// Number of positions per block (last block may have fewer).
    pub block_size: usize,
}

impl BlockCompressedPositions {
    /// Number of blocks.
    pub fn num_blocks(&self) -> usize {
        self.skip_offsets.len()
    }

    /// Total number of positions.
    pub fn num_positions(&self) -> usize {
        if self.num_blocks() == 0 {
            return 0;
        }
        let full_blocks = self.num_blocks() - 1;
        full_blocks * self.block_size + self.block_positions(self.num_blocks() - 1)
    }

    /// How many positions are in a specific block.
    pub fn block_positions(&self, block_idx: usize) -> usize {
        let start = self.skip_offsets[block_idx];
        let end = if block_idx + 1 < self.skip_offsets.len() {
            self.skip_offsets[block_idx + 1]
        } else {
            self.data.len()
        };
        let block_data = &self.data[start..end];
        let mut count = 0u32;
        let mut offset = 0;
        while offset < block_data.len() {
            if let Some((_, consumed)) = vli_decode(&block_data[offset..]) {
                count += 1;
                offset += consumed;
            } else {
                break;
            }
        }
        count as usize
    }

    /// Decompress a single block independently.
    pub fn decompress_block(&self, block_idx: usize) -> Vec<u32> {
        if block_idx >= self.num_blocks() {
            return Vec::new();
        }
        let start = self.skip_offsets[block_idx];
        let end = if block_idx + 1 < self.skip_offsets.len() {
            self.skip_offsets[block_idx + 1]
        } else {
            self.data.len()
        };
        let block_data = &self.data[start..end];

        let mut values = Vec::new();
        let mut offset = 0;
        while offset < block_data.len() {
            if let Some((val, consumed)) = vli_decode(&block_data[offset..]) {
                values.push(val);
                offset += consumed;
            } else {
                break;
            }
        }

        let mut positions = Vec::with_capacity(values.len());
        if values.is_empty() {
            return positions;
        }

        positions.push(values[0]);
        let mut current = values[0];
        for &val in &values[1..] {
            current += val;
            positions.push(current);
        }
        positions
    }

    /// Check if a target absolute position exists in this block without
    /// materializing all positions. Walks the VLI stream and exits early
    /// once the target is found or surpassed.
    pub fn contains_position(&self, block_idx: usize, target: u32) -> bool {
        if block_idx >= self.num_blocks() {
            return false;
        }
        let start = self.skip_offsets[block_idx];
        let end = if block_idx + 1 < self.skip_offsets.len() {
            self.skip_offsets[block_idx + 1]
        } else {
            self.data.len()
        };
        let block_data = &self.data[start..end];
        if block_data.is_empty() {
            return false;
        }

        let mut offset = 0;
        let first_val = match vli_decode(&block_data[offset..]) {
            Some((v, c)) => {
                offset = c;
                v
            }
            None => return false,
        };

        if target == first_val {
            return true;
        }
        if target < first_val {
            return false;
        }

        let mut current = first_val;
        while offset < block_data.len() {
            let (delta, consumed) = match vli_decode(&block_data[offset..]) {
                Some(d) => d,
                None => break,
            };
            offset += consumed;
            current += delta;
            if current == target {
                return true;
            }
            if current > target {
                return false;
            }
        }
        false
    }

    /// Find the block that could contain the target position, or None.
    /// Uses binary search on first-position-per-block for O(log B) lookup.
    pub fn find_block_for(&self, target: u32) -> Option<usize> {
        if self.num_blocks() == 0 {
            return None;
        }

        // Binary search: find the largest block index where first_val <= target
        let mut best: Option<usize> = None;
        let mut lo = 0usize;
        let mut hi = self.num_blocks();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let start = self.skip_offsets[mid];
            let block_data = &self.data[start..];
            if block_data.is_empty() {
                hi = mid;
                continue;
            }

            let first_val = match vli_decode(block_data) {
                Some((v, _)) => v,
                None => {
                    hi = mid;
                    continue;
                }
            };

            if first_val <= target {
                best = Some(mid);
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }

        let candidate = best?;

        if self.contains_position(candidate, target) {
            return Some(candidate);
        }

        None
    }

    /// Decompress all positions.
    pub fn decompress_all(&self) -> Vec<u32> {
        let mut result = Vec::new();
        for block_idx in 0..self.num_blocks() {
            result.extend(self.decompress_block(block_idx));
        }
        result
    }

    /// Binary search for a target absolute position across blocks.
    /// Returns (block_idx, position_in_block) if found, or None.
    pub fn find_position(&self, target: u32) -> Option<usize> {
        for block_idx in 0..self.num_blocks() {
            let block = self.decompress_block(block_idx);
            if block.is_empty() {
                continue;
            }
            if target < block[0] || target > *block.last().unwrap() {
                continue;
            }
            if let Ok(pos) = block.binary_search(&target) {
                return Some(pos);
            }
        }
        None
    }
}

/// Compress positions into independent blocks with skip index.
pub fn compress_positions_blocked(
    positions: &[u32],
    block_size: usize,
) -> BlockCompressedPositions {
    if positions.is_empty() {
        return BlockCompressedPositions {
            data: Vec::new(),
            skip_offsets: Vec::new(),
            block_size,
        };
    }

    let num_blocks = positions.len().div_ceil(block_size);
    let mut skip_offsets = Vec::with_capacity(num_blocks);
    let mut data = Vec::new();

    for block_idx in 0..num_blocks {
        skip_offsets.push(data.len());
        let start = block_idx * block_size;
        let end = (start + block_size).min(positions.len());
        let block = &positions[start..end];

        if block.is_empty() {
            continue;
        }

        // First value is absolute position
        data.extend_from_slice(&vli_encode(block[0]));

        // Rest are deltas within the block
        for w in block.windows(2) {
            let delta = w[1] - w[0];
            data.extend_from_slice(&vli_encode(delta));
        }
    }

    BlockCompressedPositions {
        data,
        skip_offsets,
        block_size,
    }
}

/// u64 iterator that walks through all positions block by block,
/// yielding absolute positions without full decompression.
pub struct BlockedDeltaIteratorU64<'a> {
    compressed: &'a BlockCompressedPositionsU64,
    current_block: usize,
    block_data_offset: usize,
    current_absolute: u64,
    first_in_block: bool,
}

impl<'a> BlockedDeltaIteratorU64<'a> {
    pub fn new(compressed: &'a BlockCompressedPositionsU64) -> Self {
        Self {
            compressed,
            current_block: 0,
            block_data_offset: 0,
            current_absolute: 0,
            first_in_block: true,
        }
    }
}

impl<'a> Iterator for BlockedDeltaIteratorU64<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        let mut block_end = if self.current_block + 1 < self.compressed.skip_offsets.len() {
            self.compressed.skip_offsets[self.current_block + 1]
        } else {
            self.compressed.data.len()
        };

        while self.block_data_offset >= block_end {
            self.current_block += 1;
            if self.current_block >= self.compressed.num_blocks() {
                return None;
            }
            self.block_data_offset = self.compressed.skip_offsets[self.current_block];
            self.first_in_block = true;
            block_end = if self.current_block + 1 < self.compressed.skip_offsets.len() {
                self.compressed.skip_offsets[self.current_block + 1]
            } else {
                self.compressed.data.len()
            };
        }

        let remaining = &self.compressed.data[self.block_data_offset..block_end];
        let (val, consumed) = vli_decode_u64(remaining)?;
        self.block_data_offset += consumed;

        if self.first_in_block {
            self.current_absolute = val;
            self.first_in_block = false;
        } else {
            self.current_absolute += val;
        }

        Some(self.current_absolute)
    }
}

/// Iterator that walks through all positions block by block,
/// yielding absolute positions without full decompression.
pub struct BlockedDeltaIterator<'a> {
    compressed: &'a BlockCompressedPositions,
    current_block: usize,
    block_data_offset: usize,
    current_absolute: u32,
    block_data_start: usize,
    first_in_block: bool,
}

impl<'a> BlockedDeltaIterator<'a> {
    pub fn new(compressed: &'a BlockCompressedPositions) -> Self {
        Self {
            compressed,
            current_block: 0,
            block_data_offset: 0,
            current_absolute: 0,
            block_data_start: 0,
            first_in_block: true,
        }
    }

    fn load_block(&mut self) -> bool {
        if self.current_block >= self.compressed.num_blocks() {
            return false;
        }
        self.block_data_start = self.compressed.skip_offsets[self.current_block];
        self.block_data_offset = self.block_data_start;
        self.first_in_block = true;
        true
    }
}

impl<'a> Iterator for BlockedDeltaIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        let mut block_end = if self.current_block + 1 < self.compressed.skip_offsets.len() {
            self.compressed.skip_offsets[self.current_block + 1]
        } else {
            self.compressed.data.len()
        };

        while self.block_data_offset >= block_end {
            self.current_block += 1;
            if !self.load_block() {
                return None;
            }
            block_end = if self.current_block + 1 < self.compressed.skip_offsets.len() {
                self.compressed.skip_offsets[self.current_block + 1]
            } else {
                self.compressed.data.len()
            };
        }

        let remaining = &self.compressed.data[self.block_data_offset..block_end];
        let (val, consumed) = vli_decode(remaining)?;
        self.block_data_offset += consumed;

        if self.first_in_block {
            self.current_absolute = val;
            self.first_in_block = false;
        } else {
            self.current_absolute += val;
        }

        Some(self.current_absolute)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vli_encode_decode_zero() {
        let encoded = vli_encode(0);
        assert_eq!(encoded, vec![0]);
        assert_eq!(vli_decode(&encoded), Some((0, 1)));
    }

    #[test]
    fn test_vli_encode_decode_small() {
        for i in 1..=127u32 {
            let encoded = vli_encode(i);
            assert_eq!(vli_decode(&encoded), Some((i, 1)));
        }
    }

    #[test]
    fn test_vli_encode_decode_medium() {
        let val: u32 = 300;
        let encoded = vli_encode(val);
        assert_eq!(encoded.len(), 2);
        assert_eq!(vli_decode(&encoded), Some((val, 2)));
    }

    #[test]
    fn test_vli_encode_decode_large() {
        let val: u32 = 0xFFFFFFFF;
        let encoded = vli_encode(val);
        assert_eq!(encoded.len(), 5);
        assert_eq!(vli_decode(&encoded), Some((val, 5)));
    }

    #[test]
    fn test_delta_encode_decode() {
        let positions = vec![10, 12, 13, 15, 20];
        let encoded = delta_encode(&positions);
        assert_eq!(encoded, vec![10, 2, 1, 2, 5]);

        let decoded = delta_decode(&encoded);
        assert_eq!(decoded, positions);
    }

    #[test]
    fn test_compress_decompress() {
        let positions = vec![0, 1, 2, 3, 100, 101, 102, 5000];
        let compressed = compress_positions(&positions);
        let decompressed = decompress_positions(&compressed);
        assert_eq!(decompressed, positions);

        // Compression ratio check
        let original_bytes = positions.len() * 4; // u32 per position
        assert!(
            compressed.len() < original_bytes,
            "Compressed {} bytes should be less than original {} bytes",
            compressed.len(),
            original_bytes
        );
    }

    #[test]
    fn test_consecutive_positions_compress_well() {
        // Simulate 1000 consecutive positions (common in text)
        let positions: Vec<u32> = (0..1000).collect();
        let compressed = compress_positions(&positions);
        let original_bytes = positions.len() * 4;

        // First value takes ~2 bytes, rest are +1 (1 byte each)
        let ratio = compressed.len() as f64 / original_bytes as f64;
        assert!(
            ratio < 0.3,
            "Consecutive positions should compress to <30%%, got {:.1}%%",
            ratio * 100.0
        );
    }

    #[test]
    fn test_delta_vli_iterator() {
        let positions = vec![10, 12, 13, 15];
        let compressed = compress_positions(&positions);
        let iter = DeltaVliIterator::new(&compressed);
        let collected: Vec<u32> = iter.collect();
        assert_eq!(collected, positions);
    }

    #[test]
    fn test_empty() {
        assert!(compress_positions(&[]).is_empty());
        assert!(decompress_positions(&[]).is_empty());
    }

    #[test]
    fn test_block_compress_decompress_all() {
        let positions: Vec<u32> = (0..100).collect();
        let blocked = compress_positions_blocked(&positions, 10);
        let decompressed = blocked.decompress_all();
        assert_eq!(decompressed, positions);
    }

    #[test]
    fn test_block_independent_decompression() {
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 4);
        assert_eq!(blocked.num_blocks(), 2);

        let block0 = blocked.decompress_block(0);
        assert_eq!(block0, vec![0, 1, 2, 5]);

        let block1 = blocked.decompress_block(1);
        assert_eq!(block1, vec![100, 101, 103, 200]);
    }

    #[test]
    fn test_block_iterator() {
        let positions: Vec<u32> = (0..50).collect();
        let blocked = compress_positions_blocked(&positions, 7);
        let iter = BlockedDeltaIterator::new(&blocked);
        let collected: Vec<u32> = iter.collect();
        assert_eq!(collected, positions);
    }

    #[test]
    fn test_block_find_position() {
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 4);

        assert_eq!(blocked.find_position(0), Some(0));
        assert_eq!(blocked.find_position(5), Some(3));
        assert_eq!(blocked.find_position(100), Some(0));
        assert_eq!(blocked.find_position(200), Some(3));
        assert_eq!(blocked.find_position(999), None);
    }

    #[test]
    fn test_block_compression_savings() {
        let positions: Vec<u32> = (0..1000).collect();
        let sequential = compress_positions(&positions);
        let blocked = compress_positions_blocked(&positions, 64);

        // Block compression has slight overhead from skip_offsets and repeated absolute positions,
        // but should be within reasonable range of sequential compression
        let overhead_ratio = blocked.data.len() as f64 / sequential.len() as f64;
        assert!(
            overhead_ratio < 2.0,
            "Block compression {} should not be more than 2x sequential {}",
            blocked.data.len(),
            sequential.len()
        );
    }

    #[test]
    fn test_block_random_positions() {
        let positions = vec![10, 55, 120, 200, 350, 410, 500, 750, 900, 1024];
        let blocked = compress_positions_blocked(&positions, 3);
        assert_eq!(blocked.num_blocks(), 4);
        assert_eq!(blocked.decompress_all(), positions);
    }

    #[test]
    fn test_block_empty() {
        let blocked = compress_positions_blocked(&[], 10);
        assert_eq!(blocked.num_blocks(), 0);
        assert!(blocked.decompress_all().is_empty());
    }

    #[test]
    fn test_block_single_element() {
        let blocked = compress_positions_blocked(&[42], 10);
        assert_eq!(blocked.num_blocks(), 1);
        assert_eq!(blocked.decompress_all(), vec![42]);
    }

    #[test]
    fn test_block_num_positions() {
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 3);
        assert_eq!(blocked.num_blocks(), 3);
        assert_eq!(blocked.num_positions(), 8);
    }

    #[test]
    fn test_contains_position_found() {
        // positions [0,1,2,5,100,101,103,200] with block_size=4 → 2 blocks:
        // Block 0: [0, 1, 2, 5], Block 1: [100, 101, 103, 200]
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 4);
        assert_eq!(blocked.num_blocks(), 2);

        assert!(blocked.contains_position(0, 0));
        assert!(blocked.contains_position(0, 1));
        assert!(blocked.contains_position(0, 5));
        assert!(!blocked.contains_position(0, 3));
        assert!(blocked.contains_position(1, 100));
        assert!(blocked.contains_position(1, 103));
        assert!(blocked.contains_position(1, 200));
    }

    #[test]
    fn test_contains_position_out_of_range() {
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 4);

        assert!(!blocked.contains_position(0, 100));
        assert!(!blocked.contains_position(1, 5));
    }

    #[test]
    fn test_find_block_for() {
        let positions = vec![0, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked(&positions, 4);

        assert_eq!(blocked.find_block_for(0), Some(0));
        assert_eq!(blocked.find_block_for(5), Some(0));
        assert_eq!(blocked.find_block_for(100), Some(1));
        assert_eq!(blocked.find_block_for(103), Some(1));
        assert_eq!(blocked.find_block_for(200), Some(1));
        assert_eq!(blocked.find_block_for(999), None);
        assert_eq!(blocked.find_block_for(50), None);
    }

    #[test]
    fn test_contains_vs_decompress_consistency() {
        let positions: Vec<u32> = (0..1000).step_by(7).collect();
        let blocked = compress_positions_blocked(&positions, 16);

        for &pos in &positions {
            let candidate = blocked.find_block_for(pos);
            assert!(candidate.is_some(), "Should find block for {}", pos);
            if let Some(bi) = candidate {
                assert!(blocked.contains_position(bi, pos), "Should contain {}", pos);
            }
        }

        // Check some non-existent positions
        assert!(!blocked.contains_position(0, 1));
        assert!(!blocked.contains_position(0, 15));
    }

    #[test]
    fn test_vli_encode_u64_zero() {
        let encoded = vli_encode_u64(0);
        assert_eq!(encoded, vec![0]);
        assert_eq!(vli_decode_u64(&encoded), Some((0, 1)));
    }

    #[test]
    fn test_vli_encode_u64_small() {
        for i in 1..=127u64 {
            let encoded = vli_encode_u64(i);
            assert_eq!(vli_decode_u64(&encoded), Some((i, 1)));
        }
    }

    #[test]
    fn test_vli_encode_u64_medium() {
        let val: u64 = 300;
        let encoded = vli_encode_u64(val);
        assert_eq!(encoded.len(), 2);
        assert_eq!(vli_decode_u64(&encoded), Some((val, 2)));
    }

    #[test]
    fn test_vli_encode_u64_large() {
        let val: u64 = 0xFFFFFFFFFFFFFFFF;
        let encoded = vli_encode_u64(val);
        assert_eq!(encoded.len(), 10);
        assert_eq!(vli_decode_u64(&encoded), Some((val, 10)));
    }

    #[test]
    fn test_vli_encode_u64_genome_scale() {
        // Human genome scale: 3 billion
        let val: u64 = 3_000_000_000;
        let encoded = vli_encode_u64(val);
        assert_eq!(vli_decode_u64(&encoded), Some((val, 5)));
    }

    #[test]
    fn test_block_compress_u64_decompress_all() {
        let positions: Vec<u64> = (0..100).map(|x| x * 1000).collect();
        let blocked = compress_positions_blocked_u64(&positions, 10);
        let decompressed = blocked.decompress_all();
        assert_eq!(decompressed, positions);
    }

    #[test]
    fn test_block_u64_independent_decompression() {
        let positions = vec![0u64, 1, 2, 5, 100, 101, 103, 200];
        let blocked = compress_positions_blocked_u64(&positions, 4);
        assert_eq!(blocked.num_blocks(), 2);

        let block0 = blocked.decompress_block(0);
        assert_eq!(block0, vec![0u64, 1, 2, 5]);

        let block1 = blocked.decompress_block(1);
        assert_eq!(block1, vec![100u64, 101, 103, 200]);
    }

    #[test]
    fn test_block_u64_compression_savings() {
        let positions: Vec<u64> = (0..1000).collect();
        let blocked = compress_positions_blocked_u64(&positions, 64);

        let raw_bytes = positions.len() * std::mem::size_of::<u64>();
        let ratio = blocked.data.len() as f64 / raw_bytes as f64;
        assert!(
            ratio < 0.15,
            "Sequential u64 positions should compress to <15%%, got {:.1}%%",
            ratio * 100.0
        );
    }

    #[test]
    fn test_block_u64_empty() {
        let blocked = compress_positions_blocked_u64(&[], 10);
        assert_eq!(blocked.num_blocks(), 0);
        assert!(blocked.decompress_all().is_empty());
    }

    #[test]
    fn test_block_u64_single_element() {
        let blocked = compress_positions_blocked_u64(&[42u64], 10);
        assert_eq!(blocked.num_blocks(), 1);
        assert_eq!(blocked.decompress_all(), vec![42u64]);
    }

    #[test]
    fn test_block_u64_large_positions() {
        let positions: Vec<u64> = (0..100).map(|i| i as u64 * 10_000_000).collect();
        let blocked = compress_positions_blocked_u64(&positions, 16);
        assert_eq!(blocked.decompress_all(), positions);
    }

    #[test]
    fn test_block_u64_iterator() {
        let positions: Vec<u64> = (0..50).map(|x| x * 100).collect();
        let blocked = compress_positions_blocked_u64(&positions, 7);
        let iter = BlockedDeltaIteratorU64::new(&blocked);
        let collected: Vec<u64> = iter.collect();
        assert_eq!(collected, positions);
    }

    #[test]
    fn test_block_u64_num_positions_o1() {
        let positions: Vec<u64> = (0..1234).map(|x| x * 1000).collect();
        let blocked = compress_positions_blocked_u64(&positions, 256);
        assert_eq!(blocked.num_positions(), 1234);
        assert_eq!(blocked.num_positions(), blocked.decompress_all().len());
    }

    #[test]
    fn test_block_u64_num_positions_zero() {
        let blocked = compress_positions_blocked_u64(&[], 10);
        assert_eq!(blocked.num_positions(), 0);
    }
}
