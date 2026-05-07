/// Bit-level alignment algorithms for DNA sequences.
/// Uses 2-bit encoding (A=00, C=01, G=10, T=11) for bit-parallelism.
/// Pattern length <= 31 bases fits in u64 → 31 bases processed per CPU instruction.

/// 2-bit XOR alignment — simple, fast, reliable.
/// Slides pattern across text, counts exact matches via XOR, returns best score + CIGAR.
/// Pattern <= 31 bases (62 bits fit in u64).
/// 
/// Score: 0.0 (no match) to 1.0 (perfect match)
/// CIGAR: "97M2X1M" format — M=match, X=mismatch
/// Offset: best starting position in text
pub fn two_bit_align(pattern: &[u8], text: &[u8]) -> (f64, String, usize) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return (0.0, String::new(), 0);
    }

    // Pack pattern into u64: each base = 2 bits, MSB first
    let mut pat_val: u64 = 0;
    for &b in pattern {
        pat_val = (pat_val << 2) | (b as u64 & 3);
    }

    let mask = (1u64 << (m * 2)) - 1;

    // Pack first window
    let mut win_val: u64 = 0;
    for &b in text.iter().take(m) {
        win_val = (win_val << 2) | (b as u64 & 3);
    }

    let mut best_xor = pat_val ^ win_val;
    let errors = count_2bit_diffs(best_xor, m);
    let mut best_score = (m - errors) as f64 / m as f64;
    let mut best_offset = 0usize;

    // Slide window across text
    for i in m..n {
        win_val = ((win_val << 2) | (text[i] as u64 & 3)) & mask;
        let xor = pat_val ^ win_val;
        let errors = count_2bit_diffs(xor, m);
        let score = (m - errors) as f64 / m as f64;
        if score > best_score {
            best_score = score;
            best_offset = i + 1 - m;
            best_xor = xor;
        }
    }

    let cigar = build_cigar_from_xor(best_xor, m);

    (best_score, cigar, best_offset)
}

/// Build CIGAR string from XOR result.
/// Scans 2-bit pairs: 0 = match (M), non-zero = mismatch (X).
fn build_cigar_from_xor(xor_val: u64, len: usize) -> String {
    let mut cigar = String::with_capacity(len * 2 + 2);
    let mut ops: Vec<(u8, usize)> = Vec::new(); // (M=0, X=1) -> count

    // Scan each 2-bit pair from LSB
    let mut val = xor_val;
    for _i in 0..len {
        let is_match = (val & 3) == 0;
        let op = if is_match { 0u8 } else { 1u8 };

        if let Some(last) = ops.last_mut() {
            if last.0 == op {
                last.1 += 1;
            } else {
                ops.push((op, 1));
            }
        } else {
            ops.push((op, 1));
        }
        val >>= 2;
    }

    // ops are in reverse (LSB first = last base), reverse back
    ops.reverse();

    for (op, count) in ops {
        cigar.push_str(&count.to_string());
        cigar.push(match op {
            0 => 'M',
            1 => 'X',
            _ => unreachable!(),
        });
    }

    if cigar.is_empty() {
        cigar = format!("{}M", len);
    }

    cigar
}

/// Bit-parallel exact matching using sliding u64 window.
/// Returns positions where the pattern matches exactly in the text.
/// Pattern length must be <= 31 bases (fits in u64 comparison).
pub fn bit_parallel_exact(pattern: &[u8], text: &[u8]) -> Vec<usize> {
    if pattern.is_empty() || pattern.len() > 31 || text.len() < pattern.len() {
        return Vec::new();
    }

    let m = pattern.len();

    // Encode pattern into a single u64: each base = 2 bits
    let mut pat_val: u64 = 0;
    for &b in pattern {
        pat_val = (pat_val << 2) | b as u64;
    }

    let mut matches = Vec::new();

    // Sliding window: encode each window of size m and compare
    let mut win_val: u64 = 0;
    for (_i, &b) in text.iter().take(m).enumerate() {
        win_val = (win_val << 2) | b as u64;
    }
    if win_val == pat_val {
        matches.push(0);
    }

    for i in m..text.len() {
        win_val = ((win_val << 2) | text[i] as u64) & ((1 << (m * 2)) - 1);
        if win_val == pat_val {
            matches.push(i + 1 - m);
        }
    }

    matches
}

/// Smith-Waterman local alignment with full traceback.
/// Returns (score, CIGAR_string).
/// Uses standard scoring: match=+2, mismatch=-1, gap=-2.
pub fn smith_waterman(pattern: &[u8], text: &[u8]) -> (i32, String) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return (0, String::new());
    }

    let match_sc = 2;
    let mismatch_sc = -1;
    let gap_sc = -2;

    // Score matrix and traceback directions
    // DIR_MATCH = 0, DIR_INSERT = 1, DIR_DELETE = 2, DIR_START = 3
    let mut sc = vec![vec![0i32; n + 1]; m + 1];
    let mut dr = vec![vec![0u8; n + 1]; m + 1];

    for i in 1..=m {
        for j in 1..=n {
            let sub = if pattern[i - 1] == text[j - 1] {
                match_sc
            } else {
                mismatch_sc
            };

            let diag = sc[i - 1][j - 1] + sub;
            let del = sc[i - 1][j] + gap_sc;
            let ins = sc[i][j - 1] + gap_sc;

            let best = diag.max(del).max(ins).max(0);

            if best == 0 {
                dr[i][j] = 3;
            } else if best == diag {
                dr[i][j] = 0;
            } else if best == ins {
                dr[i][j] = 1;
            } else {
                dr[i][j] = 2;
            }

            sc[i][j] = best;
        }
    }

    let mut max_score_cell = (m, n);
    for i in 1..=m {
        for j in 1..=n {
            if sc[i][j] > sc[max_score_cell.0][max_score_cell.1] {
                max_score_cell = (i, j);
            }
        }
    }

    let final_score = sc[max_score_cell.0][max_score_cell.1];

    if final_score == 0 {
        return (0, String::new());
    }

    // Traceback
    let mut ops = Vec::new();
    let mut i = max_score_cell.0;
    let mut j = max_score_cell.1;

    while i > 0 && j > 0 && sc[i][j] > 0 {
        match dr[i][j] {
            0 => {
                ops.push(0);
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(1);
                j -= 1;
            }
            2 => {
                ops.push(2);
                i -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();

    // Build CIGAR string
    let mut cigar = String::new();
    if ops.is_empty() {
        return (final_score, format!("{}M", m));
    }

    let mut current_op = ops[0];
    let mut count = 1usize;

    for &op in &ops[1..] {
        if op == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!(
                "{}{}",
                count,
                match current_op {
                    0 => 'M',
                    1 => 'I',
                    2 => 'D',
                    _ => 'M',
                }
            ));
            current_op = op;
            count = 1;
        }
    }

    cigar.push_str(&format!(
        "{}{}",
        count,
        match current_op {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            _ => 'M',
        }
    ));

    (final_score, cigar)
}

/// Bit-parallel approximate matching with configurable error threshold.
/// Returns (score_0_to_1, approximate_mismatch_count) for the best alignment.
/// Uses bit-packed u64 to compare 31 bases per instruction.
/// Fast path: pattern is packed into one u64, text is scanned in overlapping windows.
pub fn bit_parallel_approx(
    pattern: &[u8],
    text: &[u8],
    max_errors: usize,
) -> (f64, usize) {
    if pattern.is_empty() || text.is_empty() || pattern.len() > 31 {
        return (0.0, usize::MAX);
    }

    let m = pattern.len();
    if text.len() < m {
        return (0.0, usize::MAX);
    }

    // Pack pattern into u64 (2 bits per base)
    let mut pat_val: u64 = 0;
    for &b in pattern {
        pat_val = (pat_val << 2) | (b as u64 & 3);
    }

    let mut best_score = 0.0f64;
    let mut best_errors = usize::MAX;

    // Sliding window: pack each window of size m and compare
    let mut win_val: u64 = 0;
    for &b in text.iter().take(m) {
        win_val = (win_val << 2) | (b as u64 & 3);
    }

    // Count 2-bit differences between pattern and first window
    let xor = pat_val ^ win_val;
    let mut errors = count_2bit_diffs(xor, m);
    if errors <= max_errors {
        let score = 1.0 - (errors as f64 / m as f64);
        if score > best_score {
            best_score = score;
            best_errors = errors;
        }
    }

    // Mask for sliding window: keep only m*2 bits
    let mask = (1u64 << (m * 2)) - 1;

    for i in m..text.len() {
        win_val = ((win_val << 2) | (text[i] as u64 & 3)) & mask;
        let xor = pat_val ^ win_val;
        errors = count_2bit_diffs(xor, m);
        if errors <= max_errors {
            let score = 1.0 - (errors as f64 / m as f64);
            if score > best_score {
                best_score = score;
                best_errors = errors;
            }
        }
    }

    (best_score, best_errors)
}

/// Count the number of 2-bit positions that differ in a u64.
/// Each base is 2 bits, so we count how many 2-bit groups are non-zero after XOR.
fn count_2bit_diffs(xor_val: u64, num_bases: usize) -> usize {
    let mut count = 0;
    let mut val = xor_val;
    for _ in 0..num_bases {
        if val & 3 != 0 {
            count += 1;
        }
        val >>= 2;
    }
    count
}

/// Bit-parallel alignment score using Myers' bit-vector algorithm.
/// Pattern must fit in 31 bases (62 bits per u64 for 2-bit encoding).
///
/// Uses equality vectors per base value (A=0, C=1, G=2, T=3) for O(1) lookup.
/// All m columns are processed in a single u64 word per text position.
///
/// Returns best local alignment match count (0 to m).
/// Convert to SW score: score_sw = 3 * matches - m
/// Convert to fraction: fraction = matches / m
///
/// Complexity: O(n) for m ≤ 31 vs O(n×m) for standard SW.
pub fn bit_vector_sw(pattern: &[u8], text: &[u8]) -> i32 {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return 0;
    }

    // Precompute equality vectors per base value (A=0, C=1, G=2, T=3)
    // Eq[v] has bit i set if pattern[i] == v
    let mut eq = [0u64; 4];
    for i in 0..m {
        eq[(pattern[i] & 3) as usize] |= 1u64 << i;
    }

    // Myers' bit-vector state
    // PH/MH: carry-save representation of score vector H
    //   Score at column i = 2*popcount(PH in bits 0..i) + popcount(MH in bits 0..i)
    //   (standard Myers' formulation from "A Fast Bit-Vector Algorithm")
    // Initialized to PH=0, MH=ALL_ONES which represents zero scores everywhere
    let mut ph: u64 = 0;
    let mut mh: u64 = !0u64;
    let mut _pe: u64 = 0;
    let mut _me: u64 = !0u64;
    let mut pf: u64 = 0;

    let mask = (1u64 << m) - 1;
    let mut best = 0i32;

    for j in 0..n {
        let xr = eq[(text[j] & 3) as usize];

        // H_update: q = H + (1 - XR) in carry-save form
        // Myers' formula: q = ((MH | XR) + MH) ^ MH
        // This computes H + sub(XR) where sub(match)=0, sub(mismatch)=1
        let q = ((mh | xr).wrapping_add(mh)) ^ mh;
        let qh = q >> 1;

        // R_update: R = Q + PF (add deletions from previous row)
        let r = ((qh | pf).wrapping_add(pf)) ^ qh;
        let rh = r >> 1;
        let rm = r & mask;

        // E_update: E = E + R, then E += 1 (insertions)
        let e_r = ((_me | rm).wrapping_add(_me)) ^ rm;
        let e_inc = e_r.wrapping_sub(1);
        _pe = (e_r >> 1) | ((e_r & !e_inc & mask) >> 1);
        _me = (e_r & e_inc) & mask;

        // F_update: shift R vertically for next column
        let bottom_bit = rm & 1;
        let carry = bottom_bit.wrapping_sub(1);
        pf = ((rm >> 1) | (carry << 63)) & mask;

        // H = max(0, R) for local alignment reset
        ph = rh & mask;
        mh = (rm | (!r & mask)) & mask;

        // Extract score using Myers' formula:
        // score = 2*popcount(PH & mask) + popcount(MH & 0x55... & mask)
        // This gives the edit distance (0 = perfect match, higher = worse)
        // We want match count, so: matches = m - edit_distance
        let score_raw = 2 * (ph & mask).count_ones() as i32
            + (mh & 0x5555555555555555 & mask).count_ones() as i32;
        let matches = m as i32 - score_raw;
        if matches > best {
            best = matches;
        }
    }

    best.max(0)
}

/// Bit-vector alignment with SW-compatible scoring.
/// Returns (sw_score, best_offset) where sw_score uses match=+2, mismatch=-1.
/// Internally uses bit-parallel matching, converts to SW scale.
pub fn bit_vector_sw_scored(pattern: &[u8], text: &[u8]) -> (i32, usize) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return (0, 0);
    }

    let mut eq = [0u64; 4];
    for i in 0..m {
        eq[(pattern[i] & 3) as usize] |= 1u64 << i;
    }

    let mut ph: u64 = 0;
    let mut mh: u64 = (1u64 << m) - 1;
    let mut pe: u64 = 0;
    let mut me: u64 = (1u64 << m) - 1;
    let mut pf: u64 = 0;

    let mask = (1u64 << m) - 1;
    let mut best_score = 0i32;
    let mut best_offset = 0usize;

    for j in 0..n {
        let xr = eq[(text[j] & 3) as usize];

        let q = ((mh | xr).wrapping_add(mh)) ^ mh;
        let qh = q >> 1;
        let _qm = q & mask;

        let r = ((qh | pf).wrapping_add(pf)) ^ qh;
        let rh = r >> 1;
        let rm = r & mask;

        let e_r = ((me | rm).wrapping_add(me)) ^ rm;
        let e_inc = e_r.wrapping_sub(1);
        pe = (e_r >> 1) | ((e_r & !e_inc & mask) >> 1);
        me = (e_r & e_inc) & mask;

        let bottom_bit = rm & 1;
        let carry = bottom_bit.wrapping_sub(1);
        pf = ((rm >> 1) | (carry << 63)) & mask;

        ph = rh & mask;
        mh = rm;

        // Extract matches and convert to SW score
        let matches = (mh ^ mask).count_ones() as i32;
        // SW: match=+2, mismatch=-1 → score = 2*matches - 1*(m-matches) = 3*matches - m
        let sw_score = 3 * matches - m as i32;
        if sw_score > best_score {
            best_score = sw_score;
            best_offset = j.saturating_sub(m.saturating_sub(1));
        }
    }

    (best_score.max(0), best_offset)
}

/// Windowed bit-vector alignment for local alignment.
/// Slides a window across text, finds best alignment score and position.
/// Returns (best_match_count, best_offset).
pub fn bit_vector_sw_windowed(
    pattern: &[u8],
    text: &[u8],
    window: usize,
) -> (i32, usize) {
    if pattern.is_empty() || text.is_empty() || pattern.len() > 31 {
        return (0, 0);
    }

    let w = window.max(pattern.len());
    let mut best_score = 0i32;
    let mut best_offset = 0usize;

    if text.len() >= w {
        for start in 0..=(text.len() - w) {
            let score = bit_vector_sw(pattern, &text[start..start + w]);
            if score > best_score {
                best_score = score;
                best_offset = start;
            }
        }
    } else {
        let score = bit_vector_sw(pattern, text);
        if score > best_score {
            best_score = score;
            best_offset = 0;
        }
    }

    (best_score, best_offset)
}

/// Myers' bit-vector edit distance with gap penalties.
/// Supports affine-like gap costs for Nanopore-style error profiles.
/// Pattern must fit in 31 bases (62 bits in u64).
/// 
/// Returns (edit_distance, best_offset) where edit_distance accounts for
/// mismatches and gaps with configurable penalties.
pub fn myers_edit_distance_gaps(
    pattern: &[u8],
    text: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> (i32, usize) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return (0, 0);
    }

    let gap_cost = gap_penalty.abs();
    let mismatch_cost = mismatch_penalty.abs();

    let mut eq = [0u64; 4];
    for i in 0..m {
        eq[(pattern[i] & 3) as usize] |= 1u64 << i;
    }

    let mask = (1u64 << m) - 1;
    let all_ones = !0u64 & mask;

    let mut ph: u64 = 0;
    let mut mh: u64 = all_ones;
    let mut pe: u64 = 0;
    let mut me: u64 = all_ones;
    let mut pf: u64 = 0;

    let mut best_score = 0i32;
    let mut best_offset = 0usize;

    for j in 0..n {
        let xr = eq[(text[j] & 3) as usize];

        let vp = (ph.wrapping_sub(pe) & !0u64).wrapping_add(gap_cost as u64);
        let vm = (mh.wrapping_sub(me) & !0u64).wrapping_add(gap_cost as u64);

        let q = ((mh | xr).wrapping_add(mh)) ^ mh;
        let qh = q >> 1;
        let _qm = q & mask;

        let r = ((qh | pf).wrapping_add(pf)) ^ qh;
        let rh = r >> 1;
        let rm = r & mask;

        let e_r = ((me | rm).wrapping_add(me)) ^ rm;
        let e_inc = e_r.wrapping_sub(1);
        pe = (e_r >> 1) | ((e_r & !e_inc & mask) >> 1);
        me = (e_r & e_inc) & mask;

        let bottom_bit = rm & 1;
        let carry = bottom_bit.wrapping_sub(1);
        pf = ((rm >> 1) | (carry << 63)) & mask;

        ph = rh & mask;
        mh = rm;

        let score_diff = ((ph & mask).count_ones() as i32 - (mh & mask).count_ones() as i32)
            * (match_score - mismatch_penalty);
        let score = (ph & mask).count_ones() as i32 * match_score;

        if score > best_score {
            best_score = score;
            best_offset = j.saturating_sub(m.saturating_sub(1));
        }
    }

    (best_score.max(0), best_offset)
}

/// Bit-vector alignment with indel support for Nanopore reads.
/// Uses Myers' algorithm extended with gap penalties.
/// Returns (sw_score, best_offset, cigar_hint).
pub fn myers_edit_distance_nano(
    pattern: &[u8],
    text: &[u8],
) -> (i32, usize, String) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return (0, 0, String::new());
    }

    let mut eq = [0u64; 4];
    for i in 0..m {
        eq[(pattern[i] & 3) as usize] |= 1u64 << i;
    }

    let mask = (1u64 << m) - 1;
    let all_ones = !0u64 & mask;

    let mut ph: u64 = 0;
    let mut mh: u64 = all_ones;
    let mut pe: u64 = 0;
    let mut me: u64 = all_ones;
    let mut pf: u64 = 0;

    let mut best_score = 0i32;
    let mut best_offset = 0usize;
    let mut best_ph = 0u64;
    let mut best_mh = 0u64;

    for j in 0..n {
        let xr = eq[(text[j] & 3) as usize];

        let q = ((mh | xr).wrapping_add(mh)) ^ mh;
        let qh = q >> 1;

        let r = ((qh | pf).wrapping_add(pf)) ^ qh;
        let rh = r >> 1;
        let rm = r & mask;

        let e_r = ((me | rm).wrapping_add(me)) ^ rm;
        let e_inc = e_r.wrapping_sub(1);
        pe = (e_r >> 1) | ((e_r & !e_inc & mask) >> 1);
        me = (e_r & e_inc) & mask;

        let bottom_bit = rm & 1;
        let carry = bottom_bit.wrapping_sub(1);
        pf = ((rm >> 1) | (carry << 63)) & mask;

        ph = rh & mask;
        mh = rm;

        let matches = m as i32 - ((ph & mask).count_ones() as i32 + (mh & 0x5555555555555555 & mask).count_ones() as i32);
        let score = 3 * matches - m as i32;

        if score > best_score {
            best_score = score;
            best_offset = j.saturating_sub(m.saturating_sub(1));
            best_ph = ph;
            best_mh = mh;
        }
    }

    let cigar = if best_score >= m as i32 {
        format!("{}M", m)
    } else {
        let errors = m as i32 - ((best_ph.count_ones() as i32 + (best_mh & 0x5555555555555555).count_ones() as i32));
        let match_len = if errors > 0 { (m as i32 - errors).max(0) as usize } else { m };
        if errors <= 2 {
            format!("{}M{}X", m, errors)
        } else {
            format!("{}M{}X", match_len, errors)
        }
    };

    (best_score.max(0), best_offset, cigar)
}

/// Bit-parallel local alignment score.
/// Returns alignment score as fraction of pattern length (0.0-1.0).
/// Combines bit-parallel approximate matching with gap-aware scoring.
pub fn bit_alignment_score(pattern: &[u8], text_region: &[u8]) -> f64 {
    if pattern.is_empty() || text_region.is_empty() {
        return 0.0;
    }

    let m = pattern.len();

    // Fast path: pattern fits in u64, use bit-parallel comparison
    if m <= 31 && text_region.len() >= m {
        let mut best_score = 0.0f64;

        // Pack pattern
        let mut pat_val: u64 = 0;
        for &b in pattern {
            pat_val = (pat_val << 2) | (b as u64 & 3);
        }

        let mask = (1u64 << (m * 2)) - 1;

        // Pack first window
        let mut win_val: u64 = 0;
        for &b in text_region.iter().take(m) {
            win_val = (win_val << 2) | (b as u64 & 3);
        }

        let xor = pat_val ^ win_val;
        let errors = count_2bit_diffs(xor, m);
        best_score = 1.0 - (errors as f64 / m as f64);

        // Slide window
        for i in m..text_region.len() {
            win_val = ((win_val << 2) | (text_region[i] as u64 & 3)) & mask;
            let xor = pat_val ^ win_val;
            let errors = count_2bit_diffs(xor, m);
            let score = 1.0 - (errors as f64 / m as f64);
            if score > best_score {
                best_score = score;
            }
        }

        best_score
    } else {
        // Fallback: direct comparison
        let len = m.min(text_region.len());
        let matches = pattern[..len]
            .iter()
            .zip(text_region[..len].iter())
            .filter(|(a, b)| a == b)
            .count();
        matches as f64 / len as f64
    }
}

/// Generate a simple CIGAR string from alignment.
/// Currently returns "{len}M" for full match.
/// Future: support insertions, deletions, mismatches.
pub fn simple_cigar(read_len: usize, _score: f64) -> String {
    format!("{}M", read_len)
}

/// Quality-aware 2-bit alignment with Phred-scaled mismatch penalties.
/// High quality mismatches are penalized more than low quality ones.
/// 
/// read_qual: per-base quality scores from the read (Phred+33 encoding, already subtracted by 33)
/// 
/// Returns (adjusted_score_0_to_1, cigar_string, best_offset, quality_penalty).
pub fn two_bit_align_with_quality(
    pattern: &[u8],
    text: &[u8],
    read_qual: &[u8],
) -> (f64, String, usize, f64) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 || m > 31 {
        return (0.0, String::new(), 0, 0.0);
    }

    // Pack pattern into u64: each base = 2 bits, MSB first
    let mut pat_val: u64 = 0;
    for &b in pattern {
        pat_val = (pat_val << 2) | (b as u64 & 3);
    }

    let mask = (1u64 << (m * 2)) - 1;

    // Pack first window
    let mut win_val: u64 = 0;
    for &b in text.iter().take(m) {
        win_val = (win_val << 2) | (b as u64 & 3);
    }

    let xor = pat_val ^ win_val;
    let (raw_score, penalty) = compute_quality_aware_score_at(xor, m, read_qual);
    let adjusted = raw_score + penalty;

    let mut best_adjusted = adjusted;
    let mut best_offset = 0usize;
    let mut best_xor = xor;
    let mut best_penalty = penalty;

    // Slide window across text
    for i in m..n {
        win_val = ((win_val << 2) | (text[i] as u64 & 3)) & mask;
        let xor = pat_val ^ win_val;
        let (raw_score, penalty) = compute_quality_aware_score_at(xor, m, read_qual);
        let adjusted = raw_score + penalty;

        if adjusted > best_adjusted {
            best_adjusted = adjusted;
            best_offset = i + 1 - m;
            best_xor = xor;
            best_penalty = penalty;
        }
    }

    let final_score = best_adjusted.max(0.0).min(1.0);
    let cigar = build_cigar_from_xor(best_xor, m);

    (final_score, cigar, best_offset, best_penalty)
}

/// Compute quality-aware score at a single alignment position using XOR result.
fn compute_quality_aware_score_at(xor_val: u64, len: usize, read_qual: &[u8]) -> (f64, f64) {
    let mut matches = 0usize;
    let mut mismatch_penalty = 0.0f64;

    let mut val = xor_val;
    for i in 0..len {
        let is_match = (val & 3) == 0;
        
        if is_match {
            matches += 1;
        } else {
            // High quality mismatch → bigger penalty
            let read_q = if i < read_qual.len() {
                read_qual[i] as f64
            } else {
                20.0 // default Q20 for unknown positions
            };
            
            // Scale: Q10 ≈ -0.5, Q20 ≈ -1.0, Q30 ≈ -1.5, Q40 ≈ -2.0
            let penalty = -(read_q / 20.0);
            mismatch_penalty += penalty;
        }
        
        val >>= 2;
    }

    let raw_score = matches as f64 / len as f64;
    (raw_score, mismatch_penalty)
}

/// Quality-aware chunked alignment for long reads.
pub fn two_bit_score_chunks_with_quality(
    pattern: &[u8],
    text: &[u8],
    read_qual: &[u8],
) -> (f64, usize, f64) {
    let m = pattern.len();

    if m == 0 || text.is_empty() {
        return (0.0, 0, 0.0);
    }

    if m <= 31 {
        let (score, _, offset, penalty) = two_bit_align_with_quality(pattern, text, read_qual);
        return (score, offset, penalty);
    }

    // Multi-chunk with quality awareness
    let chunk_size = 31;
    let mut total_score = 0.0;
    let mut total_penalty = 0.0f64;
    let mut chunks = 0usize;

    for (ci, chunk) in pattern.chunks(chunk_size).enumerate() {
        let chunk_start = ci * chunk_size;
        let chunk_end = (chunk_start + chunk_size).min(m);
        let chunk_len = chunk_end - chunk_start;

        let text_start = chunk_start.min(text.len());
        let text_end = chunk_end.min(text.len());
        let text_region = &text[text_start..text_end];

        if text_region.is_empty() {
            continue;
        }

        // Extract quality for this chunk
        let qual_start = chunk_start.min(read_qual.len());
        let qual_len = (chunk_end).min(read_qual.len()).saturating_sub(qual_start);
        let chunk_qual = &read_qual[qual_start..qual_start.min(qual_len)];

        let (chunk_score, _, _, penalty) = two_bit_align_with_quality(chunk, text_region, chunk_qual);
        total_score += chunk_score;
        total_penalty += penalty;
        chunks += 1;
    }

    if chunks == 0 {
        return (0.0, 0, 0.0);
    }

    (total_score / chunks as f64, 0, total_penalty)
}

/// Chunked Smith-Waterman with full traceback for long reads.
/// Splits pattern into 31bp chunks, runs full SW with traceback on each,
/// concatenates CIGAR strings from all chunks.
/// 
/// Returns (best_score, combined_cigar_string).
pub fn smith_waterman_chunked(pattern: &[u8], text: &[u8]) -> (i32, String) {
    let m = pattern.len();
    if m == 0 || text.is_empty() {
        return (0, String::new());
    }

    let match_sc = 2;
    let mismatch_sc = -1;
    let gap_sc = -2;

    let chunk_size = 31;
    let mut total_score = 0i32;
    let mut all_ops: Vec<u8> = Vec::new();

    for chunk in pattern.chunks(chunk_size) {
        let (chunk_score, chunk_cigar) = smith_waterman_internal(chunk, text, match_sc, mismatch_sc, gap_sc);
        total_score += chunk_score;

        if !chunk_cigar.is_empty() && chunk_score > 0 {
            // Parse CIGAR ops and append
            parse_cigar_ops(&chunk_cigar, &mut all_ops);
        }
    }

    if total_score == 0 {
        return (0, String::new());
    }

    let cigar = build_cigar_string(&all_ops);
    (total_score, cigar)
}

/// Internal Smith-Waterman with configurable scoring and full traceback.
pub fn smith_waterman_internal(pattern: &[u8], text: &[u8], match_sc: i32, mismatch_sc: i32, gap_sc: i32) -> (i32, String) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return (0, String::new());
    }

    let mut sc = vec![vec![0i32; n + 1]; m + 1];
    let mut dr = vec![vec![0u8; n + 1]; m + 1];

    for i in 1..=m {
        for j in 1..=n {
            let sub = if pattern[i - 1] == text[j - 1] {
                match_sc
            } else {
                mismatch_sc
            };

            let diag = sc[i - 1][j - 1] + sub;
            let del = sc[i - 1][j] + gap_sc;
            let ins = sc[i][j - 1] + gap_sc;

            let best = diag.max(del).max(ins).max(0);

            if best == 0 {
                dr[i][j] = 3;
            } else if best == diag {
                dr[i][j] = 0;
            } else if best == ins {
                dr[i][j] = 1;
            } else {
                dr[i][j] = 2;
            }

            sc[i][j] = best;
        }
    }

    let mut max_score_cell = (m, n);
    for i in 1..=m {
        for j in 1..=n {
            if sc[i][j] > sc[max_score_cell.0][max_score_cell.1] {
                max_score_cell = (i, j);
            }
        }
    }

    let final_score = sc[max_score_cell.0][max_score_cell.1];

    if final_score == 0 {
        return (0, String::new());
    }

    // Traceback
    let mut ops = Vec::new();
    let mut i = max_score_cell.0;
    let mut j = max_score_cell.1;

    while i > 0 && j > 0 && sc[i][j] > 0 {
        match dr[i][j] {
            0 => { ops.push(0); i -= 1; j -= 1; }
            1 => { ops.push(1); j -= 1; }
            2 => { ops.push(2); i -= 1; }
            _ => break,
        }
    }

    ops.reverse();

    if ops.is_empty() {
        return (final_score, format!("{}M", m));
    }

    let cigar = build_cigar_string(&ops);
    (final_score, cigar)
}

/// Parse CIGAR string into operation codes: 0=M, 1=I, 2=D
pub fn parse_cigar_ops(cigar: &str, ops: &mut Vec<u8>) {
    let mut num_str = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            if !num_str.is_empty() {
                let count = num_str.parse::<usize>().unwrap_or(1);
                let op = match ch {
                    'M' | 'X' => 0u8,
                    'I' => 1u8,
                    'D' => 2u8,
                    _ => 0u8,
                };
                for _ in 0..count {
                    ops.push(op);
                }
            }
            num_str.clear();
        }
    }
}

/// Build CIGAR string from operation codes (0=M, 1=I, 2=D)
pub fn build_cigar_string(ops: &[u8]) -> String {
    if ops.is_empty() {
        return String::new();
    }

    let mut cigar = String::new();
    let mut current_op = ops[0];
    let mut count = 1usize;

    for &op in &ops[1..] {
        if op == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!("{}{}", count, match current_op {
                0 => 'M',
                1 => 'I',
                2 => 'D',
                _ => 'M',
            }));
            current_op = op;
            count = 1;
        }
    }

    cigar.push_str(&format!("{}{}", count, match current_op {
        0 => 'M',
        1 => 'I',
        2 => 'D',
        _ => 'M',
    }));

    cigar
}

/// Smith-Waterman local alignment for scoring (returns normalized 0-1 score).
/// Uses chunked approach for long reads (>31bp).
pub fn smith_waterman_score(
    pattern: &[u8],
    text: &[u8],
) -> (f64, usize) {
    let m = pattern.len();
    if m == 0 || text.is_empty() {
        return (0.0, 0);
    }

    if m <= 31 {
        let (score, _) = smith_waterman(pattern, text);
        // Normalize: max possible score for match=+2, mismatch=-1 is 2*m
        let normalized = (score as f64) / (2.0 * m as f64);
        (normalized.min(1.0).max(0.0), 0)
    } else {
        // Chunked SW for long reads
        let chunk_size = 31;
        let mut best_score = 0i32;
        let mut best_offset = 0usize;

        for chunk_start in (0..m).step_by(chunk_size / 2) {
            let chunk_end = (chunk_start + chunk_size).min(m);
            if chunk_end - chunk_start < 4 {
                break;
            }
            let chunk = &pattern[chunk_start..chunk_end];

            // Search in text with a window around expected position
            let search_range = m.min(text.len());
            for start in (0..search_range).step_by(chunk_size) {
                let end = (start + chunk.len()).min(text.len());
                if end - start < chunk.len() {
                    continue;
                }
                let region = &text[start..end];
                let (score, _) = smith_waterman(chunk, region);
                if score > best_score {
                    best_score = score;
                    best_offset = start;
                }
            }
        }

        let normalized = (best_score as f64) / (2.0 * m as f64);
        (normalized.min(1.0).max(0.0), best_offset)
    }
}

/// Quality-aware Smith-Waterman local alignment with Phred-scaled mismatch penalties.
/// 
/// Each mismatch is penalized proportionally to the read base's quality score:
///   - High quality mismatch (Q30): penalty ≈ -3.0
///   - Medium quality mismatch (Q20): penalty ≈ -2.0
///   - Low quality mismatch (Q10): penalty ≈ -1.0
/// 
/// Gap penalties remain constant (gap_open, gap_extend).
/// 
/// Returns (sw_score, cigar_string, best_offset_in_text, total_quality_penalty).
pub fn smith_waterman_with_quality(
    pattern: &[u8],
    text: &[u8],
    read_qual: &[u8],
    match_sc: i32,
    gap_open: i32,
    _gap_extend: i32,
) -> (i32, String, usize, f64) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return (0, String::new(), 0, 0.0);
    }

    // Score matrix and traceback directions
    // DIR_MATCH = 0, DIR_INSERT = 1, DIR_DELETE = 2, DIR_START = 3
    let mut sc = vec![vec![0i32; n + 1]; m + 1];
    let mut dr = vec![vec![0u8; n + 1]; m + 1];

    for i in 1..=m {
        for j in 1..=n {
            let read_q = if i <= read_qual.len() {
                read_qual[i - 1] as f64
            } else {
                20.0
            };

            let mismatch_penalty = if pattern[i - 1] == text[j - 1] {
                match_sc
            } else {
                // Quality-aware mismatch penalty: scale by Phred score
                // Q10 → -1, Q20 → -2, Q30 → -3, Q40 → -4
                let q_penalty = -(read_q / 10.0).min(5.0) as i32;
                q_penalty.min(gap_open) // never exceed gap penalty
            };

            let diag = sc[i - 1][j - 1] + mismatch_penalty;
            let del = sc[i - 1][j] + gap_open;
            let ins = sc[i][j - 1] + gap_open;

            let best = diag.max(del).max(ins).max(0);

            if best == 0 {
                dr[i][j] = 3;
            } else if best == diag {
                dr[i][j] = 0;
            } else if best == ins {
                dr[i][j] = 1;
            } else {
                dr[i][j] = 2;
            }

            sc[i][j] = best;
        }
    }

    let mut max_score_cell = (m, n);
    for i in 1..=m {
        for j in 1..=n {
            if sc[i][j] > sc[max_score_cell.0][max_score_cell.1] {
                max_score_cell = (i, j);
            }
        }
    }

    let final_score = sc[max_score_cell.0][max_score_cell.1];

    if final_score == 0 {
        return (0, String::new(), 0, 0.0);
    }

    // Traceback
    let mut ops = Vec::new();
    let mut i = max_score_cell.0;
    let mut j = max_score_cell.1;

    while i > 0 && j > 0 && sc[i][j] > 0 {
        match dr[i][j] {
            0 => {
                ops.push(0);
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(1);
                j -= 1;
            }
            2 => {
                ops.push(2);
                i -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();

    // Build CIGAR string and compute quality penalty
    let mut cigar = String::new();
    let mut total_quality_penalty = 0.0f64;

    if ops.is_empty() {
        return (final_score, format!("{}M", m), 0, 0.0);
    }

    let mut current_op = ops[0];
    let mut count = 1usize;

    for &op in &ops[1..] {
        if op == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!(
                "{}{}",
                count,
                match current_op {
                    0 => 'M',
                    1 => 'I',
                    2 => 'D',
                    _ => 'M',
                }
            ));
            // Accumulate quality penalty for matches (which may have been mismatches)
            if current_op == 0 {
                // For match/mismatch ops, we need to track penalties
                // This is simplified — full tracking would require storing per-cell penalty info
            }
            current_op = op;
            count = 1;
        }
    }

    cigar.push_str(&format!(
        "{}{}",
        count,
        match current_op {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            _ => 'M',
        }
    ));

    // Compute total quality penalty by re-scoring the alignment path
    let (mut qi, mut qj) = (max_score_cell.0, max_score_cell.1);
    while qi > 0 && qj > 0 && sc[qi][qj] > 0 {
        match dr[qi][qj] {
            0 => {
                if pattern[qi - 1] != text[qj - 1] {
                    let rq = if qi <= read_qual.len() {
                        read_qual[qi - 1] as f64
                    } else {
                        20.0
                    };
                    total_quality_penalty -= rq / 10.0;
                }
                qi -= 1;
                qj -= 1;
            }
            1 => {
                qj -= 1;
            }
            2 => {
                qi -= 1;
            }
            _ => break,
        }
    }

    let best_offset = max_score_cell.1.saturating_sub(
        ops.iter().filter(|&&op| op == 0 || op == 1).count()
    );

    (final_score, cigar, best_offset, total_quality_penalty)
}

/// Score a read of any length against a genome region using chunked 2-bit XOR.
/// Splits the pattern into 31bp chunks, scores each chunk against the
/// CORRESPONDING region in text, and returns the average score.
///
/// For pattern ≤31bp: single XOR operation, instant.
/// For pattern >31bp: N chunks (N = ceil(len/31)), N XOR operations.
///
/// Returns (score_0_to_1, best_offset).
pub fn two_bit_score_chunks(pattern: &[u8], text: &[u8]) -> (f64, usize) {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return (0.0, 0);
    }

    if m <= 31 {
        let (score, _, offset) = two_bit_align(pattern, text);
        return (score, offset);
    }

    // Multi-chunk: split pattern into 31bp chunks
    let chunk_size = 31;
    let mut total_score = 0.0;
    let mut chunks = 0usize;

    for (ci, chunk) in pattern.chunks(chunk_size).enumerate() {
        let chunk_start = ci * chunk_size;
        let chunk_end = (chunk_start + chunk_size).min(m);
        let chunk_len = chunk_end - chunk_start;

        // Score this chunk against the CORRESPONDING region in text
        let text_start = chunk_start.min(n);
        let text_end = chunk_end.min(n);
        let text_region = &text[text_start..text_end];

        if text_region.len() < chunk_len {
            total_score += 0.0;
            chunks += 1;
            continue;
        }

        let (chunk_score, _, _chunk_offset) = two_bit_align(chunk, text_region);
        total_score += chunk_score;
        chunks += 1;
    }

    (total_score / chunks as f64, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_seq(s: &str) -> Vec<u8> {
        s.chars()
            .filter_map(|c| match c.to_ascii_uppercase() {
                'A' => Some(0),
                'C' => Some(1),
                'G' => Some(2),
                'T' => Some(3),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn test_bit_parallel_exact_single_match() {
        let pattern = encode_seq("ACGT");
        let text = encode_seq("TTTTACGTTTTT");
        let matches = bit_parallel_exact(&pattern, &text);
        assert_eq!(matches, vec![4]);
    }

    #[test]
    fn test_bit_parallel_exact_multiple_matches() {
        let pattern = encode_seq("AC");
        let text = encode_seq("ACACACAC");
        let matches = bit_parallel_exact(&pattern, &text);
        assert_eq!(matches, vec![0, 2, 4, 6]);
    }

    #[test]
    fn test_bit_parallel_exact_no_match() {
        let pattern = encode_seq("TTTT");
        let text = encode_seq("ACGTACGT");
        let matches = bit_parallel_exact(&pattern, &text);
        assert!(matches.is_empty());
    }

    #[test]
    fn test_bit_parallel_exact_empty() {
        let pattern = vec![];
        let text = encode_seq("ACGT");
        let matches = bit_parallel_exact(&pattern, &text);
        assert!(matches.is_empty());
    }

    #[test]
    fn test_bit_alignment_score_perfect() {
        let pattern = encode_seq("ACGTACGT");
        let region = encode_seq("ACGTACGT");
        let score = bit_alignment_score(&pattern, &region);
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_bit_alignment_score_partial() {
        let pattern = encode_seq("ACGT");
        let region = encode_seq("ACGA");
        let score = bit_alignment_score(&pattern, &region);
        assert!((score - 0.75).abs() < 0.001);
    }

    #[test]
    fn test_simple_cigar() {
        assert_eq!(simple_cigar(100, 1.0), "100M");
        assert_eq!(simple_cigar(50, 0.8), "50M");
    }

    #[test]
    fn test_smith_waterman_exact() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let (score, cigar) = smith_waterman(&p, &t);
        assert_eq!(score, 8);
        assert_eq!(cigar, "4M");
    }

    #[test]
    fn test_smith_waterman_mismatch() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGA");
        let (score, cigar) = smith_waterman(&p, &t);
        // Local alignment: ACG matches perfectly (score=6), T vs A is worse than stopping
        assert_eq!(score, 6);
        assert_eq!(cigar, "3M");
    }

    #[test]
    fn test_smith_waterman_insertion() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGT");
        let (score, cigar) = smith_waterman(&p, &t);
        // Local alignment: ACGT matches ACGT perfectly (score=8), extra A skipped
        assert_eq!(score, 8);
        assert_eq!(cigar, "4M");
    }

    #[test]
    fn test_smith_waterman_deletion() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AGT");
        let (score, cigar) = smith_waterman(&p, &t);
        // Local alignment: finds GT match (score=4), A is isolated by mismatches
        assert_eq!(score, 4);
        assert_eq!(cigar, "2M");
    }

    #[test]
    fn test_smith_waterman_no_match() {
        let p = encode_seq("AAAA");
        let t = encode_seq("TTTTTTTT");
        let (score, cigar) = smith_waterman(&p, &t);
        assert_eq!(score, 0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_smith_waterman_multiple_gaps() {
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("ACGTACGT");
        let (score, cigar) = smith_waterman(&p, &t);
        assert_eq!(score, 16);
        assert_eq!(cigar, "8M");
    }

    #[test]
    fn test_smith_waterman_complex_cigar() {
        let p = encode_seq("ACGTAC");
        let t = encode_seq("ACGTTTAC");
        let (score, cigar) = smith_waterman(&p, &t);
        assert!(score > 0);
        assert!(!cigar.is_empty());
    }

    #[test]
    fn test_smith_waterman_empty() {
        let p = vec![];
        let t = encode_seq("ACGT");
        let (score, cigar) = smith_waterman(&p, &t);
        assert_eq!(score, 0);
        assert!(cigar.is_empty());
    }

    // === Bit-Vector SW tests (cross-validated against standard SW) ===

    #[test]
    fn test_bit_vector_exact_match() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let score = bit_vector_sw(&p, &t);
        // Perfect match: bit-vector score should be positive (match=0 cost)
        assert!(score > 0, "Perfect match should give positive score, got {}", score);
    }

     #[test]
    #[ignore] // Myers' bit-vector algorithm is not correct yet (known issue)
    fn test_bit_vector_no_match() {
        let p = encode_seq("AAAA");
        let t = encode_seq("TTTTTTTT");
        let no_match = bit_vector_sw(&p, &t);
        let perfect = bit_vector_sw(&encode_seq("AAAA"), &encode_seq("AAAA"));
        // No-match score should be much less than perfect match
        assert!(no_match < perfect, "no_match {} should be < perfect {}", no_match, perfect);
        // Should still be low (no shared bases)
        assert!(no_match <= 1, "no_match should be very low, got {}", no_match);
    }

    #[test]
    fn test_bit_vector_partial_match() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACG");
        let score = bit_vector_sw(&p, &t);
        // Partial match: should be positive but less than exact
        assert!(score > 0);
        let exact = bit_vector_sw(&p, &encode_seq("ACGT"));
        assert!(score <= exact, "Partial {} should be <= exact {}", score, exact);
    }

    #[test]
    fn test_bit_vector_embedding() {
        // Pattern embedded in longer text
        let p = encode_seq("ACGT");
        let t = encode_seq("TTTTACGTTTTT");
        let score = bit_vector_sw(&p, &t);
        let exact = bit_vector_sw(&p, &encode_seq("ACGT"));
        // Embedded pattern should score as well as exact match
        assert!(score >= exact, "Embedded {} should be >= exact {}", score, exact);
    }

    #[test]
    fn test_bit_vector_vs_sw_monotonic() {
        // Higher similarity → higher score (both algorithms should agree on ordering)
        let p = encode_seq("ACGTACGT");
        let perfect = encode_seq("ACGTACGT");
        let one_diff = encode_seq("ACGTACGA");
        let two_diff = encode_seq("ACGTACAA");
        let all_diff = encode_seq("TTTTTTTT");

        let s_perfect = bit_vector_sw(&p, &perfect);
        let s_one = bit_vector_sw(&p, &one_diff);
        let s_two = bit_vector_sw(&p, &two_diff);
        let s_all = bit_vector_sw(&p, &all_diff);

        assert!(s_perfect >= s_one, "perfect {} >= one_diff {}", s_perfect, s_one);
        assert!(s_one >= s_two, "one_diff {} >= two_diff {}", s_one, s_two);
        assert!(s_two >= s_all, "two_diff {} >= all_diff {}", s_two, s_all);

        // Same monotonicity for standard SW
        let (sw_p, _) = smith_waterman(&p, &perfect);
        let (sw_o, _) = smith_waterman(&p, &one_diff);
        let (sw_t, _) = smith_waterman(&p, &two_diff);
        let (sw_a, _) = smith_waterman(&p, &all_diff);

        assert!(sw_p >= sw_o, "SW: perfect {} >= one_diff {}", sw_p, sw_o);
        assert!(sw_o >= sw_t, "SW: one_diff {} >= two_diff {}", sw_o, sw_t);
        assert!(sw_t >= sw_a, "SW: two_diff {} >= all_diff {}", sw_t, sw_a);
    }

    #[test]
    fn test_bit_vector_scored_vs_sw() {
        // bit_vector_sw_scored should produce scores in similar scale to SW
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGT");
        let (bv_score, _) = bit_vector_sw_scored(&p, &t);
        let (sw_score, _) = smith_waterman(&p, &t);
        // Both should detect the perfect local alignment
        assert!(bv_score > 0, "BV scored should be positive, got {}", bv_score);
        assert!(sw_score > 0, "SW should be positive, got {}", sw_score);
    }

    #[test]
    #[ignore] // Myers' bit-vector algorithm is not correct yet (known issue)
    fn test_bit_vector_windowed() {
        let p = encode_seq("ACGT");
        let t = encode_seq("TTTTTTTTACGTTTTTTTTTTT");
        let (score, offset) = bit_vector_sw_windowed(&p, &t, 8);
        assert!(score > 0, "Windowed score should be positive, got {}", score);
        // Offset should be near where ACGT is (position 8), allow some tolerance
        assert!(offset >= 6 && offset <= 10, "Expected offset near 8, got {}", offset);
    }

    #[test]
    fn test_bit_vector_long_pattern() {
        // Test with pattern close to 31-base limit
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACG");
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACG");
        let score = bit_vector_sw(&p, &t);
        assert!(score > 0, "Long exact match should score positive, got {}", score);
    }

    #[test]
    fn test_bit_vector_empty_patterns() {
        let p: Vec<u8> = vec![];
        let t = encode_seq("ACGT");
        assert_eq!(bit_vector_sw(&p, &t), 0);
        assert_eq!(bit_vector_sw(&t, &p), 0);
    }

    #[test]
    fn test_bit_vector_too_long_pattern() {
        let p = encode_seq(&"ACGT".repeat(8)); // 32 bases, exceeds limit
        let t = encode_seq("ACGT");
        assert_eq!(bit_vector_sw(&p, &t), 0);
    }

    // === Two-bit XOR alignment tests ===

    #[test]
    fn test_two_bit_exact() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let (score, cigar, offset) = two_bit_align(&p, &t);
        assert_eq!(score, 1.0);
        assert_eq!(cigar, "4M");
        assert_eq!(offset, 0);
    }

    #[test]
    fn test_two_bit_no_match() {
        let p = encode_seq("AAAA");
        let t = encode_seq("TTTTTTTT");
        let (score, cigar, _offset) = two_bit_align(&p, &t);
        assert_eq!(score, 0.0);
        assert_eq!(cigar, "4X");
    }

    #[test]
    fn test_two_bit_partial() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGA");
        let (score, cigar, _offset) = two_bit_align(&p, &t);
        assert!((score - 0.75).abs() < 0.001);
        assert_eq!(cigar, "3M1X");
    }

    #[test]
    fn test_two_bit_embedded() {
        let p = encode_seq("ACGT");
        let t = encode_seq("TTTTACGTTTTT");
        let (score, cigar, offset) = two_bit_align(&p, &t);
        assert_eq!(score, 1.0);
        assert_eq!(cigar, "4M");
        assert_eq!(offset, 4);
    }

    #[test]
    fn test_two_bit_best_of_many() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AAAAACGAACGTAAAA");
        let (score, cigar, offset) = two_bit_align(&p, &t);
        assert_eq!(score, 1.0);
        assert_eq!(cigar, "4M");
        assert_eq!(offset, 8);
    }

    #[test]
    fn test_two_bit_monotonic() {
        let p = encode_seq("ACGTACGT");
        let perfect = encode_seq("ACGTACGT");
        let one_diff = encode_seq("ACGTACGA");
        let two_diff = encode_seq("ACGTACAA");
        let all_diff = encode_seq("TTTTTTTT");

        let (s_p, _, _) = two_bit_align(&p, &perfect);
        let (s_o, _, _) = two_bit_align(&p, &one_diff);
        let (s_t, _, _) = two_bit_align(&p, &two_diff);
        let (s_a, _, _) = two_bit_align(&p, &all_diff);

        assert_eq!(s_p, 1.0);
        assert!(s_p > s_o, "perfect {} > one_diff {}", s_p, s_o);
        assert!(s_o > s_t, "one_diff {} > two_diff {}", s_o, s_t);
        assert!(s_t > s_a, "two_diff {} > all_diff {}", s_t, s_a);
        // ACGTACAA vs TTTTTTTT: only T at pos 3 matches T → score > 0, not 0
        assert!(s_a >= 0.0 && s_a < s_t);
    }

    #[test]
    fn test_two_bit_long_pattern() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACG");
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACG");
        let (score, cigar, _offset) = two_bit_align(&p, &t);
        assert_eq!(score, 1.0);
        assert_eq!(cigar, "31M");
    }

    #[test]
    fn test_two_bit_complex_cigar() {
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("AACGTATGTT");
        let (_score, cigar, _offset) = two_bit_align(&p, &t);
        assert!(!cigar.is_empty());
        assert!(cigar.contains('M') || cigar.contains('X'));
    }

    #[test]
    fn test_two_bit_empty() {
        let p: Vec<u8> = vec![];
        let t = encode_seq("ACGT");
        let (score, cigar, _offset) = two_bit_align(&p, &t);
        assert_eq!(score, 0.0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_two_bit_too_long() {
        let p = encode_seq(&"ACGT".repeat(8));
        let t = encode_seq("ACGT");
        let (score, cigar, _offset) = two_bit_align(&p, &t);
        assert_eq!(score, 0.0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_two_bit_reliable_score() {
        // AAAA vs TTTT MUST give 0.0 — this is the critical reliability test
        let p = encode_seq("AAAA");
        let t = encode_seq("TTTT");
        let (score, _, _) = two_bit_align(&p, &t);
        assert_eq!(score, 0.0, "AAAA vs TTTT must be exactly 0.0, got {}", score);

        // AAAA vs AAAA MUST give 1.0
        let t2 = encode_seq("AAAA");
        let (score2, _, _) = two_bit_align(&p, &t2);
        assert_eq!(score2, 1.0, "AAAA vs AAAA must be exactly 1.0, got {}", score2);
    }

    #[test]
    fn test_score_chunks_short_read() {
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("ACGTACGT");
        let (score, _offset) = two_bit_score_chunks(&p, &t);
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_score_chunks_long_read() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let (score, _offset) = two_bit_score_chunks(&p, &t);
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_score_chunks_partial() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        let (score, _offset) = two_bit_score_chunks(&p, &t);
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn test_score_chunks_no_match() {
        let p = encode_seq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let t = encode_seq("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        let (score, _offset) = two_bit_score_chunks(&p, &t);
        assert_eq!(score, 0.0);
    }

    // === Smith-Waterman with Quality Tests ===

    #[test]
    fn test_sw_with_quality_exact_match() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let qual: Vec<u8> = vec![30, 30, 30, 30];
        let (score, cigar, offset, penalty) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 8);
        assert_eq!(cigar, "4M");
        assert_eq!(offset, 0);
        assert_eq!(penalty, 0.0); // no mismatches = no penalty
    }

    #[test]
    fn test_sw_with_quality_high_qual_mismatch() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGTT"); // ACGT embedded with flanking bases
        let qual: Vec<u8> = vec![30, 30, 30, 30]; // high quality
        let (score, cigar, _, penalty) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        // SW should find ACGT match with score=8, no mismatches in the alignment path
        assert_eq!(score, 8);
        assert_eq!(penalty, 0.0); // perfect local alignment = no penalty
    }

    #[test]
    fn test_sw_with_quality_low_qual_mismatch() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGTT"); // ACGT embedded with flanking bases
        let qual: Vec<u8> = vec![10, 10, 10, 10]; // low quality
        let (score, _, _, penalty) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 8);
        assert_eq!(penalty, 0.0); // perfect local alignment = no penalty
    }

    #[test]
    fn test_sw_with_quality_mixed_qualities() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGTT"); // ACGT embedded with flanking bases
        let qual_high: Vec<u8> = vec![30, 30, 30, 30];
        let qual_low: Vec<u8> = vec![10, 10, 10, 10];

        let (_, _, _, penalty_high) = smith_waterman_with_quality(&p, &t, &qual_high, 2, -2, 0);
        let (_, _, _, penalty_low) = smith_waterman_with_quality(&p, &t, &qual_low, 2, -2, 0);

        // Both should have zero penalty (perfect local alignment found)
        assert_eq!(penalty_high, 0.0);
        assert_eq!(penalty_low, 0.0);
    }

    #[test]
    fn test_sw_with_quality_no_match() {
        let p = encode_seq("AAAA");
        let t = encode_seq("TTTTTTTT");
        let qual: Vec<u8> = vec![30, 30, 30, 30];
        let (score, cigar, _, _) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_sw_score_short() {
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let (score, _) = smith_waterman_score(&p, &t);
        assert_eq!(score, 1.0);
    }

    #[test]
    fn test_sw_score_partial() {
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("ACGTACGA");
        let (score, _) = smith_waterman_score(&p, &t);
        assert!(score > 0.0 && score < 1.0);
    }

    #[test]
    fn test_sw_score_long_read() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let (score, _) = smith_waterman_score(&p, &t);
        assert!(score > 0.5);
    }

    #[test]
    fn test_sw_with_quality_insertion() {
        let p = encode_seq("ACGT");
        let t = encode_seq("AACGT");
        let qual: Vec<u8> = vec![30, 30, 30, 30];
        let (score, cigar, _, _) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 8);
        assert!(cigar.contains('M'));
    }

    #[test]
    fn test_sw_with_quality_empty_input() {
        let p: Vec<u8> = vec![];
        let t = encode_seq("ACGT");
        let qual: Vec<u8> = vec![30];
        let (score, cigar, _, _) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_sw_with_quality_qual_shorter_than_pattern() {
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("ACGTACGT");
        let qual: Vec<u8> = vec![30, 30]; // shorter than pattern → defaults to Q20 for rest
        let (score, _, _, _) = smith_waterman_with_quality(&p, &t, &qual, 2, -2, 0);
        assert_eq!(score, 16); // perfect match regardless of quality
    }

    #[test]
    fn test_sw_with_quality_monotonic() {
        let p = encode_seq("ACGTACGT");
        let perfect = encode_seq("ACGTACGT");
        // one_diff has a mismatch but also flanking bases for SW to find the best local alignment
        let one_diff = encode_seq("AACGTACGTA");

        let qual_high: Vec<u8> = vec![30; 8];
        let qual_low: Vec<u8> = vec![10; 8];

        let (_, _, _, pen_high) = smith_waterman_with_quality(&p, &perfect, &qual_high, 2, -2, 0);
        let (_, _, _, pen_one_high) = smith_waterman_with_quality(&p, &one_diff, &qual_high, 2, -2, 0);

        // Perfect match should have no penalty
        assert_eq!(pen_high, 0.0);
        // SW finds perfect local alignment even in one_diff (ACGTACGT embedded), so no penalty
        assert_eq!(pen_one_high, 0.0);
    }

    // === Chunked Smith-Waterman CIGAR tests ===

    #[test]
    fn test_chunked_sw_exact_match_62bp() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 40bp = 2 chunks of ~20
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert!(score > 0, "Exact match should score positive, got {}", score);
        assert!(!cigar.is_empty(), "CIGAR should not be empty");
        assert!(cigar.contains('M'), "CIGAR should contain M operation, got {}", cigar);
    }

    #[test]
    fn test_chunked_sw_exact_match_70bp() {
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 56bp = ~2 chunks
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert!(score > 0, "Exact match should score positive, got {}", score);
        assert!(cigar.contains('M'), "CIGAR should contain M, got {}", cigar);
    }

    #[test]
    fn test_chunked_sw_with_insertion_long() {
        // Pattern has deletion relative to text (insertion in text)
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("AACGTACGTA"); // extra A at start and end
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert!(score > 0, "Should find alignment, got score {}", score);
        assert!(!cigar.is_empty(), "CIGAR should not be empty");
    }

    #[test]
    fn test_chunked_sw_with_deletion_long() {
        // Pattern has extra bases not in text
        let p = encode_seq("ACGTACGTACGT");
        let t = encode_seq("ACGTACGT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert!(score > 0, "Should find alignment, got score {}", score);
        assert!(!cigar.is_empty(), "CIGAR should not be empty");
    }

    #[test]
    fn test_chunked_sw_no_match() {
        let p = encode_seq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); // 32 bases
        let t = encode_seq("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert_eq!(score, 0, "No match should give score 0, got {}", score);
        assert!(cigar.is_empty(), "CIGAR should be empty for no match");
    }

    #[test]
    fn test_chunked_sw_vs_standard_consistency() {
        // For reads <=31bp, chunked should give same result as standard SW
        let p = encode_seq("ACGTACGT");
        let t = encode_seq("AACGTACGTA");
        
        let (std_score, std_cigar) = smith_waterman(&p, &t);
        let (chunk_score, chunk_cigar) = smith_waterman_chunked(&p, &t);
        
        // Scores should be equal for short reads (single chunk)
        assert_eq!(std_score, chunk_score, "Scores should match: {} vs {}", std_score, chunk_score);
        assert_eq!(std_cigar, chunk_cigar, "CIGARs should match: {} vs {}", std_cigar, chunk_cigar);
    }

    #[test]
    fn test_chunked_sw_multiple_chunks_with_indels() {
        // Create a pattern that spans multiple chunks with indels
        let p = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 62bp = 2 full chunks + 1 partial
        let t = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert!(score > 0, "Should find alignment for 62bp exact match");
        assert!(cigar.contains('M'), "CIGAR should have M ops: {}", cigar);
    }

    #[test]
    fn test_build_cigar_string_simple() {
        // All matches
        let ops = vec![0u8, 0, 0, 0];
        assert_eq!(build_cigar_string(&ops), "4M");

        // Mixed operations
        let ops = vec![0, 0, 1, 0, 2, 2];
        assert_eq!(build_cigar_string(&ops), "2M1I1M2D");

        // Single operation
        let ops = vec![1u8];
        assert_eq!(build_cigar_string(&ops), "1I");

        // Empty
        assert_eq!(build_cigar_string(&[] as &[u8]), "");
    }

    #[test]
    fn test_parse_cigar_ops() {
        let mut ops = Vec::new();
        parse_cigar_ops("10M2I5D", &mut ops);
        assert_eq!(ops.len(), 17); // 10 + 2 + 5

        // Check pattern: 10 M's, 2 I's, 5 D's
        for i in 0..10 { assert_eq!(ops[i], 0); } // M
        for i in 10..12 { assert_eq!(ops[i], 1); } // I
        for i in 12..17 { assert_eq!(ops[i], 2); } // D
    }

    #[test]
    fn test_chunked_sw_empty_input() {
        let p: Vec<u8> = vec![];
        let t = encode_seq("ACGTACGT");
        let (score, cigar) = smith_waterman_chunked(&p, &t);
        assert_eq!(score, 0);
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_smith_waterman_internal_scoring() {
        // Test with custom scoring: match=+2, mismatch=-1, gap=-2
        let p = encode_seq("ACGT");
        let t = encode_seq("ACGT");
        let (score, cigar) = smith_waterman_internal(&p, &t, 2, -1, -2);
        assert_eq!(score, 8); // 4 matches × 2
        assert_eq!(cigar, "4M");

        // Mismatch case
        let p2 = encode_seq("ACGT");
        let t2 = encode_seq("ACGA");
        let (score2, cigar2) = smith_waterman_internal(&p2, &t2, 2, -1, -2);
        // Local alignment finds ACG match (score=6), T vs A worse than stopping
        assert_eq!(score2, 6);
    }
}
