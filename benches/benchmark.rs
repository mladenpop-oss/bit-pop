use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use bit_pop::BitPop;
use bit_pop::align::{two_bit_align, smith_waterman, bit_parallel_exact, myers_edit_distance_gaps, myers_edit_distance_nano};
use bit_pop::fm::FmIndex;
use bit_pop::serialize::{serialize_bitpop, deserialize_bitpop};

// --- Helpers ---

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

fn encode_dna(s: &str) -> Vec<u8> {
    encode_seq(s)
}

fn pseudo_random_sequence(len: usize) -> String {
    let mut s = String::with_capacity(len);
    for i in 0..len {
        s.push(match (i * 7 + i.wrapping_mul(31) * 13) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    s
}

// --- Alignment benches (unchanged, proven fast) ---

fn bench_two_bit_align(c: &mut Criterion) {
    let mut group = c.benchmark_group("two_bit_align");

    let read_50 = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let genome_region = encode_seq("TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
    group.bench_function("50bp_exact", |b| b.iter(|| two_bit_align(&read_50, &genome_region)));

    let read_100 = encode_seq(&pseudo_random_sequence(100));
    let genome_100 = encode_seq(&format!("TTTTT{}", pseudo_random_sequence(100)));
    group.bench_function("100bp_exact", |b| b.iter(|| two_bit_align(&read_100, &genome_100)));

    let read_150 = encode_seq(&pseudo_random_sequence(150));
    let genome_150 = encode_seq(&format!("TTTTT{}", pseudo_random_sequence(150)));
    group.bench_function("150bp_exact", |b| b.iter(|| two_bit_align(&read_150, &genome_150)));

    let read_mismatch = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAACA");
    group.bench_function("50bp_1pct_mismatch", |b| b.iter(|| two_bit_align(&read_mismatch, &genome_region)));

    let read_none = encode_seq(&"T".repeat(50));
    group.bench_function("50bp_no_match", |b| b.iter(|| two_bit_align(&read_none, &genome_region)));

    group.finish();
}

fn bench_smith_waterman(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman");

    let read_50 = encode_seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let genome_region = encode_seq("TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
    group.bench_function("50bp_exact", |b| b.iter(|| smith_waterman(&read_50, &genome_region)));

    let read_100 = encode_seq(&pseudo_random_sequence(100));
    let genome_100 = encode_seq(&format!("TTTTT{}", pseudo_random_sequence(100)));
    group.bench_function("100bp_exact", |b| b.iter(|| smith_waterman(&read_100, &genome_100)));

    group.finish();
}

fn bench_bit_parallel_exact(c: &mut Criterion) {
    let mut group = c.benchmark_group("bit_parallel_exact");

    let text = encode_seq(&pseudo_random_sequence(10000));
    let pattern = encode_seq(&pseudo_random_sequence(8));
    group.bench_function("8bp_in_10kb", |b| b.iter(|| bit_parallel_exact(&pattern, &text)));

    let pattern_31 = encode_seq(&pseudo_random_sequence(31));
    group.bench_function("31bp_in_10kb", |b| b.iter(|| bit_parallel_exact(&pattern_31, &text)));

    group.finish();
}

fn bench_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("comparison");

    for size in [31, 50, 100, 150, 300] {
        let read = encode_seq(&pseudo_random_sequence(size));
        let genome = encode_seq(&format!("TTTTT{}", pseudo_random_sequence(size)));

        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(format!("{size}bp_two_bit"), &(), |b, _| b.iter(|| two_bit_align(&read, &genome)));
        group.bench_with_input(format!("{size}bp_smith_waterman"), &(), |b, _| b.iter(|| smith_waterman(&read, &genome)));
    }

    group.finish();
}

// --- Myers vs Smith-Waterman benchmarks ---

fn bench_myers_vs_sw(c: &mut Criterion) {
    let mut group = c.benchmark_group("myers_vs_sw");

    for size in [15, 20, 30] {
        let pattern = encode_seq(&pseudo_random_sequence(size));
        let text = encode_seq(&format!("TTTTT{}", pseudo_random_sequence(size + 20)));

        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(format!("{size}bp_myers_gaps"), &pattern.clone(), |b, p| {
            b.iter(|| myers_edit_distance_gaps(p, &text, 2, -1, -2))
        });
        group.bench_with_input(format!("{size}bp_myers_nano"), &pattern.clone(), |b, p| {
            b.iter(|| myers_edit_distance_nano(p, &text))
        });
        group.bench_with_input(format!("{size}bp_smith_waterman"), &pattern.clone(), |b, p| {
            b.iter(|| smith_waterman(p, &text))
        });
    }

    group.finish();
}

fn bench_myers_accuracy(c: &mut Criterion) {
    let mut group = c.benchmark_group("myers_accuracy");

    let test_cases = [
        ("perfect", "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"),
        ("1_mismatch", "ACGTACGTACGTACGTACGT", "ACGTACGTTCGTACGTACGT"),
        ("2_mismatches", "ACGTACGTACGTACGTACGT", "ACGTTCGTTCAAGTACGTACGT"),
        ("1_ins", "ACGTACGTACGTACGTACGT", "ACGTACGTAACGTACGTACGT"),
        ("1_del", "ACGTACGTACGTACGTACGT", "ACGTACGTCGTACGTACGT"),
    ];

    for (name, pat_str, text_str) in test_cases {
        let pattern = encode_seq(pat_str);
        let text = encode_seq(text_str);

        group.bench_function(name, |b| {
            b.iter(|| {
                let myers_score = myers_edit_distance_gaps(&pattern, &text, 2, -1, -2);
                let (sw_score, _) = smith_waterman(&pattern, &text);
                (myers_score, sw_score)
            })
        });
    }

    group.finish();
}

// --- FM-Index benches (NEW) ---

fn bench_fm_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("fm_build");

    let sizes = [
        ("10kb", 10_000),
        ("100kb", 100_000),
        ("1mb", 1_000_000),
        ("10mb", 10_000_000),
    ];

    for (label, size) in &sizes {
        let seq = encode_dna(&pseudo_random_sequence(*size));
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_function(format!("{label}"), |b| {
            b.iter(|| FmIndex::build(&[("test", &seq)]))
        });
    }

    group.finish();
}

fn bench_fm_backward_search(c: &mut Criterion) {
    let mut group = c.benchmark_group("fm_backward_search");

    let genome = encode_dna(&pseudo_random_sequence(1_000_000));
    let idx = FmIndex::build(&[("1mb", &genome)]);

    let pattern_lens = [
        ("4bp", 4),
        ("8bp", 8),
        ("15bp", 15),
        ("31bp", 31),
    ];

    for (label, plen) in &pattern_lens {
        let pattern = encode_dna(&pseudo_random_sequence(*plen));
        group.throughput(Throughput::Elements(*plen as u64));
        group.bench_function(format!("{label}_1mb_genome"), |b| {
            b.iter(|| idx.backward_search(&pattern))
        });
    }

    group.finish();
}

fn bench_fm_count_occurrences(c: &mut Criterion) {
    let mut group = c.benchmark_group("fm_count_occurrences");

    let genome = encode_dna(&pseudo_random_sequence(1_000_000));
    let idx = FmIndex::build(&[("1mb", &genome)]);

    let pattern = encode_dna(&pseudo_random_sequence(15));
    group.bench_function("15bp_1mb_genome", |b| {
        b.iter(|| idx.count_occurrences(&pattern))
    });

    group.finish();
}

// --- kmer_filter benches (FIXED: pseudo-random instead of repetitive) ---

fn bench_kmer_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_filter");

    // Small genome (10KB)
    let mut bp = BitPop::new(6);
    bp.add_genome("small", &pseudo_random_sequence(10_000));
    bp.build();
    let read = &pseudo_random_sequence(20);
    group.bench_function("10kb_built", |b| b.iter(|| bp.kmer_filter(read)));

    // Medium genome (100KB) - NOTE: find_positions is O(n) for each k-mer
    // This is expected to be slower than the old inverted index for large genomes
    // The anchor_filter is the optimized path that avoids full position materialization
    let mut bp2 = BitPop::new(8);
    bp2.add_genome("medium", &pseudo_random_sequence(100_000));
    bp2.build();
    let read2 = &pseudo_random_sequence(50);
    group.bench_function("100kb_built", |b| b.iter(|| bp2.kmer_filter(read2)));

    group.finish();
}

// --- kmer_filter with threshold ---

fn bench_kmer_filter_threshold(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_filter_threshold");

    let mut bp = BitPop::new(8);
    bp.add_genome("test", &pseudo_random_sequence(100_000));
    bp.build();

    let read = &pseudo_random_sequence(50);

    group.bench_function("100kb_no_threshold", |b| b.iter(|| bp.kmer_filter(read)));
    group.bench_function("100kb_threshold_10", |b| b.iter(|| bp.kmer_filter_with_threshold(read, 10)));
    group.bench_function("100kb_threshold_100", |b| b.iter(|| bp.kmer_filter_with_threshold(read, 100)));

    group.finish();
}

// --- anchor_filter benches (NEW) ---

fn bench_anchor_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("anchor_filter");

    // 1MB single genome
    let mut bp = BitPop::new(8);
    bp.add_genome("1mb", &pseudo_random_sequence(1_000_000));
    bp.build();
    let read = &pseudo_random_sequence(50);
    group.bench_function("50bp_1mb", |b| b.iter(|| bp.anchor_filter(read, 0.5)));

    // Multi-genome (5 x 200KB)
    let mut bp_multi = BitPop::new(8);
    for i in 0..5 {
        bp_multi.add_genome(&format!("g{i}"), &pseudo_random_sequence(200_000));
    }
    bp_multi.build();
    group.bench_function("50bp_5genomes_1mb_total", |b| b.iter(|| bp_multi.anchor_filter(read, 0.5)));

    group.finish();
}

fn bench_anchor_filter_threshold(c: &mut Criterion) {
    let mut group = c.benchmark_group("anchor_filter_threshold");

    let mut bp = BitPop::new(8);
    bp.add_genome("test", &pseudo_random_sequence(1_000_000));
    bp.build();
    let read = &pseudo_random_sequence(50);

    group.bench_function("1mb_no_threshold", |b| b.iter(|| bp.anchor_filter(read, 0.5)));
    group.bench_function("1mb_threshold_100", |b| b.iter(|| bp.anchor_filter_with_threshold(read, 0.5, 100)));
    group.bench_function("1mb_threshold_1000", |b| b.iter(|| bp.anchor_filter_with_threshold(read, 0.5, 1000)));

    group.finish();
}

// --- full_pipeline benches (FIXED: pseudo-random) ---

fn bench_full_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_pipeline");

    // Single genome (10KB)
    let mut bp = BitPop::new(6);
    bp.add_genome("genome1", &pseudo_random_sequence(10_000));
    bp.build();
    let read = &pseudo_random_sequence(20);
    group.bench_function("20bp_1gen_10kb", |b| b.iter(|| bp.map_read(read, 3)));

    // Multi-genome (3 x 10KB)
    let mut bp3 = BitPop::new(6);
    for i in 0..3 {
        bp3.add_genome(&format!("g{i}"), &pseudo_random_sequence(10_000));
    }
    bp3.build();
    group.bench_function("20bp_3genomes", |b| b.iter(|| bp3.map_read(read, 3)));

    // 1MB genome with longer read
    let mut bp_big = BitPop::new(8);
    bp_big.add_genome("large", &pseudo_random_sequence(1_000_000));
    bp_big.build();
    let read_long = &pseudo_random_sequence(100);
    group.bench_function("100bp_1mb", |b| b.iter(|| bp_big.map_read(read_long, 3)));

    // Multi-genome realistic (5 x 200KB)
    let mut bp_multi = BitPop::new(8);
    for i in 0..5 {
        bp_multi.add_genome(&format!("g{i}"), &pseudo_random_sequence(200_000));
    }
    bp_multi.build();
    group.bench_function("100bp_5genomes_1mb_total", |b| b.iter(|| bp_multi.map_read(read_long, 3)));

    group.finish();
}

// --- Serialize/Deserialize benches (NEW) ---

fn bench_serialize(c: &mut Criterion) {
    let mut group = c.benchmark_group("serialize");

    let sizes = [
        ("10kb", 10_000),
        ("100kb", 100_000),
    ];

    for (label, sz) in &sizes {
        let mut bp = BitPop::new(8);
        bp.add_genome("test", &pseudo_random_sequence(*sz));
        bp.build();
        let bp_ref = bp;

        group.throughput(Throughput::Elements(*sz as u64));
        group.bench_function(format!("{}_genome", label), |b| {
            b.iter(|| serialize_bitpop(&bp_ref))
        });
    }

    group.finish();
}

fn bench_deserialize(c: &mut Criterion) {
    let mut group = c.benchmark_group("deserialize");

    let sizes = [
        ("10kb", 10_000),
        ("100kb", 100_000),
    ];

    for (label, sz) in &sizes {
        let mut bp = BitPop::new(8);
        bp.add_genome("test", &pseudo_random_sequence(*sz));
        bp.build();
        let bytes = serialize_bitpop(&bp).unwrap();

        group.throughput(Throughput::Elements(*sz as u64));
        group.bench_function(format!("{}_genome", label), |b| {
            b.iter(|| deserialize_bitpop(&bytes))
        });
    }

    group.finish();
}

fn bench_serialize_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("serialize_roundtrip");

    let mut bp = BitPop::new(8);
    bp.add_genome("test", &pseudo_random_sequence(100_000));
    bp.build();
    let bytes = serialize_bitpop(&bp).unwrap();

    group.bench_function("100kb_roundtrip", |b| {
        b.iter(|| {
            let loaded = deserialize_bitpop(&bytes).unwrap();
            serialize_bitpop(&loaded)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_two_bit_align,
    bench_smith_waterman,
    bench_bit_parallel_exact,
    bench_comparison,
    bench_myers_vs_sw,
    bench_myers_accuracy,
    bench_fm_build,
    bench_fm_backward_search,
    bench_fm_count_occurrences,
    bench_kmer_filter,
    bench_kmer_filter_threshold,
    bench_anchor_filter,
    bench_anchor_filter_threshold,
    bench_full_pipeline,
    bench_serialize,
    bench_deserialize,
    bench_serialize_roundtrip,
);
criterion_main!(benches);
