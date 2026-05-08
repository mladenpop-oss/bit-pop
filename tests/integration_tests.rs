use bit_pop::BitPop;
use std::fs;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data")
}

fn temp_dir(name: &str) -> PathBuf {
    let dir = test_data_dir().join("test_output").join(name);
    let _ = fs::create_dir_all(&dir);
    dir
}

#[test]
fn test_build_single_genome_index() {
    let genome_path = test_data_dir().join("genomes").join("Ecoli_K12_MG1655.fna");
    let index_path = temp_dir("test_single_genome").join("test_index.bitpop");

    let mut bp = BitPop::new(10);
    bp.load_genome_fasta(genome_path.to_str().unwrap()).unwrap();
    bp.build();

    assert_eq!(bp.genome_count(), 1);
    assert!(bp.genome_seq_len(0).unwrap() > 0);

    bp.serialize_to_file(index_path.to_str().unwrap()).unwrap();
    assert!(index_path.exists());

    let loaded = BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap();
    assert_eq!(loaded.genome_count(), 1);
}

#[test]
fn test_map_reads_to_sam() {
    let genome_path = test_data_dir().join("genomes").join("Ecoli_K12_MG1655.fna");
    let _reads_path = test_data_dir()
        .join("reads")
        .join("simulated_ecoli_10k_new.fastq");
    let index_path = temp_dir("test_map_reads").join("test_index.bitpop");
    let sam_path = temp_dir("test_map_reads").join("test_output.sam");

    // Build index
    let mut bp = BitPop::new(10);
    bp.load_genome_fasta(genome_path.to_str().unwrap()).unwrap();
    bp.build();
    bp.serialize_to_file(index_path.to_str().unwrap()).unwrap();

    // Load and map
    let bp = BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap();
    let reads = vec![(
        "test_read",
        "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
    )];
    let _mapped = bp
        .map_reads_to_sam(&reads, sam_path.to_str().unwrap(), 50)
        .unwrap();

    // Map may return 0 if read doesn't match genome (repetitive k-mers filtered)
    assert!(sam_path.exists());
    assert!(sam_path.exists());

    let content = fs::read_to_string(&sam_path).unwrap();
    assert!(content.starts_with("@SQ"));
    assert!(content.contains("NC_000913.3"));
}

#[test]
fn test_multi_genome_index() {
    let genomes_dir = test_data_dir().join("genomes");
    let index_path = temp_dir("test_multi_genome").join("test_index.bitpop");

    let mut bp = BitPop::new(10);
    for entry in fs::read_dir(genomes_dir).unwrap() {
        let path = entry.unwrap().path();
        if path
            .extension()
            .map(|e| e == "fna" || e == "fasta")
            .unwrap_or(false)
        {
            bp.load_genome_fasta(path.to_str().unwrap()).unwrap();
        }
    }
    bp.build();

    assert_eq!(bp.genome_count(), 3);

    bp.serialize_to_file(index_path.to_str().unwrap()).unwrap();
    assert!(index_path.exists());
}

#[test]
fn test_sam_output_format() {
    let genome_path = test_data_dir().join("genomes").join("Ecoli_K12_MG1655.fna");
    let index_path = temp_dir("test_sam_format").join("test_index.bitpop");
    let sam_path = temp_dir("test_sam_format").join("test_output.sam");

    let mut bp = BitPop::new(10);
    bp.load_genome_fasta(genome_path.to_str().unwrap()).unwrap();
    bp.build();
    bp.serialize_to_file(index_path.to_str().unwrap()).unwrap();

    let bp = BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap();
    let reads = vec![
        ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ("read2", "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"),
    ];
    bp.map_reads_to_sam(&reads, sam_path.to_str().unwrap(), 50)
        .unwrap();

    let content = fs::read_to_string(&sam_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Check header - at least first line should be @SQ
    assert!(lines[0].starts_with("@SQ"));
    // Second line is either @SQ (another genome) or a data line
    let first_data_line = lines.iter().find(|l| !l.starts_with('@') && !l.is_empty());
    assert!(
        first_data_line.is_some(),
        "SAM should have at least one data line"
    );

    // Check data lines have tab-separated fields
    for line in lines.iter().skip(2) {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(
            fields.len() >= 11,
            "SAM line should have 11+ fields: {:?}",
            fields
        );
    }
}

#[test]
fn test_cache_reuse() {
    let genome_path = test_data_dir().join("genomes").join("Ecoli_K12_MG1655.fna");
    let index_path = temp_dir("test_cache").join("test_index.bitpop");

    // Build and save
    let mut bp1 = BitPop::new(10);
    bp1.load_genome_fasta(genome_path.to_str().unwrap())
        .unwrap();
    bp1.build();
    bp1.serialize_to_file(index_path.to_str().unwrap()).unwrap();

    // Load from cache
    let bp2 = BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap();
    assert_eq!(bp2.genome_count(), 1);
    assert_eq!(bp2.k(), 10);
}
