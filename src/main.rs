use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};
use std::time::Instant;

use bit_pop::{BitPop, AlignMode};
use bit_pop::fastq::{parse_reads, ReadsFormat};

#[derive(Parser)]
#[command(name = "bit-pop", about = "Multi-genome DNA read mapper", long_about = None)]
struct Cli {
    /// Verbose output
    #[arg(short, long)]
    verbose: bool,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build FM-Index from FASTA file(s)
    Build(BuildArgs),
    /// Map reads to indexed genomes
    Map(MapArgs),
    /// Add genomes to existing index (incremental)
    Load(LoadArgs),
    /// Show index statistics
    Stats(StatsArgs),
}

#[derive(clap::Args)]
struct BuildArgs {
    /// Input FASTA file(s)
    #[arg(short, long, required = true)]
    fasta: Vec<PathBuf>,

    /// Output index path
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// K-mer size (default: 8)
    #[arg(short, long, default_value = "8")]
    k: usize,

    /// Number of threads
    #[arg(short, long, default_value = "1")]
    threads: usize,
}

#[derive(clap::Args)]
struct MapArgs {
    /// Input index path
    #[arg(short, long, required = true)]
    index: PathBuf,

    /// Input reads file (FASTA or FASTQ). Use -1 and -2 for paired-end.
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// Output SAM file path
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// Minimum alignment score (0.0-1.0)
    #[arg(short, long, default_value = "0.7")]
    min_score: f64,

    /// Alignment mode: xor (fast), sw (accurate), hybrid (balanced)
    #[arg(short, long, default_value = "xor")]
    align_mode: String,

    /// Minimum average quality score for FASTQ reads (default: 0 = no filter)
    #[arg(short = 'q', long, default_value = "0")]
    min_quality: u8,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    reads_threads: usize,

    /// Number of top rarest k-mers to try as anchors (default: 1)
    #[arg(long, default_value = "1")]
    top_n: usize,
}

#[derive(clap::Args)]
struct LoadArgs {
    /// Existing index path
    #[arg(short, long, required = true)]
    index: PathBuf,

    /// New FASTA file(s) to add
    #[arg(short, long, required = true)]
    fasta: Vec<PathBuf>,

    /// Updated index output path
    #[arg(short, long, required = true)]
    output: PathBuf,
}

#[derive(clap::Args)]
struct StatsArgs {
    /// Index path
    #[arg(short, long, required = true)]
    index: PathBuf,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Build(args) => cmd_build(&args, cli.verbose),
        Commands::Map(args) => cmd_map(&args, cli.verbose),
        Commands::Load(args) => cmd_load(&args, cli.verbose),
        Commands::Stats(args) => cmd_stats(&args, cli.verbose),
    }
}

fn cmd_build(args: &BuildArgs, verbose: bool) {
    let start = Instant::now();

    println!("Building FM-Index...");

    let mut bp = BitPop::new(args.k);
    let mut total_bases: usize = 0;

    for fasta_path in &args.fasta {
        let path_str = fasta_path.to_string_lossy().to_string();
        if verbose {
            println!("  Loading: {}", path_str);
        }

        match bp.load_genome_fasta(&path_str) {
            Ok(ids) => {
                for gid in ids {
                    let seq_len = bp.genome_seq_len(gid).unwrap_or(0);
                    total_bases += seq_len;
                    if verbose {
                        if let Some(name) = bp.genome_name(gid) {
                            println!("    Added genome: {} ({} bases)", name, seq_len);
                        }
                    }
                }
            }
            Err(e) => {
                eprintln!("Error loading {}: {}", path_str, e);
                std::process::exit(1);
            }
        }
    }

    if verbose {
        println!("Building index...");
    }
    let build_start = Instant::now();
    bp.build();
    let build_time = build_start.elapsed();

    if verbose {
        println!("Saving index...");
    }

    match bp.serialize_to_file(args.output.to_str().unwrap()) {
        Ok(_) => {
            let elapsed = start.elapsed();
            println!(
                "Index built successfully: {} genomes, {} total bases, {} bytes",
                bp.genome_count(),
                total_bases,
                std::fs::metadata(&args.output).map(|m| m.len()).unwrap_or(0),
            );
            println!("  Build time: {:.2}s", build_time.as_secs_f64());
            println!("  Total time: {:.2}s", elapsed.as_secs_f64());
        }
        Err(e) => {
            eprintln!("Error saving index: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_map(args: &MapArgs, verbose: bool) {
    let start = Instant::now();

    println!("Loading index: {}", args.index.to_string_lossy());
    let load_start = Instant::now();

    let mut bp = match BitPop::deserialize_from_file(args.index.to_str().unwrap()) {
        Ok(bp) => bp,
        Err(e) => {
            eprintln!("Error loading index: {}", e);
            std::process::exit(1);
        }
    };

    if args.top_n > 1 {
        bp.set_top_n(args.top_n);
    }

    let load_time = load_start.elapsed();

    let align_mode = match args.align_mode.as_str() {
        "sw" => AlignMode::Sw,
        "hybrid" => AlignMode::Hybrid,
        _ => AlignMode::Xor,
    };

    println!("Index loaded in {:.3}s ({})\n", load_time.as_secs_f64(), bp.genome_count());
    if verbose {
        println!("Alignment mode: {}\n", align_mode);
    }

    // Check for paired-end mode
    if let (Some(r1_path), Some(r2_path)) = (&args.reads_1, &args.reads_2) {
        cmd_map_paired(&bp, r1_path, r2_path, &args.output, args.min_quality);
        return;
    }

    // Single-end mode
    let reads_path = match &args.reads {
        Some(p) => p.to_string_lossy().to_string(),
        None => {
            eprintln!("Error: --reads (-r) required for single-end mode, or use --reads-1/-2 for paired-end");
            std::process::exit(1);
        }
    };

    let reads_format = match parse_reads(&reads_path) {
        Ok(format) => format,
        Err(e) => {
            eprintln!("Error parsing reads: {}", e);
            std::process::exit(1);
        }
    };

    if verbose {
        match &reads_format {
            ReadsFormat::Fasta(_) => println!("FASTA detected"),
            ReadsFormat::Fastq(reads) => {
                println!("FASTQ detected ({} reads with quality scores)", reads.len());
            }
        }
    }

    println!("Loaded {} reads", reads_format.count());

    let has_quality = reads_format.has_quality();

    let filtered_reads_fasta: Vec<(String, String)> = if args.min_quality > 0 {
        match &reads_format {
            ReadsFormat::Fastq(reads) => {
                let passed = bit_pop::fastq::filter_by_quality(reads, args.min_quality);
                println!("Quality filter (min Q{}): {}/{} reads passed",
                    args.min_quality, passed.len(), reads.len());
                passed.iter()
                    .map(|&i| (reads[i].0.clone(), reads[i].1.clone()))
                    .collect()
            }
            ReadsFormat::Fasta(_) => {
                println!("Warning: quality filtering ignored for FASTA input");
                reads_format.iter_fasta().map(|(n, s)| (n.to_string(), s.to_string())).collect()
            }
        }
    } else {
        reads_format.iter_fasta().map(|(n, s)| (n.to_string(), s.to_string())).collect()
    };

    let map_start = Instant::now();

    let mapped_count = if has_quality && args.min_quality > 0 {
        match &reads_format {
            ReadsFormat::Fastq(reads) => {
                let genomes_owned: Vec<(String, usize)> = (0..bp.genome_count() as u32)
                    .filter_map(|gid| {
                        bp.genome_name(gid)
                            .map(|name| (name.to_string(), bp.genome_seq_len(gid).unwrap_or(0)))
                    })
                    .collect();

                let genome_name_refs: Vec<&str> = genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
                let genome_header: Vec<(&str, usize)> = genomes_owned.iter()
                    .map(|(n, l)| (n.as_str(), *l)).collect();

                let name_refs: Vec<&str> = genome_name_refs.clone();

                let mapped: Vec<(String, String, Vec<bit_pop::QualityMappingResult>)> = reads.iter()
                    .map(|(name, seq, qual)| {
                        let results = bp.map_read_with_quality_mode(seq, qual, align_mode, args.min_quality, 50);
                        (name.clone(), seq.clone(), results)
                    })
                    .collect();

                let mut writer = bit_pop::sam::SamWriter::new(args.output.to_str().unwrap()).unwrap();
                writer.write_header(&genome_header).unwrap();

                let mut mapped_count = 0;
                for (name, seq, results) in &mapped {
                    writer.write_quality_mappings(name, seq, results, &name_refs).unwrap();
                    if !results.is_empty() {
                        mapped_count += 1;
                    }
                }

                mapped_count
            }
            ReadsFormat::Fasta(_) => {
                let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
                    .map(|(name, seq)| (name.as_str(), seq.as_str()))
                    .collect();
                if args.reads_threads > 1 {
                    bp.map_reads_parallel(&reads_refs, args.output.to_str().unwrap(), 50).unwrap_or(0)
                } else {
                    bp.map_reads_to_sam(&reads_refs, args.output.to_str().unwrap(), 50).unwrap_or(0)
                }
            }
        }
    } else if args.reads_threads > 1 {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        bp.map_reads_parallel(&reads_refs, args.output.to_str().unwrap(), 50)
            .unwrap_or(0)
    } else {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        bp.map_reads_to_sam(&reads_refs, args.output.to_str().unwrap(), 50)
            .unwrap_or(0)
    };

    let elapsed = start.elapsed();

    println!("\nMapping complete: {}/{} reads mapped", mapped_count, filtered_reads_fasta.len());
    println!("  Alignment mode: {}", align_mode);
    println!("  Load time:  {:.3}s", load_time.as_secs_f64());
    println!("  Map time:   {:.2}s", map_start.elapsed().as_secs_f64());
    println!("  Total time: {:.2}s", elapsed.as_secs_f64());
}

fn cmd_map_paired(bp: &BitPop, r1_path: &Path, r2_path: &Path, output: &Path, min_quality: u8) {
    let map_start = Instant::now();

    println!("Paired-end mapping mode");
    println!("  R1: {}", r1_path.to_string_lossy());
    println!("  R2: {}", r2_path.to_string_lossy());

    let pairs = match bit_pop::fastq::parse_paired_fastq(r1_path.to_str().unwrap(), r2_path.to_str().unwrap()) {
        Ok(pairs) => pairs,
        Err(e) => {
            eprintln!("Error parsing paired FASTQ: {}", e);
            std::process::exit(1);
        }
    };

    println!("Loaded {} read pairs", pairs.len());

    let mapped_count = if min_quality > 0 {
        bp.map_paired_reads_parallel_quality(
            &pairs,
            output.to_str().unwrap(),
            min_quality,
            50,
        ).unwrap_or(0)
    } else {
        bp.map_paired_reads_parallel(
            &pairs,
            output.to_str().unwrap(),
            50,
        ).unwrap_or(0)
    };

    let elapsed = map_start.elapsed();

    println!("\nPaired-end mapping complete: {} pairs processed", mapped_count);
    println!("  Map time:   {:.2}s", elapsed.as_secs_f64());
    println!("  Total time: {:.2}s", map_start.elapsed().as_secs_f64());
}

fn cmd_load(args: &LoadArgs, verbose: bool) {
    let start = Instant::now();

    println!("Loading existing index...");
    let mut bp = match BitPop::deserialize_from_file(args.index.to_str().unwrap()) {
        Ok(bp) => bp,
        Err(e) => {
            eprintln!("Error loading index: {}", e);
            std::process::exit(1);
        }
    };

    let old_count = bp.genome_count();

    for fasta_path in &args.fasta {
        let path_str = fasta_path.to_string_lossy().to_string();
        if verbose {
            println!("  Adding: {}", path_str);
        }
        match bp.load_genome_fasta(&path_str) {
            Ok(ids) => {
                for gid in ids {
                    if verbose {
                        if let Some(name) = bp.genome_name(gid) {
                            println!("    Added: {}", name);
                        }
                    }
                }
            }
            Err(e) => {
                eprintln!("Error loading {}: {}", path_str, e);
                std::process::exit(1);
            }
        }
    }

    let new_count = bp.genome_count();
    println!("Added {} new genomes ({} -> {})", new_count - old_count, old_count, new_count);

    if verbose {
        println!("Rebuilding index...");
    }
    bp.build();

    match bp.serialize_to_file(args.output.to_str().unwrap()) {
        Ok(_) => {
            println!("Index updated: {}", args.output.to_string_lossy());
            println!("  Total time: {:.2}s", start.elapsed().as_secs_f64());
        }
        Err(e) => {
            eprintln!("Error saving index: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_stats(args: &StatsArgs, _verbose: bool) {
    let bp = match BitPop::deserialize_from_file(args.index.to_str().unwrap()) {
        Ok(bp) => bp,
        Err(e) => {
            eprintln!("Error loading index: {}", e);
            std::process::exit(1);
        }
    };

    let total_bases: usize = (0..bp.genome_count())
        .filter_map(|i| bp.genome_seq_len(i as u32))
        .sum();

    let file_size = std::fs::metadata(&args.index).map(|m| m.len()).unwrap_or(0);

    println!("=== Bit-Pop Index Statistics ===\n");
    println!("File size:     {} bytes ({:.1} MB)", file_size, file_size as f64 / 1_000_000.0);
    println!("Genomes:       {}", bp.genome_count());
    println!("Total bases:   {}", total_bases);
    println!("K-mer size:    {}", bp.k());
    println!("BWT length:    {}", bp.bwt_len());
    println!();

    println!("Genomes:");
    let names = bp.genome_names_ordered();
    for (i, name) in names.iter().enumerate() {
        let len = bp.genome_seq_len(i as u32).unwrap_or(0);
        println!("  {}. {} ({} bases)", i + 1, name, len);
    }
}
