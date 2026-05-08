use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle};

use bit_pop::{BitPop, AlignMode};
use bit_pop::fastq::{parse_reads, ReadsFormat};
use bit_pop::ncbi::{NcbiClient, NcbiConfig};
use bit_pop::cache::CacheManager;

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
    /// One-command workflow: build index (if needed) + map reads
    Run(RunArgs),
    /// Build FM-Index from FASTA file(s)
    Build(BuildArgs),
    /// Map reads to indexed genomes
    Map(MapArgs),
    /// Add genomes to existing index (incremental)
    Load(LoadArgs),
    /// Show index statistics
    Stats(StatsArgs),
    /// Search NCBI for genome accessions by organism name
    Search(SearchArgs),
    /// Fetch genome sequences from NCBI and build index
    Fetch(FetchArgs),
    /// Update cached genomes with latest versions from NCBI
    Update(UpdateArgs),
}

#[derive(clap::Args)]
struct RunArgs {
    /// Genome source: FASTA file, folder of FASTA files, or organism name with --ncbi
    genome: Option<String>,

    /// Use existing .bitpop index file (instead of building from genomes)
    #[arg(short, long)]
    index: Option<PathBuf>,

    /// Reads file (FASTQ/FASTA) for single-end mode
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// Fetch genome from NCBI instead of using local file
    #[arg(short, long)]
    ncbi: bool,

    /// Output SAM file (default: <reads_name>.sam)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// K-mer size (default: 10)
    #[arg(short, long, default_value = "10")]
    k: usize,

    /// Auto-scale k-mer size based on genome size
    #[arg(long)]
    auto_k: bool,

    /// Use spaced seed pattern for k-mer matching (better for high-error long reads)
    #[arg(short = 's', long)]
    spaced_seed: bool,

    /// Read type: short (Illumina, k=10) or long (Nanopore/PacBio, auto k)
    #[arg(long, default_value = "short")]
    read_type: String,

    /// Alignment mode: xor (fast), sw (accurate), hybrid (balanced)
    #[arg(short, long, default_value = "hybrid")]
    align_mode: String,

    /// Minimum alignment score (0.0-1.0)
    #[arg(short, long, default_value = "0.7")]
    min_score: f64,

    /// Minimum average quality score for FASTQ reads
    #[arg(short = 'q', long, default_value = "0")]
    min_quality: u8,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// NCBI API key
    #[arg(long)]
    api_key: Option<String>,

    /// Email for NCBI request tracking
    #[arg(long)]
    email: Option<String>,

    /// Force rebuild index even if cached
    #[arg(long)]
    force: bool,

    /// Number of top rarest k-mers to try as anchors (default: 1)
    #[arg(long, default_value = "1")]
    top_n: usize,

    /// Use memory-mapped I/O for FASTA file loading (reduces memory usage)
    #[cfg(feature = "mmap")]
    #[arg(long)]
    mmap: bool,
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

    /// Auto-scale k-mer size based on genome size
    #[arg(long)]
    auto_k: bool,

    /// Read type: short (Illumina, k=10) or long (Nanopore/PacBio, auto k)
    #[arg(long, default_value = "short")]
    read_type: String,

    /// Number of threads
    #[arg(short, long, default_value = "1")]
    threads: usize,

    /// Use memory-mapped I/O for FASTA file loading (reduces memory usage)
    #[cfg(feature = "mmap")]
    #[arg(long)]
    mmap: bool,
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

    /// Use memory-mapped I/O for FASTA file loading (reduces memory usage)
    #[cfg(feature = "mmap")]
    #[arg(long)]
    mmap: bool,

    /// Number of threads for parallel index build
    #[arg(short, long, default_value = "1")]
    threads: usize,
}

#[derive(clap::Args)]
struct StatsArgs {
    /// Index path
    #[arg(short, long, required = true)]
    index: PathBuf,
}

#[derive(clap::Args)]
struct SearchArgs {
    /// Organism name to search (e.g., "Escherichia coli")
    #[arg(short, long, required = true)]
    organism: String,

    /// Filter by molecule type (e.g., "genomic DNA")
    #[arg(short, long, default_value = "genomic DNA")]
    molecule_type: String,

    /// Maximum number of results to return
    #[arg(short = 'n', long, default_value = "10")]
    max_results: usize,

    /// NCBI API key for higher rate limit
    #[arg(long)]
    api_key: Option<String>,

    /// Email for NCBI request tracking
    #[arg(long)]
    email: Option<String>,
}

#[derive(clap::Args)]
struct FetchArgs {
    /// Accession ID(s) to fetch (e.g., NC_000913.3)
    #[arg(short, long, required = true)]
    accession: Vec<String>,

    /// Output index path
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// K-mer size
    #[arg(short, long, default_value = "10")]
    k: usize,

    /// Auto-scale k-mer size based on genome size
    #[arg(long)]
    auto_k: bool,

    /// Output FASTA file instead of building index
    #[arg(short, long)]
    fasta_only: bool,

    /// Number of threads for index build
    #[arg(short, long, default_value = "1")]
    threads: usize,

    /// NCBI API key for higher rate limit
    #[arg(long)]
    api_key: Option<String>,

    /// Email for NCBI request tracking
    #[arg(long)]
    email: Option<String>,

    /// Custom cache directory
    #[arg(long)]
    cache_dir: Option<PathBuf>,

    /// Force re-download even if cached
    #[arg(long)]
    force: bool,
}

#[derive(clap::Args)]
struct UpdateArgs {
    /// Index path to check for updates
    #[arg(short, long)]
    index: Option<PathBuf>,

    /// NCBI API key for higher rate limit
    #[arg(long)]
    api_key: Option<String>,

    /// Email for NCBI request tracking
    #[arg(long)]
    email: Option<String>,

    /// Custom cache directory
    #[arg(long)]
    cache_dir: Option<PathBuf>,

    /// Force update all genomes
    #[arg(long)]
    force: bool,
}

#[tokio::main]
async fn main() {
    let cli = Cli::parse();

    if let Err(e) = match cli.command {
        Commands::Run(args) => cmd_run(&args).await,
        Commands::Build(args) => { cmd_build(&args, cli.verbose); Ok(()) }
        Commands::Map(args) => { cmd_map(&args, cli.verbose); Ok(()) }
        Commands::Load(args) => { cmd_load(&args, cli.verbose); Ok(()) }
        Commands::Stats(args) => { cmd_stats(&args, cli.verbose); Ok(()) }
        Commands::Search(args) => { cmd_search(&args, cli.verbose).await; Ok(()) }
        Commands::Fetch(args) => cmd_fetch(&args, cli.verbose).await,
        Commands::Update(args) => cmd_update(&args, cli.verbose).await,
    } {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn cmd_build(args: &BuildArgs, verbose: bool) {
    let start = Instant::now();

    println!("Building FM-Index...");

    let mut bp = BitPop::new(args.k);
    bp.set_auto_k(args.auto_k);
    bp.set_read_type(&args.read_type);
    let mut total_bases: usize = 0;

    let pb = ProgressBar::new(args.fasta.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner} {msg}: [{elapsed_precise} {bar:40} {pos}/{len}]")
        .unwrap());

    for fasta_path in &args.fasta {
        let path_str = fasta_path.to_string_lossy().to_string();
        pb.set_message(format!("Loading: {}", path_str));

        #[cfg(feature = "mmap")]
        let ids = if args.mmap {
            bp.load_genome_fasta_mmap(&path_str)
        } else {
            bp.load_genome_fasta(&path_str)
        };
        #[cfg(not(feature = "mmap"))]
        let ids = bp.load_genome_fasta(&path_str);

        match ids {
            Ok(gids) => {
                for gid in gids {
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
                pb.finish_with_message(format!("Error loading {}: {}", path_str, e));
                eprintln!("Error loading {}: {}", path_str, e);
                std::process::exit(1);
            }
        }
        pb.inc(1);
    }
    pb.finish_with_message("Genomes loaded");

    if verbose {
        println!("Building index...");
    }
    let build_start = Instant::now();
    if args.threads > 1 {
        if verbose {
            println!("  Using {} threads for parallel build", args.threads);
        }
        bp.build_parallel();
    } else {
        bp.build();
    }
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
                let total = reads.len();

                let pb = ProgressBar::new(total as u64);
                pb.set_style(ProgressStyle::default_bar()
                    .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                    .unwrap());

                let mapped: Vec<(String, String, Vec<bit_pop::QualityMappingResult>)> = reads.iter()
                    .enumerate()
                    .map(|(i, (name, seq, qual))| {
                        let results = bp.map_read_with_quality_mode(seq, qual, align_mode, args.min_quality, 50);
                        if (i + 1) % 10 == 0 || i + 1 == total {
                            pb.set_position((i + 1) as u64);
                            pb.set_message(format!("{}/{} reads", i + 1, total));
                        }
                        (name.clone(), seq.clone(), results)
                    })
                    .collect();

                pb.finish_with_message("Mapping complete");

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
                let total = reads_refs.len();

                if args.reads_threads > 1 {
                    let pb = ProgressBar::new(total as u64);
                    pb.set_style(ProgressStyle::default_bar()
                        .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                        .unwrap());
                    let pb_clone = pb.clone();

                    let result = bp.map_reads_parallel_with_progress(
                        &reads_refs,
                        args.output.to_str().unwrap(),
                        50,
                        if total > 1000 { 100 } else { 10 },
                        move |completed, total| {
                            pb_clone.set_position(completed as u64);
                            pb_clone.set_message(format!("{}/{} reads", completed, total));
                        },
                    ).unwrap_or(0);

                    pb.finish_with_message("Mapping complete");
                    result
                } else {
                    let pb = ProgressBar::new(total as u64);
                    pb.set_style(ProgressStyle::default_bar()
                        .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                        .unwrap());
                    let pb_clone = pb.clone();

                    let result = bp.map_reads_to_sam_with_progress(
                        &reads_refs,
                        args.output.to_str().unwrap(),
                        50,
                        if total > 1000 { 100 } else { 10 },
                        move |completed, total| {
                            pb_clone.set_position(completed as u64);
                            pb_clone.set_message(format!("{}/{} reads", completed, total));
                        },
                    ).unwrap_or(0);

                    pb.finish_with_message("Mapping complete");
                    result
                }
            }
        }
    } else if args.reads_threads > 1 {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        let total = reads_refs.len();

        let pb = ProgressBar::new(total as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
            .unwrap());
        let pb_clone = pb.clone();

        let result = bp.map_reads_parallel_with_progress(
            &reads_refs,
            args.output.to_str().unwrap(),
            50,
            if total > 1000 { 100 } else { 10 },
            move |completed, total| {
                pb_clone.set_position(completed as u64);
                pb_clone.set_message(format!("{}/{} reads", completed, total));
            },
        ).unwrap_or(0);

        pb.finish_with_message("Mapping complete");
        result
    } else {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        let total = reads_refs.len();

        let pb = ProgressBar::new(total as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
            .unwrap());
        let pb_clone = pb.clone();

        let result = bp.map_reads_to_sam_with_progress(
            &reads_refs,
            args.output.to_str().unwrap(),
            50,
            if total > 1000 { 100 } else { 10 },
            move |completed, total| {
                pb_clone.set_position(completed as u64);
                pb_clone.set_message(format!("{}/{} reads", completed, total));
            },
        ).unwrap_or(0);

        pb.finish_with_message("Mapping complete");
        result
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

    let total_pairs = pairs.len();
    let pb = ProgressBar::new(total_pairs as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner} Mapping pairs: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
        .unwrap());

    let mapped_count = if min_quality > 0 {
        let result = bp.map_paired_reads_parallel_quality(
            &pairs,
            output.to_str().unwrap(),
            min_quality,
            50,
        ).unwrap_or(0);

        pb.set_position(total_pairs as u64);
        pb.set_message(format!("{} pairs", total_pairs));
        pb.finish_with_message("Mapping complete");
        result
    } else {
        let result = bp.map_paired_reads_parallel(
            &pairs,
            output.to_str().unwrap(),
            50,
        ).unwrap_or(0);

        pb.set_position(total_pairs as u64);
        pb.set_message(format!("{} pairs", total_pairs));
        pb.finish_with_message("Mapping complete");
        result
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
        #[cfg(feature = "mmap")]
        let ids = if args.mmap {
            bp.load_genome_fasta_mmap(&path_str)
        } else {
            bp.load_genome_fasta(&path_str)
        };
        #[cfg(not(feature = "mmap"))]
        let ids = bp.load_genome_fasta(&path_str);

        match ids {
            Ok(gids) => {
                for gid in gids {
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
    if args.threads > 1 {
        if verbose {
            println!("  Using {} threads for parallel build", args.threads);
        }
        bp.build_parallel();
    } else {
        bp.build();
    }

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

async fn cmd_search(args: &SearchArgs, _verbose: bool) {
    let _start = Instant::now();

    let mut config = NcbiConfig::new();
    if let Some(ref key) = args.api_key {
        config = config.with_api_key(key.clone());
    }
    if let Some(ref email) = args.email {
        config = config.with_email(email.clone());
    }

    let mut client = NcbiClient::new(config);

    println!("Searching NCBI for: {} ({})", args.organism, args.molecule_type);

    let search_start = Instant::now();
    let result = match client.search(&format!("{}[Organism]", args.organism)).await {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Search failed: {}", e);
            std::process::exit(1);
        }
    };
    println!("Search completed in {:.2}s", search_start.elapsed().as_secs_f64());
    println!("Found {} results", result.count);

    if result.idlist.is_empty() {
        println!("No results found. Try a different organism name or molecule type.");
        return;
    }

    let display_count = result.idlist.len().min(args.max_results);
    println!("\nTop {} results:", display_count);
    println!("{:<25} {:<50} {:<10} Type", "Accession", "Description", "Length");
    println!("{:-<100}", "");

    if result.idlist.len() > display_count {
        println!("  ... and {} more results (use -n to increase)", result.idlist.len() - display_count);
    }

    // Fetch summaries for all IDs
    let batch_size = 200;
    let mut all_docsums: Vec<bit_pop::ncbi::DocSum> = Vec::new();

    for chunk in result.idlist.chunks(batch_size) {
        let ids: Vec<&str> = chunk.iter().map(|s| s.as_str()).collect();
        match client.summary(&ids).await {
            Ok(docsums) => all_docsums.extend(docsums),
            Err(e) => {
                eprintln!("  Warning: failed to fetch summaries: {}", e);
                break;
            }
        }
    }

    // Filter for RefSeq genomic sequences
    let filtered: Vec<&bit_pop::ncbi::DocSum> = all_docsums.iter()
        .filter(|ds| {
            let is_refseq = ds.title.as_ref().map(|t| t.contains("RefSeq")).unwrap_or(false);
            let is_genomic = ds.nuc_genesim.as_ref().map(|n| n.contains("Genomic DNA")).unwrap_or(false);
            is_refseq || is_genomic
        })
        .take(display_count)
        .collect();

    if filtered.is_empty() {
        // Fall back to showing all results
        for ds in &all_docsums[..display_count.min(all_docsums.len())] {
            let accession = ds.id.clone();
            let title = ds.title.as_deref().unwrap_or("N/A");
            let pavg = ds.pavg.as_deref().unwrap_or("?");
            let title_display = title.chars().take(50).collect::<String>();
            println!("{:<25} {:<50} {:<10} -", accession, title_display, pavg);
        }
    } else {
        for ds in &filtered {
            let accession = ds.id.clone();
            let title = ds.title.as_deref().unwrap_or("N/A");
            let pavg = ds.pavg.as_deref().unwrap_or("?");
            let title_display = title.chars().take(50).collect::<String>();
            println!("{:<25} {:<50} {:<10} RefSeq", accession, title_display, pavg);
        }
    }
}

async fn cmd_fetch(args: &FetchArgs, _verbose: bool) -> Result<(), String> {
    let start = Instant::now();

    let mut config = NcbiConfig::new();
    if let Some(ref key) = args.api_key {
        config = config.with_api_key(key.clone());
    }
    if let Some(ref email) = args.email {
        config = config.with_email(email.clone());
    }

    let mut client = NcbiClient::new(config);
    let mut cache = CacheManager::new(args.cache_dir.clone()).map_err(|e| e.to_string())?;

    println!("Fetching {} genome(s) from NCBI...", args.accession.len());

    let mut genomes: Vec<(String, String)> = Vec::new();
    let mut failed = Vec::new();

    for accession in &args.accession {
        if args.force && cache.has_sequence(accession) {
            if let Err(e) = cache.remove_genome(accession) {
                eprintln!("  Warning: failed to remove cached {}: {}", accession, e);
            }
        }

        let fasta = if cache.has_sequence(accession) {
            None
        } else {
            match client.fetch_by_accession_version(accession).await {
                Ok(f) => Some(f),
                Err(e) => {
                    eprintln!("  Error fetching {}: {}", accession, e);
                    failed.push(accession.clone());
                    continue;
                }
            }
        };

        if let Some(f) = fasta {
            let parts: Vec<&str> = accession.split('.').collect();
            let version = if parts.len() >= 2 { parts[1] } else { "1" };
            let base = if parts.len() >= 2 { parts[0] } else { accession };
            cache.cache_sequence(accession, version, base, &f)
                .map_err(|e| e.to_string())?;
        }

        let result = if cache.has_sequence(accession) {
            let _genome = cache.manifest().get(accession).unwrap();
            let path = cache.get_fasta_path(accession);
            Some(path)
        } else {
            None
        };

        match result {
            Some(path) => {
                let content = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
                let genome = cache.manifest().get(accession).unwrap();
                genomes.push((genome.accession.clone(), content));
                println!("  Fetched: {} ({} bases)", accession, genome.genome_size);
            }
            None => {
                failed.push(accession.clone());
            }
        }
    }

    if genomes.is_empty() {
        eprintln!("No genomes were successfully fetched.");
        if !failed.is_empty() {
            eprintln!("Failed accessions: {}", failed.join(", "));
        }
        std::process::exit(1);
    }

    if !args.fasta_only {
        println!("\nBuilding index...");
        let build_start = Instant::now();

        let mut bp = BitPop::new(args.k);
        bp.set_auto_k(args.auto_k);
        for (name, seq) in &genomes {
            bp.add_genome(name, seq);
        }
        bp.build();

        let build_time = build_start.elapsed();

        println!("Saving index to {}...", args.output.to_string_lossy());
        match bp.serialize_to_file(args.output.to_str().unwrap()) {
            Ok(_) => {
                let file_size = std::fs::metadata(&args.output).map(|m| m.len()).unwrap_or(0);
                for (name, _) in &genomes {
                    if let Some(_genome) = cache.manifest().get(name) {
                        let _ = cache.cache_index(name, &args.output, args.k);
                    }
                }

                println!("\nDone!");
                println!("  Genomes:    {}", genomes.len());
                println!("  Index size: {} bytes ({:.1} MB)", file_size, file_size as f64 / 1_000_000.0);
                println!("  Build time: {:.2}s", build_time.as_secs_f64());
                println!("  Total time: {:.2}s", start.elapsed().as_secs_f64());

                if !failed.is_empty() {
                    println!("\n  Failed: {}", failed.join(", "));
                }
            }
            Err(e) => {
                eprintln!("Error saving index: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        println!("\nFASTQ-only mode: {} genomes cached", genomes.len());
        if !failed.is_empty() {
            println!("  Failed: {}", failed.join(", "));
        }
        println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
    }

    Ok(())
}

async fn cmd_update(args: &UpdateArgs, _verbose: bool) -> Result<(), String> {
    let start = Instant::now();

    let mut config = NcbiConfig::new();
    if let Some(ref key) = args.api_key {
        config = config.with_api_key(key.clone());
    }
    if let Some(ref email) = args.email {
        config = config.with_email(email.clone());
    }

    let mut client = NcbiClient::new(config);
    let mut cache = CacheManager::new(args.cache_dir.clone()).unwrap_or_else(|e| {
        eprintln!("Failed to initialize cache: {}", e);
        std::process::exit(1);
    });

    println!("Checking for updates in {} genome(s)...", cache.manifest().len());

    if cache.manifest().is_empty() {
        println!("No genomes cached. Use 'fetch' to download genomes first.");
        return Ok(());
    }

    let mut updated = Vec::new();
    let mut already_current = Vec::new();

    let genomes_list: Vec<(String, String, String, String)> = cache.list_genomes()
        .iter()
        .map(|g| (g.accession.clone(), g.version.clone(), g.base_accession.clone(), g.checksum.clone()))
        .collect();

    for (acc, version, base_accession, checksum) in genomes_list {
        print!("  Checking {}... ", acc);
        match client.fetch_by_accession_version(&acc).await {
            Ok(fasta) => {
                use sha2::{Sha256, Digest};
                let mut hasher = Sha256::new();
                hasher.update(fasta.as_bytes());
                let new_checksum = format!("{:x}", hasher.finalize());

                if args.force || checksum != new_checksum {
                    let _ = cache.cache_sequence(&acc, &version, &base_accession, &fasta);
                    println!("UPDATED");
                    updated.push(acc.clone());
                } else {
                    println!("up to date");
                    already_current.push(acc.clone());
                }
            }
            Err(e) => {
                println!("error: {}", e);
            }
        }
    }

    println!("\nUpdate complete:");
    println!("  Updated:   {}", updated.len());
    println!("  Up to date: {}", already_current.len());
    println!("  Total time: {:.2}s", start.elapsed().as_secs_f64());

    if !updated.is_empty() {
        println!("\nUpdated genomes:");
        for acc in &updated {
            println!("  - {}", acc);
        }
    }

    Ok(())
}

fn sha256_file(path: &Path) -> Result<String, String> {
    use sha2::{Sha256, Digest};
    let data = std::fs::read(path).map_err(|e| e.to_string())?;
    let mut hasher = Sha256::new();
    hasher.update(&data);
    Ok(format!("{:x}", hasher.finalize()))
}

fn default_output_path(reads_path: &Path) -> PathBuf {
    let stem = reads_path.file_stem().unwrap_or_default();
    let parent = reads_path.parent().unwrap_or_else(|| Path::new("."));
    let mut name = stem.to_string_lossy().to_string();
    if name.ends_with(".fastq") || name.ends_with(".fasta") {
        name = name.trim_end_matches(".fastq").trim_end_matches(".fasta").to_string();
    }
    parent.join(format!("{}.sam", name))
}

fn find_or_build_index(genome_paths: &[PathBuf], k: usize, auto_k: bool, force: bool) -> Result<BitPop, String> {
    if genome_paths.is_empty() {
        return Err("No genome files provided".to_string());
    }

    if genome_paths.len() == 1 {
        let genome_path = &genome_paths[0];
        let index_path = genome_path.with_extension("bitpop");

        if !force && index_path.exists() {
            let _genome_hash = sha256_file(genome_path)?;
            let meta = std::fs::metadata(&index_path).map_err(|e| e.to_string())?;
            let index_mtime = meta.modified().map_err(|e| e.to_string())?;
            let genome_mtime = std::fs::metadata(genome_path).map_err(|e| e.to_string())?.modified().map_err(|e| e.to_string())?;

            if genome_mtime <= index_mtime {
                println!("  Using cached index: {}", index_path.display());
                match BitPop::deserialize_from_file(index_path.to_str().unwrap()) {
                    Ok(bp) => {
                        if bp.genome_count() > 0 {
                            return Ok(bp);
                        }
                    }
                    Err(_) => {
                        println!("  Cache corrupted, rebuilding...");
                    }
                }
            }
        }
    }

    println!("  Building index ({} genomes, k={})...", genome_paths.len(), k);
    let build_start = Instant::now();

    let mut bp = BitPop::new(k);
    bp.set_auto_k(auto_k);
    for path in genome_paths {
        let path_str = path.to_string_lossy();
        let ids = bp.load_genome_fasta(&path_str)
            .map_err(|e| format!("Failed to load {}: {}", path.display(), e))?;
        if let Some(name) = ids.first().and_then(|&gid| bp.genome_name(gid)) {
            let seq_len = bp.genome_seq_len(ids[0]).unwrap_or(0);
            println!("    Loaded: {} ({} bases)", name, seq_len);
        }
    }

    bp.build();
    let build_time = build_start.elapsed();
    println!("  Index built in {:.2}s", build_time.as_secs_f64());

    if genome_paths.len() == 1 {
        let index_path = genome_paths[0].with_extension("bitpop");
        if let Err(e) = bp.serialize_to_file(index_path.to_str().unwrap()) {
            eprintln!("  Warning: failed to cache index: {}", e);
        }
    }

    Ok(bp)
}

#[cfg(feature = "mmap")]
fn find_or_build_index_mmap(genome_paths: &[PathBuf], k: usize, auto_k: bool, force: bool) -> Result<BitPop, String> {
    if genome_paths.is_empty() {
        return Err("No genome files provided".to_string());
    }

    if genome_paths.len() == 1 {
        let genome_path = &genome_paths[0];
        let index_path = genome_path.with_extension("bitpop");

        if !force && index_path.exists() {
            let genome_hash = sha256_file(genome_path)?;
            let meta = std::fs::metadata(&index_path).map_err(|e| e.to_string())?;
            let index_mtime = meta.modified().map_err(|e| e.to_string())?;
            let genome_mtime = std::fs::metadata(genome_path).map_err(|e| e.to_string())?.modified().map_err(|e| e.to_string())?;

            if genome_mtime <= index_mtime {
                println!("  Using cached index: {}", index_path.display());
                match BitPop::deserialize_from_file(index_path.to_str().unwrap()) {
                    Ok(bp) => {
                        if bp.genome_count() > 0 {
                            return Ok(bp);
                        }
                    }
                    Err(_) => {
                        println!("  Cache corrupted, rebuilding...");
                    }
                }
            }
        }
    }

    println!("  Building index (mmap, {} genomes, k={})...", genome_paths.len(), k);
    let build_start = Instant::now();

    let mut bp = BitPop::new(k);
    bp.set_auto_k(auto_k);
    for path in genome_paths {
        let path_str = path.to_string_lossy();
        let ids = bp.load_genome_fasta_mmap(&path_str)
            .map_err(|e| format!("Failed to load {}: {}", path.display(), e))?;
        if let Some(name) = ids.first().and_then(|&gid| bp.genome_name(gid)) {
            let seq_len = bp.genome_seq_len(ids[0]).unwrap_or(0);
            println!("    Loaded: {} ({} bases)", name, seq_len);
        }
    }

    bp.build();
    let build_time = build_start.elapsed();
    println!("  Index built in {:.2}s", build_time.as_secs_f64());

    if genome_paths.len() == 1 {
        let index_path = genome_paths[0].with_extension("bitpop");
        if let Err(e) = bp.serialize_to_file(index_path.to_str().unwrap()) {
            eprintln!("  Warning: failed to cache index: {}", e);
        }
    }

    Ok(bp)
}

async fn cmd_run(args: &RunArgs) -> Result<(), String> {
    let start = Instant::now();
    println!("Bit-Pop run");
    println!("═══════════");

    let use_index = args.index.is_some();
    let total_steps = if use_index { 2 } else { 3 };

    // Validate: --index and --ncbi are mutually exclusive
    if use_index && args.ncbi {
        return Err("--index and --ncbi cannot be used together".to_string());
    }

    // Validate: need either --index or genome source
    if !use_index && args.genome.is_none() {
        return Err("Either --index or genome path required".to_string());
    }

    // Step 1: Resolve genome source (only if not using --index)
    let genome_paths: Vec<PathBuf> = if !use_index {
        let genome = args.genome.clone().unwrap();
        if args.ncbi {
            println!("\n[1/{}] Fetching '{}' from NCBI...", total_steps, genome);
            let mut config = NcbiConfig::new();
            if let Some(ref key) = args.api_key {
                config = config.with_api_key(key.clone());
            }
            if let Some(ref email) = args.email {
                config = config.with_email(email.clone());
            }
            let mut client = NcbiClient::new(config);
            let mut cache = CacheManager::new(None).map_err(|e| e.to_string())?;

            let accessions = if genome.starts_with("NC_") || genome.starts_with("AC_") || genome.contains('.') {
                vec![genome.clone()]
            } else {
                let search_result = client.search(&format!("{}[Organism]", genome)).await
                    .map_err(|e| format!("NCBI search failed: {}", e))?;
                if search_result.idlist.is_empty() {
                    return Err(format!("No genomes found for '{}'", genome));
                }
                vec![search_result.idlist[0].clone()]
            };

            let mut paths = Vec::new();
            for acc in &accessions {
                print!("  Fetching {}... ", acc);
                let _fasta = if !args.force && cache.has_sequence(acc) {
                    println!("(cached)");
                    None
                } else {
                    let f = client.fetch_by_accession_version(acc).await
                        .map_err(|e| format!("Failed to fetch {}: {}", acc, e))?;
                    let parts: Vec<&str> = acc.split('.').collect();
                    let version = if parts.len() >= 2 { parts[1] } else { "1" };
                    let base = if parts.len() >= 2 { parts[0] } else { acc };
                    cache.cache_sequence(acc, version, base, &f)
                        .map_err(|e| e.to_string())?;
                    println!("({} bases)", f.lines().filter(|l| !l.starts_with('>')).map(|l| l.len()).sum::<usize>() / 2);
                    Some(f)
                };
                let path = cache.get_fasta_path(acc);
                paths.push(path);
            }
            paths
        } else {
            println!("\n[1/{}] Resolving genome source...", total_steps);
            let path = PathBuf::from(&genome);

            if path.is_dir() {
                let entries: Vec<_> = std::fs::read_dir(&path)
                    .map_err(|e| format!("Cannot read directory {}: {}", path.display(), e))?
                    .filter_map(|e| e.ok())
                    .map(|e| e.path())
                    .filter(|p| {
                        p.extension().map(|e| e == "fna" || e == "fasta" || e == "fa").unwrap_or(false)
                    })
                    .collect();
                if entries.is_empty() {
                    return Err(format!("No .fna/.fasta files found in {}", path.display()));
                }
                println!("  Found {} genome file(s) in {}", entries.len(), path.display());
                entries
            } else if path.exists() {
                println!("  Genome: {}", path.display());
                vec![path]
            } else {
                return Err(format!("Path '{}' not found", genome));
            }
        }
    } else {
        Vec::new()
    };

    // Step 2 (or 1): Build or load index
    let step_index = if use_index { 1 } else { 2 };
    println!("\n[{}] Preparing index...", step_index);
    let mut bp = if let Some(ref index_path) = args.index {
        println!("  Loading index: {}", index_path.display());
        BitPop::deserialize_from_file(index_path.to_str().unwrap())
            .map_err(|e| format!("Failed to load index: {}", e))?
    } else {
        #[cfg(feature = "mmap")]
        if args.mmap {
            find_or_build_index_mmap(&genome_paths, args.k, args.auto_k, args.force)?
        } else {
            find_or_build_index(&genome_paths, args.k, args.auto_k, args.force)?
        }
        #[cfg(not(feature = "mmap"))]
        find_or_build_index(&genome_paths, args.k, args.auto_k, args.force)?
    };

    if args.spaced_seed {
        println!("  Spaced seed: enabled (pattern: 11101001110111)");
        bp.set_spaced_seed(true);
    }

    if args.top_n > 1 {
        bp.set_top_n(args.top_n);
    }

    bp.set_read_type(&args.read_type);
    println!("  Read type: {}", args.read_type);

    // Step 3 (or 2): Map reads
    let step_map = if use_index { 2 } else { 3 };
    println!("\n[{}] Mapping reads...", step_map);

    let _mapped_count = if let (Some(r1_path), Some(r2_path)) = (&args.reads_1, &args.reads_2) {
        // Paired-end mode
        println!("  Paired-end mode");
        println!("    R1: {}", r1_path.display());
        println!("    R2: {}", r2_path.display());

        let pairs = bit_pop::fastq::parse_paired_fastq(r1_path.to_str().unwrap(), r2_path.to_str().unwrap())
            .map_err(|e| format!("Failed to parse paired FASTQ: {}", e))?;
        println!("  Loaded {} read pairs", pairs.len());

        let output_path = if let Some(ref p) = args.output {
            p.clone()
        } else {
            default_output_path(r1_path)
        };

        let total_pairs = pairs.len();
        let pb = ProgressBar::new(total_pairs as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner} Mapping pairs: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
            .unwrap());

        let mapped = if args.min_quality > 0 {
            let result = bp.map_paired_reads_parallel_quality(
                &pairs,
                output_path.to_str().unwrap(),
                args.min_quality,
                50,
            ).map_err(|e| format!("Mapping failed: {}", e))?;

            pb.set_position(total_pairs as u64);
            pb.set_message(format!("{} pairs", total_pairs));
            pb.finish_with_message("Mapping complete");
            result
        } else {
            let result = bp.map_paired_reads_parallel(
                &pairs,
                output_path.to_str().unwrap(),
                50,
            ).map_err(|e| format!("Mapping failed: {}", e))?;

            pb.set_position(total_pairs as u64);
            pb.set_message(format!("{} pairs", total_pairs));
            pb.finish_with_message("Mapping complete");
            result
        };

        // Parse SAM for results
        let genome_counts = parse_sam_summary(&output_path);
        let total = genome_counts.values().sum::<usize>();

        println!("\n═══════════");
        println!("Done!");
        println!("  Mapped:     {}/{} pairs", mapped, pairs.len());
        println!("  Output:     {}", output_path.display());

        if mapped > 0 {
            println!("\n  Results:");
            println!("  {:<60} {:>10} {:>8}", "Genome", "Count", "Percent");
            println!("  {:─<72}", "");
            for (name, count) in &genome_counts {
                let pct = *count as f64 / total as f64 * 100.0;
                let display_name = if name.len() > 58 {
                    format!("…{}", &name[name.len() - 55..])
                } else {
                    name.clone()
                };
                println!("  {:<60} {:>10} {:>7.1}%", display_name, count, pct);
            }
            println!("  {:─<72}", "");
            println!("  {:<60} {:>10} {:>7.1}%", "Total", total, 100.0);
        }

        mapped
    } else {
        // Single-end mode
        let reads_path = args.reads.as_ref()
            .ok_or("Either --reads (-r) or --reads-1/--reads-2 required")?;

        if !reads_path.exists() {
            return Err(format!("Reads file '{}' not found", reads_path.display()));
        }

        let reads_format = parse_reads(reads_path.to_str().unwrap())
            .map_err(|e| format!("Failed to parse reads: {}", e))?;
        let format_name = match &reads_format {
            ReadsFormat::Fasta(_) => "FASTA",
            ReadsFormat::Fastq(_) => "FASTQ",
        };
        println!("  Loaded {} reads ({})", reads_format.count(), format_name);

        let output_path = if let Some(ref p) = args.output {
            p
        } else {
            &default_output_path(reads_path)
        };

        let filtered_reads_fasta: Vec<(String, String)> = if args.min_quality > 0 && reads_format.has_quality() {
            if let ReadsFormat::Fastq(reads) = &reads_format {
                let passed = bit_pop::fastq::filter_by_quality(reads, args.min_quality);
                println!("  Quality filter (min Q{}): {}/{} reads passed",
                    args.min_quality, passed.len(), reads.len());
                passed.iter()
                    .map(|&i| (reads[i].0.clone(), reads[i].1.clone()))
                    .collect()
            } else {
                println!("  Warning: quality filtering ignored for FASTA input");
                reads_format.iter_fasta().map(|(n, s)| (n.to_string(), s.to_string())).collect()
            }
        } else {
            reads_format.iter_fasta().map(|(n, s)| (n.to_string(), s.to_string())).collect()
        };

        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta.iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();

        let total_reads = reads_refs.len();
        let mapped_count = if args.threads > 1 {
            let pb = ProgressBar::new(total_reads as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap());
            let pb_clone = pb.clone();

            let result = bp.map_reads_parallel_with_progress(
                &reads_refs,
                output_path.to_str().unwrap(),
                50,
                if total_reads > 1000 { 100 } else { 10 },
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                },
            ).map_err(|e| format!("Mapping failed: {}", e))?;

            pb.finish_with_message("Mapping complete");
            result
        } else {
            let pb = ProgressBar::new(total_reads as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap());
            let pb_clone = pb.clone();

            let result = bp.map_reads_to_sam_with_progress(
                &reads_refs,
                output_path.to_str().unwrap(),
                50,
                if total_reads > 1000 { 100 } else { 10 },
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                },
            ).map_err(|e| format!("Mapping failed: {}", e))?;

            pb.finish_with_message("Mapping complete");
            result
        };

        let elapsed = start.elapsed();

        // Parse SAM and show results
        println!("\n═══════════");
        println!("Done!");
        println!("  Mapped:     {}/{} reads", mapped_count, filtered_reads_fasta.len());
        println!("  Output:     {}", output_path.display());

        if mapped_count > 0 {
            let genome_counts = parse_sam_summary(output_path);
            let total = genome_counts.values().sum::<usize>();
            println!("\n  Results:");
            println!("  {:<60} {:>10} {:>8}", "Genome", "Count", "Percent");
            println!("  {:─<72}", "");
            for (name, count) in &genome_counts {
                let pct = *count as f64 / total as f64 * 100.0;
                let display_name = if name.len() > 58 {
                    format!("…{}", &name[name.len() - 55..])
                } else {
                    name.clone()
                };
                println!("  {:<60} {:>10} {:>7.1}%", display_name, count, pct);
            }
            println!("  {:─<72}", "");
            println!("  {:<60} {:>10} {:>7.1}%", "Total", total, 100.0);
        }

        println!("  Total time: {:.2}s", elapsed.as_secs_f64());

        mapped_count
    };

    Ok(())
}

fn parse_sam_summary(path: &Path) -> std::collections::HashMap<String, usize> {
    use std::collections::HashMap;
    use std::io::BufRead;
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut seen_reads: std::collections::HashSet<String> = std::collections::HashSet::new();

    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return counts,
    };
    let reader = std::io::BufReader::new(file);

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let read_name = fields[0].to_string();
        let ref_name = fields[2];
        if ref_name == "*" {
            continue;
        }
        // Count each read only once (use first occurrence)
        if seen_reads.insert(read_name) {
            let clean_name = ref_name.trim_end();
            *counts.entry(clean_name.to_string()).or_insert(0) += 1;
        }
    }

    counts
}
