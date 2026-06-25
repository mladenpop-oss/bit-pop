use bit_pop::report_atomic_progress;
use bit_pop::MappingResult;
use clap::{Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use bit_pop::cache::CacheManager;
use bit_pop::chunk_consensus::MultiChunkConsensus;
use bit_pop::concon::ConCon;
use bit_pop::consensus::MultiKConsensus;
use bit_pop::fastq::{parse_reads, ReadsFormat};
use bit_pop::ncbi::{NcbiClient, NcbiConfig};
use bit_pop::{AlignMode, BitPop, FuzzyMethod};

fn extract_cami_genome_name(basename: &str) -> String {
    if basename.starts_with("evo_") {
        basename.to_string()
    } else if basename
        .chars()
        .next()
        .map(|c| c.is_ascii_digit())
        .unwrap_or(false)
    {
        let without_run = basename.split("_run").next().unwrap_or(basename);
        let parts: Vec<&str> = without_run.splitn(2, '.').collect();
        parts[0].to_string()
    } else {
        basename.to_string()
    }
}

fn extract_pacbio_genome_name(basename: &str) -> String {
    basename.to_string()
}

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

    /// Apply EM algorithm for soft-assignment classification
    Em(EmArgs),

    /// Multi-k consensus: map reads against multiple k-indexes with voting
    Consensus(ConsensusArgs),

    /// Consensus: run `bit-pop map` for each index, then combine
    ConCon(ConConArgs),

    /// Multi chunk-% consensus: same index, different chunk sizes, voting
    ChunkConsensus(ChunkConsensusArgs),

    /// Generate taxonomic classification report from SAM output
    Tax(TaxArgs),
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

    /// Custom spaced seed pattern (e.g., "11101001110111" for 10/14). Default: "11111011111111" (13/14)
    #[arg(long, requires = "spaced_seed")]
    spaced_seed_pattern: Option<String>,

    /// Fuzzy k-mer matching method: none, fuzzy-kmer, fuzzy-seed, neighborhood
    #[arg(long, default_value = "none")]
    method: String,

    /// Maximum number of mismatches for fuzzy matching (default: 1)
    #[arg(long, default_value = "1")]
    fuzzy_mismatches: usize,

    /// Read type: short (Illumina, k=10) or long (Nanopore/PacBio, auto k)
    #[arg(long, default_value = "short")]
    read_type: String,

    /// Use golden anchor selection (quality-weighted k-mer anchors for long reads)
    #[arg(long)]
    golden_anchors: bool,

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

    /// Apply EM algorithm for soft-assignment classification (improves strain resolution)
    #[arg(long)]
    em: bool,

    /// Search radius in bp (±N around anchor position, default: 5, max: 200)
    #[arg(long, default_value = "5")]
    search_radius: isize,

    /// Chunk size for PacBio long-read mapping (0 = auto-detect, 150 recommended)
    /// Reads >1000bp are split into chunks for improved mapping rate
    #[arg(long)]
    chunk_size: Option<usize>,

    /// Chunk size as percentage of read length (0.0-1.0, e.g. 0.01 = 1%).
    /// Overrides --chunk-size when set. Enables dynamic per-read chunk sizing.
    /// Clamped to [chunk_min, chunk_max] range (default: 20-500bp).
    #[arg(long)]
    chunk_pct: Option<f64>,

    /// Minimum chunk size clamp for dynamic chunking (overrides default 20bp).
    #[arg(long)]
    chunk_min: Option<usize>,

    /// Maximum chunk size clamp for dynamic chunking (overrides default 500bp).
    #[arg(long)]
    chunk_max: Option<usize>,

    /// Minimum fraction of chunks that must agree (0.0-1.0, default: 0.0 = no threshold).
    /// Example: 0.6 requires 60% of chunks to agree before accepting a mapping.
    #[arg(long)]
    chunk_vote_threshold: Option<f64>,

    /// Number of top genomes to return per read in chunk-based mode (default: 1).
    /// Use 2-3 for multi-genome uncertainty scenarios.
    #[arg(long)]
    chunk_top_n: Option<usize>,

    /// Anchor strategy for chunk-based mapping: rarest (default), golden (quality-weighted), spaced (spaced seed)
    #[arg(long, default_value = "rarest")]
    chunk_strategy: String,

    /// Score aggregation mode for chunk-based mapping: quality (default, score*score), base (raw sum, like JNI Android)
    #[arg(long, default_value = "quality")]
    score_mode: String,

    /// Minimum anchor score threshold for chunk-based mapping (default: 0.5).
    /// Use 0.0 to match JNI Android behavior (no filtering).
    #[arg(long, default_value = "0.5")]
    anchor_min_score: f64,

    /// Use legacy anchor_filter for chunks (instead of full map_read pipeline).
    /// For testing/comparison with older behavior.
    #[arg(long)]
    anchor_filter: bool,

    /// Enable SNP-aware scoring for strain resolution.
    /// Collects mismatches across all reads, builds a SNP map, and boosts
    /// scores for genomes with known strain-specific SNPs.
    #[arg(long)]
    snp_detect: bool,

    /// Minimum support count for SNP detection (default: 3).
    /// A position must have this many reads supporting the same mismatch to be considered a SNP.
    #[arg(long, default_value = "3")]
    snp_min_support: u32,

    /// Enable homopolymer fingerprint scoring for strain resolution
    #[arg(long)]
    hf: bool,

    /// Minimum run length for homopolymer fingerprint (default: 3)
    #[arg(long, default_value = "3")]
    hf_min: usize,

    /// Output BAM format instead of SAM
    #[arg(long)]
    bam: bool,
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

    /// Fuzzy k-mer matching method: none, fuzzy-kmer, fuzzy-seed, neighborhood
    #[arg(long, default_value = "none")]
    method: String,

    /// Maximum number of mismatches for fuzzy matching (default: 1)
    #[arg(long, default_value = "1")]
    fuzzy_mismatches: usize,

    /// Number of threads
    #[arg(short, long, default_value = "1")]
    threads: usize,

    /// Use memory-mapped I/O for FASTA file loading (reduces memory usage)
    #[cfg(feature = "mmap")]
    #[arg(long)]
    mmap: bool,

    /// Use CAMI mode: extract genome name from filename (e.g., 1036554.gt1kb.fasta -> 1036554)
    #[arg(long)]
    cami: bool,

    /// Use PacBio mode: extract genome name from filename (e.g., A_baumannii_AYE_bc2001.fa -> A_baumannii_AYE_bc2001)
    #[arg(long)]
    pacbio: bool,

    /// Enable spaced seed pattern for k-mer matching (builds spaced seed hash index)
    #[arg(long)]
    spaced_seed: bool,

    /// Custom spaced seed pattern (e.g., "11101001110111" for 10/14). Default: "11111011111111" (13/14)
    #[arg(long, requires = "spaced_seed")]
    spaced_seed_pattern: Option<String>,

    /// Search radius in bp (±N around anchor position, default: 5, max: 200)
    #[arg(long, default_value = "5")]
    search_radius: isize,

    /// Enable homopolymer fingerprint scoring for strain resolution
    #[arg(long)]
    hf: bool,

    /// Minimum run length for homopolymer fingerprint (default: 3)
    #[arg(long, default_value = "3")]
    hf_min: usize,
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

    /// Fuzzy k-mer matching method: none, fuzzy-kmer, fuzzy-seed, neighborhood
    #[arg(long, default_value = "none")]
    method: String,

    /// Maximum number of mismatches for fuzzy matching (default: 1)
    #[arg(long, default_value = "1")]
    fuzzy_mismatches: usize,

    /// Apply EM algorithm for soft-assignment classification (improves strain resolution)
    #[arg(long)]
    em: bool,

    /// Search radius in bp (±N around anchor position, default: 5, max: 200)
    #[arg(long, default_value = "5")]
    search_radius: isize,

    /// Enable spaced seed matching (uses spaced seed hash index if available)
    #[arg(long)]
    spaced_seed: bool,

    /// Custom spaced seed pattern (e.g., "11101001110111" for 10/14). Default: "11111011111111" (13/14)
    #[arg(long, requires = "spaced_seed")]
    spaced_seed_pattern: Option<String>,

    /// Use golden anchor selection (quality-weighted k-mer anchors for long reads)
    #[arg(long)]
    golden_anchors: bool,

    /// Chunk size for PacBio long-read mapping (0 = auto-detect, 150 recommended)
    #[arg(long)]
    chunk_size: Option<usize>,

    /// Chunk size as percentage of read length (0.0-1.0, e.g. 0.01 = 1%).
    /// Overrides --chunk-size when set. Enables dynamic per-read chunk sizing.
    /// Clamped to [chunk_min, chunk_max] range (default: 20-500bp).
    #[arg(long)]
    chunk_pct: Option<f64>,

    /// Minimum chunk size clamp for dynamic chunking (overrides default 20bp).
    #[arg(long)]
    chunk_min: Option<usize>,

    /// Maximum chunk size clamp for dynamic chunking (overrides default 500bp).
    #[arg(long)]
    chunk_max: Option<usize>,

    /// Minimum fraction of chunks that must agree (0.0-1.0, default: 0.0 = no threshold).
    /// Example: 0.6 requires 60% of chunks to agree before accepting a mapping.
    #[arg(long)]
    chunk_vote_threshold: Option<f64>,

    /// Number of top genomes to return per read in chunk-based mode (default: 1).
    /// Use 2-3 for multi-genome uncertainty scenarios.
    #[arg(long)]
    chunk_top_n: Option<usize>,

    /// Anchor strategy for chunk-based mapping: rarest (default), golden (quality-weighted), spaced (spaced seed)
    #[arg(long, default_value = "rarest")]
    chunk_strategy: String,

    /// Score aggregation mode for chunk-based mapping: quality (default, score*score), base (raw sum, like JNI Android)
    #[arg(long, default_value = "quality")]
    score_mode: String,

    /// Minimum anchor score threshold for chunk-based mapping (default: 0.5).
    /// Use 0.0 to match JNI Android behavior (no filtering).
    #[arg(long, default_value = "0.5")]
    anchor_min_score: f64,

    /// Use full map_read pipeline for chunks (like JNI Android) instead of anchor_filter.
    #[arg(long)]
    anchor_filter: bool,

    /// Enable SNP-aware scoring for strain resolution.
    #[arg(long)]
    snp_detect: bool,

    /// Minimum support count for SNP detection (default: 3).
    #[arg(long, default_value = "3")]
    snp_min_support: u32,

    /// Enable homopolymer fingerprint scoring for strain resolution
    #[arg(long)]
    hf: bool,

    /// Minimum run length for homopolymer fingerprint (default: 3)
    #[arg(long, default_value = "3")]
    hf_min: usize,

    /// Output BAM format instead of SAM
    #[arg(long)]
    bam: bool,

    /// Stream reads in chunks (limits memory usage for large FASTQ files)
    #[arg(long)]
    stream: bool,

    /// Max RAM to use for streaming (e.g., "32G", "16GB"). Auto-calculates chunk size.
    #[arg(long)]
    max_ram: Option<String>,

    /// Diagnose unmapped reads (sample up to 1000, report why they failed)
    #[arg(long)]
    diagnose_unmapped: bool,

    /// Two-pass mapping: re-map unmapped reads with lower threshold (default: 0.5)
    #[arg(long)]
    two_pass: bool,

    /// Minimum score for second pass (default: 0.5)
    #[arg(long, default_value = "0.5")]
    second_pass_score: f64,

    /// Context window size for alignment (flanking bases, default: 50)
    #[arg(long, default_value = "50")]
    context_window: usize,
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

#[derive(clap::Args)]
struct EmArgs {
    /// Input SAM file (bit-pop mapping output)
    #[arg(short, long, required = true)]
    input: PathBuf,

    /// Output SAM file with EM improved classifications
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// Convergence threshold (KL divergence)
    #[arg(long, default_value = "0.001")]
    convergence: f64,

    /// Maximum EM iterations
    #[arg(long, default_value = "50")]
    max_iterations: usize,

    /// Softmax temperature (lower = sharper)
    #[arg(long, default_value = "0.1")]
    temperature: f64,

    /// Top-K genomes per read for EM
    #[arg(long, default_value = "10")]
    top_k: usize,

    /// Minimum probability to apply EM reassignment (0.0 = always apply, 0.75 = only high confidence)
    #[arg(long, default_value = "0.0")]
    confidence_threshold: f64,
}

#[derive(clap::Args)]
struct ConsensusArgs {
    /// List of index files (comma-separated), e.g. "index_k10.bitpop index_k20.bitpop index_k50.bitpop"
    #[arg(short, long, required = true, num_args = 1..)]
    indexes: Vec<String>,

    /// Reads file (FASTQ)
    #[arg(short, long, required = true)]
    reads: PathBuf,

    /// Output SAM file
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// Strategy: weighted_score (default), majority, best_score (union), base_score (raw sum, no k-weight)
    #[arg(long, default_value = "weighted_score")]
    strategy: String,

    /// Minimum alignment score threshold (0.0-1.0, 0 = no filter)
    #[arg(long, default_value = "0.0")]
    min_score: f64,

    /// Chunk size for long reads (0 = no chunking)
    #[arg(long, default_value = "0")]
    chunk_size: usize,

    /// Chunk size as percentage of read length (0.0-1.0, 0 = disabled)
    #[arg(long, default_value = "0.0")]
    chunk_pct: f64,

    /// Minimum chunk size in bp (default: 20)
    #[arg(long, default_value = "20")]
    chunk_min: usize,

    /// Maximum chunk size in bp (default: 500)
    #[arg(long, default_value = "500")]
    chunk_max: usize,

    /// Enable SNP detection
    #[arg(long)]
    snp_detect: bool,

    /// SNP minimum support count
    #[arg(long, default_value = "3")]
    snp_min_support: u32,

    /// SNP penalty value
    #[arg(long, default_value = "0.1")]
    snp_penalty: f64,

    /// Minimum k-values that must find a mapping (0 = any, default: 1)
    #[arg(long, default_value = "1")]
    min_k_mappings: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Number of top candidates to output per read (0 = only winner, default: 1)
    #[arg(long, default_value = "1")]
    top_n: usize,

    /// Output BAM format instead of SAM
    #[arg(long)]
    bam: bool,

    /// Stream reads in chunks (limits memory usage)
    #[arg(long)]
    stream: bool,

    /// Max RAM to use for streaming (e.g., "32G", "16G"). Auto-calculates chunk size.
    #[arg(long)]
    max_ram: Option<String>,

    /// Two-pass mode: map each k separately (faster, like Python script)
    #[arg(long)]
    two_pass: bool,

    /// Minimum anchor score threshold for chunk-based mapping (default: 0.5)
    #[arg(long, default_value = "0.5")]
    anchor_min_score: f64,

    /// Use full map_read pipeline for chunks (like JNI Android)
    #[arg(long)]
    anchor_filter: bool,
}

#[derive(clap::Args)]
struct ConConArgs {
    /// List of index files, e.g. "index_k10.bitpop index_k20.bitpop"
    #[arg(short, long, required = true, num_args = 1..)]
    indexes: Vec<String>,

    /// Reads file (FASTQ)
    #[arg(short, long, required = true)]
    reads: PathBuf,

    /// Output SAM file
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// Strategy: weighted_score (default), majority, best_score, base_score (raw sum, no k-weight)
    #[arg(long, default_value = "weighted_score")]
    strategy: String,

    /// Minimum alignment score threshold (0.0-1.0, 0 = no filter)
    #[arg(long, default_value = "0.0")]
    min_score: f64,

    /// Minimum k-values that must find a mapping (default: 1)
    #[arg(long, default_value = "1")]
    min_k_mappings: usize,

    /// Number of threads per map
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Top-N rarest k-mer anchors per map (default: 1)
    #[arg(long, default_value = "1")]
    top_n: usize,

    /// Top-N consensus candidates per read (default: 1, >1 for EM)
    #[arg(long, default_value = "1")]
    consensus_top_n: usize,

    /// Chunk size for long reads (0 = no chunking)
    #[arg(long, default_value = "0")]
    chunk_size: usize,

    /// Chunk size as percentage of read length (0.0-1.0, 0 = disabled)
    #[arg(long, default_value = "0.0")]
    chunk_pct: f64,

    /// Minimum chunk size in bp (default: 20)
    #[arg(long, default_value = "20")]
    chunk_min: usize,

    /// Maximum chunk size in bp (default: 500)
    #[arg(long, default_value = "500")]
    chunk_max: usize,

    /// Path to bit-pop executable (auto-detected)
    #[arg(long)]
    bit_pop: Option<PathBuf>,

    /// Context window size for alignment (flanking bases, default: 50)
    #[arg(long, default_value = "50")]
    context_window: usize,

    /// Minimum anchor score threshold for chunk-based mapping (default: 0.5)
    #[arg(long, default_value = "0.5")]
    anchor_min_score: f64,

    /// Use full map_read pipeline for chunks (like JNI Android)
    #[arg(long)]
    anchor_filter: bool,
}

#[derive(clap::Args)]
struct ChunkConsensusArgs {
    /// Index file (.bitpop)
    #[arg(short, long, required = true)]
    index: PathBuf,

    /// Reads file (FASTQ)
    #[arg(short, long, required = true)]
    reads: PathBuf,

    /// Output SAM file
    #[arg(short, long, required = true)]
    output: PathBuf,

    /// Chunk percentages as fraction of read length, comma-separated (e.g. "0.01,0.10,0.50")
    #[arg(short = 'c', long, required = true)]
    chunk_pcts: String,

    /// Voting strategy: majority (default), weighted_score, base_score (raw sum, no chunk-weight)
    #[arg(long, default_value = "majority")]
    strategy: String,

    /// Minimum alignment score threshold (0.0-1.0)
    #[arg(long, default_value = "0.5")]
    min_score: f64,

    /// Minimum configs that must agree (default: majority = N/2 + 1)
    #[arg(long)]
    min_agreement: Option<usize>,

    /// Minimum chunk size in bp (default: 20)
    #[arg(long, default_value = "20")]
    chunk_min: usize,

    /// Maximum chunk size in bp (default: 500)
    #[arg(long, default_value = "500")]
    chunk_max: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Number of top candidates per read (default: 1)
    #[arg(long, default_value = "1")]
    top_n: usize,

    /// Output BAM format instead of SAM
    #[arg(long)]
    bam: bool,

    /// Stream reads in chunks (limits memory usage for large FASTQ files)
    #[arg(long)]
    stream: bool,

    /// Max RAM to use for streaming (e.g., "32G", "16GB"). Auto-calculates chunk size.
    #[arg(long)]
    max_ram: Option<String>,

    /// Minimum anchor score threshold for chunk-based mapping (default: 0.5)
    #[arg(long, default_value = "0.5")]
    anchor_min_score: f64,

    /// Use full map_read pipeline for chunks (like JNI Android)
    #[arg(long)]
    anchor_filter: bool,
}

#[derive(clap::Args)]
struct TaxArgs {
    /// Input SAM file (bit-pop mapping output)
    #[arg(short, long, required = true)]
    input: PathBuf,

    /// Path to NCBI nodes.dmp file
    #[arg(long, required = true)]
    nodes_dmp: PathBuf,

    /// Path to NCBI names.dmp file
    #[arg(long, required = true)]
    names_dmp: PathBuf,

    /// Output file for taxonomic report (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of top entries per rank to display (default: 10)
    #[arg(long, default_value = "10")]
    top_n: usize,

    /// Output format: text, json
    #[arg(long, default_value = "text")]
    format: String,
}

#[tokio::main]
async fn main() {
    let cli = Cli::parse();

    if let Err(e) = match cli.command {
        Commands::Run(args) => cmd_run(&args).await,
        Commands::Build(args) => {
            cmd_build(&args, cli.verbose);
            Ok(())
        }
        Commands::Map(args) => {
            cmd_map(&args, cli.verbose);
            Ok(())
        }
        Commands::Load(args) => {
            cmd_load(&args, cli.verbose);
            Ok(())
        }
        Commands::Stats(args) => {
            cmd_stats(&args, cli.verbose);
            Ok(())
        }
        Commands::Search(args) => {
            cmd_search(&args, cli.verbose).await;
            Ok(())
        }
        Commands::Fetch(args) => cmd_fetch(&args, cli.verbose).await,
        Commands::Update(args) => cmd_update(&args, cli.verbose).await,
        Commands::Em(args) => {
            cmd_em(&args);
            Ok(())
        }
        Commands::Consensus(args) => {
            cmd_consensus(&args);
            Ok(())
        }
        Commands::ConCon(args) => {
            cmd_concon(&args);
            Ok(())
        }
        Commands::ChunkConsensus(args) => {
            cmd_chunk_consensus(&args);
            Ok(())
        }
        Commands::Tax(args) => {
            cmd_tax(&args);
            Ok(())
        }
    } {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn expand_genome_paths(paths: &[PathBuf]) -> Vec<PathBuf> {
    let mut expanded = Vec::new();
    for path in paths {
        if path.is_dir() {
            match std::fs::read_dir(path) {
                Ok(entries) => {
                    let files: Vec<_> = entries
                        .filter_map(|e| e.ok())
                        .map(|e| e.path())
                        .filter(|p| {
                            p.extension()
                                .map(|e| e == "fna" || e == "fasta" || e == "fa")
                                .unwrap_or(false)
                        })
                        .collect();
                    println!(
                        "  Found {} genome file(s) in {}",
                        files.len(),
                        path.display()
                    );
                    expanded.extend(files);
                }
                Err(e) => {
                    eprintln!("Cannot read directory {}: {}", path.display(), e);
                }
            }
        } else {
            expanded.push(path.clone());
        }
    }
    expanded
}

fn cmd_build(args: &BuildArgs, verbose: bool) {
    let start = Instant::now();

    println!("Building FM-Index...");

    let fasta_paths = expand_genome_paths(&args.fasta);
    if fasta_paths.is_empty() {
        eprintln!("No genome files found.");
        std::process::exit(1);
    }

    let mut bp = BitPop::new(args.k);
    bp.set_auto_k(args.auto_k);
    bp.set_read_type(&args.read_type);

    if args.method != "none" {
        let fuzzy_method = match args.method.as_str() {
            "fuzzy-kmer" => FuzzyMethod::FuzzyKmer,
            "fuzzy-seed" => FuzzyMethod::FuzzySeed,
            "neighborhood" => FuzzyMethod::Neighborhood,
            _ => FuzzyMethod::None,
        };
        bp.set_fuzzy_method(fuzzy_method);
        bp.set_fuzzy_mismatches(args.fuzzy_mismatches);
        if verbose {
            println!(
                "  Fuzzy method: {} (mismatches: {})",
                args.method, args.fuzzy_mismatches
            );
        }
    }

    if args.spaced_seed {
        bp.set_spaced_seed(true);
        if let Some(pattern) = &args.spaced_seed_pattern {
            bp.set_spaced_seed_pattern(pattern);
            if verbose {
                println!("  Spaced seed pattern: {}", pattern);
            }
        }
        if verbose {
            println!(
                "  Spaced seed: enabled (pattern: {})",
                bp.spaced_seed_pattern()
            );
        }
    }

    if args.hf {
        bp.set_hf(true);
        bp.set_hf_min_run(args.hf_min);
        if verbose {
            println!(
                "  Homopolymer fingerprint: enabled (min_run: {})",
                args.hf_min
            );
        }
    }

    bp.set_search_radius(args.search_radius);
    let mut total_bases: usize = 0;

    let pb = ProgressBar::new(fasta_paths.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner} {msg}: [{elapsed_precise} {bar:40} {pos}/{len}]")
            .unwrap(),
    );

    let mut loaded = 0;
    for fasta_path in &fasta_paths {
        let path_str = fasta_path.to_string_lossy().to_string();
        pb.set_message(format!("Loading: {}", path_str));

        let gids = if args.cami {
            let basename = fasta_path
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default();
            let cami_name = extract_cami_genome_name(&basename);

            if verbose {
                println!("  CAMI mode: {} -> {}", path_str, cami_name);
            }

            let seqs = bit_pop::fasta::read_all_sequences(&path_str);
            match seqs {
                Ok(sequences) => {
                    let mut ids = Vec::new();
                    for (_header, seq) in sequences {
                        let gid = bp.add_genome(&cami_name, &seq);
                        ids.push(gid);
                    }
                    ids
                }
                Err(e) => {
                    pb.finish_with_message(format!("Error reading {}: {}", path_str, e));
                    eprintln!("Error reading {}: {}", path_str, e);
                    std::process::exit(1);
                }
            }
        } else if args.pacbio {
            let basename = fasta_path
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default();
            let pacbio_name = extract_pacbio_genome_name(&basename);

            if verbose {
                println!("  PacBio mode: {} -> {}", path_str, pacbio_name);
            }

            let seqs = bit_pop::fasta::read_all_sequences(&path_str);
            match seqs {
                Ok(sequences) => {
                    let mut ids = Vec::new();
                    for (_header, seq) in sequences {
                        let gid = bp.add_genome(&pacbio_name, &seq);
                        ids.push(gid);
                    }
                    ids
                }
                Err(e) => {
                    pb.finish_with_message(format!("Error reading {}: {}", path_str, e));
                    eprintln!("Error reading {}: {}", path_str, e);
                    std::process::exit(1);
                }
            }
        } else {
            #[cfg(feature = "mmap")]
            let ids = if args.mmap {
                bp.load_genome_fasta_mmap(&path_str)
            } else {
                bp.load_genome_fasta(&path_str)
            };
            #[cfg(not(feature = "mmap"))]
            let ids = bp.load_genome_fasta(&path_str);
            ids.unwrap_or_default()
        };

        for gid in gids {
            let seq_len = bp.genome_seq_len(gid).unwrap_or(0);
            total_bases += seq_len;
            if verbose {
                if let Some(name) = bp.genome_name(gid) {
                    println!("    Added genome: {} ({} bases)", name, seq_len);
                }
            }
        }
        loaded += 1;
        pb.inc(1);
        report_atomic_progress(loaded, fasta_paths.len() as u64);
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
                std::fs::metadata(&args.output)
                    .map(|m| m.len())
                    .unwrap_or(0),
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

/// Second pass: re-map unmapped reads with lower threshold.
fn cmd_map_second_pass(
    bp: &BitPop,
    reads_path: &str,
    _output_path: &str,
    chunk_size: usize,
    _threads: usize,
    align_mode: AlignMode,
    min_score: f64,
) -> usize {
    use bit_pop::fastq::FastqChunkParser;

    let total = match FastqChunkParser::count_reads(reads_path) {
        Ok(n) => n,
        Err(_) => return 0,
    };

    let mut parser = match FastqChunkParser::new(reads_path, chunk_size) {
        Ok(p) => p,
        Err(_) => return 0,
    };

    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner} 2nd pass: [{elapsed_precise} {bar:40}] {pos}/{len}")
            .unwrap(),
    );

    let mut mapped = 0usize;

    // Collect unmapped reads and re-map
    while let Some(chunk) = parser.next_chunk().unwrap() {
        for (_name, seq, _qual) in &chunk {
            // First check if this read already mapped in first pass
            let first_results = bp.map_read_with_mode(seq, align_mode, 50);
            if !first_results.is_empty() {
                pb.inc(1);
                continue;
            }

            // Try with lower threshold
            let second_results = bp.anchor_filter_with_mode(seq, align_mode, min_score, 100);
            if !second_results.is_empty() {
                mapped += 1;
            }
            pb.inc(1);
        }
    }

    pb.finish_with_message(format!("2nd pass: {} reads mapped", mapped));
    mapped
}

/// Configuration for stream mapping.
struct StreamMapConfig {
    reads_path: String,
    output_path: String,
    chunk_size: usize,
    threads: usize,
    align_mode: AlignMode,
    golden_anchors: bool,
    min_quality: u8,
    write_bam: bool,
    use_chunking: bool,
    diagnose: bool,
    second_pass_score: f64,
}

/// Stream map: process reads in chunks to limit memory usage.
fn cmd_map_stream(bp: &BitPop, config: &StreamMapConfig) -> usize {
    use bit_pop::fastq::FastqChunkParser;

    let total = match FastqChunkParser::count_reads(&config.reads_path) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("Failed to count reads: {}", e);
            std::process::exit(1);
        }
    };
    println!(
        "Total reads: {} (streaming, chunk={})",
        total, config.chunk_size
    );

    let mut parser = match FastqChunkParser::new(&config.reads_path, config.chunk_size) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Failed to open FASTQ: {}", e);
            std::process::exit(1);
        }
    };

    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner} Mapping: [{elapsed_precise} {bar:40}] {pos}/{len} {msg}")
            .unwrap(),
    );

    let mut mapped = 0usize;
    let mut unmapped = 0usize;
    let mut chunk_num = 0usize;
    let mut processed = 0u64;

    while let Some(chunk) = parser.next_chunk().unwrap() {
        chunk_num += 1;
        if chunk.is_empty() {
            break;
        }

        let chunk_reads: Vec<(&str, &str)> = chunk
            .iter()
            .map(|(n, s, _)| (n.as_str(), s.as_str()))
            .collect();
        let chunk_start = processed;

        let chunk_output = format!("{}_chunk{}.tmp", config.output_path, chunk_num);

        if config.threads > 1 {
            let pb_clone = pb.clone();
            let result = if config.use_chunking {
                bp.map_reads_with_chunking_parallel_with_progress(
                    &chunk_reads,
                    &chunk_output,
                    50,
                    move |completed, _total| {
                        pb_clone.set_position(chunk_start + completed as u64);
                    },
                )
                .unwrap_or(0)
            } else {
                bp.map_reads_parallel_with_progress(
                    &chunk_reads,
                    &chunk_output,
                    50,
                    if chunk.len() > 1000 { 100 } else { 10 },
                    move |completed, _total| {
                        pb_clone.set_position(chunk_start + completed as u64);
                    },
                )
                .unwrap_or(0)
            };
            mapped += result;
        } else {
            for (_name, seq, qual) in &chunk {
                let results = if config.golden_anchors {
                    bp.map_read_with_golden_anchors(seq, qual, config.align_mode, 50)
                } else {
                    bp.map_read_with_quality_mode(
                        seq,
                        qual,
                        config.align_mode,
                        config.min_quality,
                        50,
                    )
                };
                if !results.is_empty() {
                    mapped += 1;
                } else {
                    unmapped += 1;
                }
                processed += 1;
                pb.set_position(processed);
                report_atomic_progress(processed, total as u64);
            }
        }

        println!("  Chunk {} done: {} reads", chunk_num, chunk.len());
        drop(chunk);
    }

    pb.finish_with_message("Mapping complete");

    // Two-pass: re-map unmapped with lower threshold
    if unmapped > 0 {
        println!(
            "\nSecond pass: remapping {} unmapped reads with lower threshold",
            unmapped
        );
        let second_mapped = cmd_map_second_pass(
            bp,
            &config.reads_path,
            &config.output_path,
            config.chunk_size,
            config.threads,
            config.align_mode,
            config.second_pass_score,
        );
        mapped += second_mapped;
        unmapped -= second_mapped;
    }

    // Diagnostic output
    if config.diagnose && unmapped > 0 {
        println!("\nUnmapped read diagnostics (sample: 1000 reads):");
        let mut reasons: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        let mut sample_count = 0;

        let mut parser = bit_pop::fastq::FastqChunkParser::new(&config.reads_path, 1000).unwrap();
        while let Some(chunk) = parser.next_chunk().unwrap() {
            for (_name, seq, _qual) in &chunk {
                let results = bp.map_read_with_mode(seq, config.align_mode, 50);
                if results.is_empty() {
                    let reason = bp.diagnose_read(seq);
                    *reasons.entry(reason).or_insert(0) += 1;
                    sample_count += 1;
                    if sample_count >= 1000 {
                        break;
                    }
                }
            }
            if sample_count >= 1000 {
                break;
            }
        }

        let mut sorted: Vec<_> = reasons.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (reason, count) in sorted {
            println!("  [{:5}] {}", count, reason);
        }
    }

    // Merge chunk outputs
    if chunk_num > 1 {
        merge_sam_chunks(&config.output_path, chunk_num, config.write_bam);
    } else if chunk_num == 1 {
        let chunk_output = format!("{}_chunk1.tmp", config.output_path);
        if std::path::Path::new(&chunk_output).exists() {
            std::fs::rename(&chunk_output, &config.output_path).unwrap();
        }
    }

    mapped
}

/// Merge SAM chunk files into final output.
fn merge_sam_chunks(output_path: &str, num_chunks: usize, write_bam: bool) {
    let mut lines: Vec<String> = Vec::new();
    let mut header_done = false;

    for i in 1..=num_chunks {
        let chunk_path = format!("{}_chunk{}.tmp", output_path, i);
        if let Ok(content) = std::fs::read_to_string(&chunk_path) {
            for line in content.lines() {
                if line.starts_with('@') {
                    if !header_done || line.starts_with("@HD") {
                        lines.push(line.to_string());
                        if line.starts_with("@HD") {
                            header_done = true;
                        }
                    }
                } else {
                    lines.push(line.to_string());
                }
            }
        }
        std::fs::remove_file(&chunk_path).ok();
    }

    if write_bam {
        let temp_sam = format!("{}.tmp.sam", output_path);
        std::fs::write(&temp_sam, lines.join("\n").as_str()).unwrap();
        if let Ok(mut subprocess) = std::process::Command::new("samtools")
            .args(["view", "-bS", &temp_sam, "-o", output_path])
            .spawn()
        {
            subprocess.wait().ok();
        }
        std::fs::remove_file(&temp_sam).ok();
    } else {
        std::fs::write(output_path, lines.join("\n").as_str()).unwrap();
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
    println!("  top_n: {}", bp.top_n());

    if args.method != "none" {
        let fuzzy_method = match args.method.as_str() {
            "fuzzy-kmer" => FuzzyMethod::FuzzyKmer,
            "fuzzy-seed" => FuzzyMethod::FuzzySeed,
            "neighborhood" => FuzzyMethod::Neighborhood,
            _ => FuzzyMethod::None,
        };
        bp.set_fuzzy_method(fuzzy_method);
        bp.set_fuzzy_mismatches(args.fuzzy_mismatches);
        if verbose {
            println!(
                "  Fuzzy method: {} (mismatches: {})",
                args.method, args.fuzzy_mismatches
            );
        }
    }

    if args.spaced_seed {
        bp.set_spaced_seed(true);
        if let Some(pattern) = &args.spaced_seed_pattern {
            bp.set_spaced_seed_pattern(pattern);
            if verbose {
                println!("  Spaced seed pattern: {}", pattern);
            }
        }
        if verbose {
            println!(
                "  Spaced seed: enabled (pattern: {})",
                bp.spaced_seed_pattern()
            );
        }
    }

    bp.set_search_radius(args.search_radius);

    if let Some(chunk_size) = args.chunk_size {
        bp.set_chunk_size(chunk_size);
        println!("  Chunk size: {}bp (PacBio mode)", chunk_size);
    }

    if let Some(chunk_pct) = args.chunk_pct {
        bp.set_chunk_pct(chunk_pct);
    }

    if let Some(chunk_min) = args.chunk_min {
        bp.set_chunk_min(chunk_min);
    }

    if let Some(chunk_max) = args.chunk_max {
        bp.set_chunk_max(chunk_max);
    }

    if bp.chunk_pct() > 0.0 {
        println!(
            "  Chunk pct: {:.2}% (dynamic, clamped {}-{}bp)",
            bp.chunk_pct() * 100.0,
            bp.chunk_min(),
            bp.chunk_max()
        );
    }

    if let Some(threshold) = args.chunk_vote_threshold {
        bp.set_chunk_vote_threshold(threshold);
        println!(
            "  Chunk vote threshold: {:.0}% (requires {:.0}% chunk agreement)",
            threshold * 100.0,
            threshold * 100.0
        );
    }

    if let Some(top_n) = args.chunk_top_n {
        bp.set_chunk_top_n(top_n);
        println!("  Chunk top-N: {} genomes per read", top_n);
    }

    let chunk_strategy = match args.chunk_strategy.as_str() {
        "golden" => bit_pop::ChunkAnchorStrategy::Golden,
        "spaced" => bit_pop::ChunkAnchorStrategy::Spaced,
        _ => bit_pop::ChunkAnchorStrategy::Rarest,
    };
    bp.set_chunk_anchor_strategy(chunk_strategy);
    if args.chunk_strategy != "rarest" {
        println!("  Chunk strategy: {}", args.chunk_strategy);
    }

    let score_mode = match args.score_mode.as_str() {
        "base" => bit_pop::ChunkScoreMode::Base,
        _ => bit_pop::ChunkScoreMode::Quality,
    };
    bp.set_chunk_score_mode(score_mode);
    if args.score_mode != "quality" {
        println!("  Score mode: {}", args.score_mode);
    }

    bp.set_chunk_anchor_min_score(args.anchor_min_score);
    if args.anchor_min_score != 0.5 {
        println!("  Anchor min score: {}", args.anchor_min_score);
    }

    bp.set_chunk_min_score(args.min_score);
    if args.min_score > 0.0 {
        println!("  Chunk min score: {}", args.min_score);
    }

    if args.anchor_filter {
        bp.set_chunk_use_anchor_filter(true);
        println!("  Anchor filter: enabled (legacy mode)");
    }

    if args.snp_detect {
        bp.set_snp_detect(true);
        bp.set_snp_min_support(args.snp_min_support);
        println!(
            "  SNP detection: enabled (min support: {}, penalty: {})",
            args.snp_min_support,
            bp.snp_penalty()
        );
    }

    if args.hf {
        bp.set_hf(true);
        bp.set_hf_min_run(args.hf_min);
        println!(
            "  Homopolymer fingerprint: enabled (min_run: {})",
            args.hf_min
        );
    }

    let load_time = load_start.elapsed();

    let align_mode = match args.align_mode.as_str() {
        "sw" => AlignMode::Sw,
        "hybrid" => AlignMode::Hybrid,
        "softclip" => AlignMode::Softclip,
        "chain" => AlignMode::Chain,
        _ => AlignMode::Xor,
    };
    bp.set_align_mode(align_mode);

    println!(
        "Index loaded in {:.3}s ({})\n",
        load_time.as_secs_f64(),
        bp.genome_count()
    );
    if verbose {
        println!("Alignment mode: {}\n", align_mode);
    }

    // Check for paired-end mode
    if let (Some(r1_path), Some(r2_path)) = (&args.reads_1, &args.reads_2) {
        cmd_map_paired(
            &bp,
            r1_path,
            r2_path,
            &args.output,
            args.min_quality,
            args.bam,
        );
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

    if args.stream {
        let chunk_size = if let Some(ref max_ram) = args.max_ram {
            parse_stream_chunk_size(&Some(max_ram.clone()))
        } else {
            20_000_000
        };
        println!("Streaming mode: chunk size = {} reads", chunk_size);

        let map_start = Instant::now();
        let mapped_count = cmd_map_stream(
            &bp,
            &StreamMapConfig {
                reads_path,
                output_path: args.output.to_str().unwrap().to_string(),
                chunk_size,
                threads: args.reads_threads,
                align_mode,
                golden_anchors: args.golden_anchors,
                min_quality: args.min_quality,
                write_bam: args.bam,
                use_chunking: args.chunk_size.is_some() || args.chunk_pct.is_some(),
                diagnose: args.diagnose_unmapped,
                second_pass_score: args.second_pass_score,
            },
        );

        let elapsed = start.elapsed();
        println!("\nMapping complete: {} reads mapped", mapped_count);
        println!("  Alignment mode: {}", align_mode);
        println!("  Load time:  {:.3}s", load_time.as_secs_f64());
        println!("  Map time:   {:.2}s", map_start.elapsed().as_secs_f64());
        println!("  Total time: {:.2}s", elapsed.as_secs_f64());
        return;
    }

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
                println!(
                    "Quality filter (min Q{}): {}/{} reads passed",
                    args.min_quality,
                    passed.len(),
                    reads.len()
                );
                passed
                    .iter()
                    .map(|&i| (reads[i].0.clone(), reads[i].1.clone()))
                    .collect()
            }
            ReadsFormat::Fasta(_) => {
                println!("Warning: quality filtering ignored for FASTA input");
                reads_format
                    .iter_fasta()
                    .map(|(n, s)| (n.to_string(), s.to_string()))
                    .collect()
            }
        }
    } else {
        reads_format
            .iter_fasta()
            .map(|(n, s)| (n.to_string(), s.to_string()))
            .collect()
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

                let genome_name_refs: Vec<&str> =
                    genomes_owned.iter().map(|(n, _)| n.as_str()).collect();
                let genome_header: Vec<(&str, usize)> = genomes_owned
                    .iter()
                    .map(|(n, l)| (n.as_str(), *l))
                    .collect();

                let name_refs: Vec<&str> = genome_name_refs.clone();
                let total = reads.len();

                let pb = ProgressBar::new(total as u64);
                pb.set_style(ProgressStyle::default_bar()
                    .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                    .unwrap());

                let ctx_window = args.context_window;
                let mapped: Vec<(String, String, Vec<bit_pop::QualityMappingResult>)> = reads
                    .iter()
                    .enumerate()
                    .map(|(i, (name, seq, qual))| {
                        let results = if args.golden_anchors {
                            bp.map_read_with_golden_anchors(seq, qual, align_mode, ctx_window)
                        } else {
                            bp.map_read_with_quality_mode(
                                seq,
                                qual,
                                align_mode,
                                args.min_quality,
                                ctx_window,
                            )
                        };
                        if (i + 1) % 10 == 0 || i + 1 == total {
                            pb.set_position((i + 1) as u64);
                            pb.set_message(format!("{}/{} reads", i + 1, total));
                            report_atomic_progress((i + 1) as u64, total as u64);
                        }
                        (name.clone(), seq.clone(), results)
                    })
                    .collect();

                pb.finish_with_message("Mapping complete");

                let mut writer =
                    bit_pop::sam::SamWriter::new(args.output.to_str().unwrap()).unwrap();
                writer.write_header(&genome_header).unwrap();

                let mut mapped_count = 0;
                for (name, seq, results) in &mapped {
                    writer
                        .write_quality_mappings(name, seq, results, &name_refs)
                        .unwrap();
                    if !results.is_empty() {
                        mapped_count += 1;
                    }
                }

                mapped_count
            }
            ReadsFormat::Fasta(_) => {
                let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta
                    .iter()
                    .map(|(name, seq)| (name.as_str(), seq.as_str()))
                    .collect();
                let total = reads_refs.len();

                if args.reads_threads > 1 {
                    let pb = ProgressBar::new(total as u64);
                    pb.set_style(ProgressStyle::default_bar()
                        .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                        .unwrap());
                    let pb_clone = pb.clone();

                    let result = bp
                        .map_reads_parallel_with_progress(
                            &reads_refs,
                            args.output.to_str().unwrap(),
                            args.context_window,
                            if total > 1000 { 100 } else { 10 },
                            move |completed, total| {
                                pb_clone.set_position(completed as u64);
                                pb_clone.set_message(format!("{}/{} reads", completed, total));
                                report_atomic_progress(completed as u64, total as u64);
                            },
                        )
                        .unwrap_or(0);

                    pb.finish_with_message("Mapping complete");
                    result
                } else {
                    let pb = ProgressBar::new(total as u64);
                    pb.set_style(ProgressStyle::default_bar()
                        .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                        .unwrap());
                    let pb_clone = pb.clone();

                    let result = bp
                        .map_reads_to_output_with_progress(
                            &reads_refs,
                            args.output.to_str().unwrap(),
                            args.context_window,
                            if total > 1000 { 100 } else { 10 },
                            args.bam,
                            move |completed, total| {
                                pb_clone.set_position(completed as u64);
                                pb_clone.set_message(format!("{}/{} reads", completed, total));
                                report_atomic_progress(completed as u64, total as u64);
                            },
                        )
                        .unwrap_or(0);

                    pb.finish_with_message("Mapping complete");
                    result
                }
            }
        }
    } else if args.reads_threads > 1 {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        let total = reads_refs.len();

        let pb = ProgressBar::new(total as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap(),
        );
        let pb_clone = pb.clone();

        let result = if bp.chunk_size() > 0 || bp.chunk_pct() > 0.0 {
            bp.map_reads_with_chunking_parallel_with_progress(
                &reads_refs,
                args.output.to_str().unwrap(),
                args.context_window,
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                    report_atomic_progress(completed as u64, total as u64);
                },
            )
            .unwrap_or(0)
        } else {
            bp.map_reads_parallel_with_progress(
                &reads_refs,
                args.output.to_str().unwrap(),
                args.context_window,
                if total > 1000 { 100 } else { 10 },
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                    report_atomic_progress(completed as u64, total as u64);
                },
            )
            .unwrap_or(0)
        };

        pb.finish_with_message("Mapping complete");
        result
    } else {
        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();
        let total = reads_refs.len();

        let pb = ProgressBar::new(total as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap(),
        );
        let pb_clone = pb.clone();

        let result = if bp.chunk_size() > 0 || bp.chunk_pct() > 0.0 {
            bp.map_reads_with_chunking_parallel_with_progress(
                &reads_refs,
                args.output.to_str().unwrap(),
                args.context_window,
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                    report_atomic_progress(completed as u64, total as u64);
                },
            )
            .unwrap_or(0)
        } else {
            bp.map_reads_to_output_with_progress(
                &reads_refs,
                args.output.to_str().unwrap(),
                args.context_window,
                if total > 1000 { 100 } else { 10 },
                args.bam,
                move |completed, total| {
                    pb_clone.set_position(completed as u64);
                    pb_clone.set_message(format!("{}/{} reads", completed, total));
                    report_atomic_progress(completed as u64, total as u64);
                },
            )
            .unwrap_or(0)
        };

        pb.finish_with_message("Mapping complete");
        result
    };

    let total_reads = filtered_reads_fasta.len();
    let unmapped_reads = total_reads - mapped_count;
    let mut mapped_count = mapped_count;

    // Diagnostic output
    if args.diagnose_unmapped && unmapped_reads > 0 {
        println!("\nUnmapped read diagnostics (sample: 1000 reads):");
        let mut reasons: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        let mut sample_count = 0;

        for (_name, seq) in &filtered_reads_fasta {
            if sample_count >= 1000 {
                break;
            }
            let results = bp.map_read_with_mode(seq, align_mode, 50);
            if results.is_empty() {
                let reason = bp.diagnose_read(seq);
                *reasons.entry(reason).or_insert(0) += 1;
                sample_count += 1;
            }
        }

        let mut sorted: Vec<_> = reasons.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (reason, count) in sorted {
            println!("  [{:5}] {}", count, reason);
        }
    }

    // Two-pass: re-map unmapped with lower threshold + EM refinement
    if args.two_pass && unmapped_reads > 0 {
        println!(
            "\nSecond pass: remapping {} unmapped reads with lower threshold ({})",
            unmapped_reads, args.second_pass_score
        );

        // Phase 1: Collect best mapping per unmapped read - parallel
        let names = bp.genome_names_ordered();
        let second_pass_score = args.second_pass_score;
        let pb = ProgressBar::new(unmapped_reads as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} 2nd pass: [{elapsed_precise} {bar:40} {pos}/{len}]")
                .unwrap(),
        );
        let pb_clone = pb.clone();
        let num_threads = args.reads_threads.max(1);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let second_mappings: Vec<(String, String, String, MappingResult)> = pool.install(|| {
            filtered_reads_fasta
                .par_iter()
                .filter_map(|(name, seq)| {
                    let first_results = bp.map_read_with_mode(seq, align_mode, 50);
                    if first_results.is_empty() {
                        let second_results =
                            bp.map_read_with_threshold(seq, align_mode, 50, second_pass_score);
                        if let Some(best_hit) = second_results.into_iter().next() {
                            let gname = if best_hit.genome_id < names.len() as u32 {
                                names[best_hit.genome_id as usize].clone()
                            } else {
                                "*".to_string()
                            };
                            pb_clone.inc(1);
                            Some((name.clone(), gname, seq.clone(), best_hit))
                        } else {
                            pb_clone.inc(1);
                            None
                        }
                    } else {
                        pb_clone.inc(1);
                        None
                    }
                })
                .collect()
        });

        pb.finish_with_message("2nd pass mapping done");

        // Phase 2: Write mappings to SAM
        let mut second_mapped = 0;
        let output_path = args.output.to_path_buf();
        let mut sam_file = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .append(true)
                .open(&output_path)
                .unwrap(),
        );

        for (name, gname, seq, hit) in &second_mappings {
            second_mapped += 1;
            let pos = hit.position + 1;
            let mapq = ((hit.score * 60.0) as u16).min(60);
            let flag = if hit.is_reverse { 16u16 } else { 0u16 };
            writeln!(
                sam_file,
                "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tMD:Z:{}\tS:f:{}",
                name, flag, gname, pos, mapq, hit.cigar, seq, hit.md_string, hit.score
            )
            .unwrap();
        }
        mapped_count += second_mapped;
        println!("  Second pass mapped: {} additional reads", second_mapped);
    }

    let elapsed = start.elapsed();

    println!(
        "\nMapping complete: {}/{} reads mapped",
        mapped_count, total_reads
    );
    println!("  Alignment mode: {}", align_mode);
    println!("  Load time:  {:.3}s", load_time.as_secs_f64());
    println!("  Map time:   {:.2}s", map_start.elapsed().as_secs_f64());
    println!("  Total time: {:.2}s", elapsed.as_secs_f64());
}

fn cmd_map_paired(
    bp: &BitPop,
    r1_path: &Path,
    r2_path: &Path,
    output: &Path,
    min_quality: u8,
    write_bam: bool,
) {
    let map_start = Instant::now();

    println!("Paired-end mapping mode");
    println!("  R1: {}", r1_path.to_string_lossy());
    println!("  R2: {}", r2_path.to_string_lossy());

    let pairs = match bit_pop::fastq::parse_paired_fastq(
        r1_path.to_str().unwrap(),
        r2_path.to_str().unwrap(),
    ) {
        Ok(pairs) => pairs,
        Err(e) => {
            eprintln!("Error parsing paired FASTQ: {}", e);
            std::process::exit(1);
        }
    };

    println!("Loaded {} read pairs", pairs.len());

    let total_pairs = pairs.len();
    let pb = ProgressBar::new(total_pairs as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner} Mapping pairs: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
            .unwrap(),
    );

    let mapped_count = if min_quality > 0 {
        let result = bp
            .map_paired_reads_parallel_quality(
                &pairs,
                output.to_str().unwrap(),
                min_quality,
                50,
                true,
                5,
                write_bam,
            )
            .unwrap_or(0);

        pb.set_position(total_pairs as u64);
        pb.set_message(format!("{} pairs", total_pairs));
        report_atomic_progress(total_pairs as u64, total_pairs as u64);
        pb.finish_with_message("Mapping complete");
        result
    } else {
        let result = bp
            .map_paired_reads_parallel(&pairs, output.to_str().unwrap(), 50, true, 5, write_bam)
            .unwrap_or(0);

        pb.set_position(total_pairs as u64);
        pb.set_message(format!("{} pairs", total_pairs));
        report_atomic_progress(total_pairs as u64, total_pairs as u64);
        pb.finish_with_message("Mapping complete");
        result
    };

    let elapsed = map_start.elapsed();

    println!(
        "\nPaired-end mapping complete: {} pairs processed",
        mapped_count
    );
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
    println!(
        "Added {} new genomes ({} -> {})",
        new_count - old_count,
        old_count,
        new_count
    );

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
    println!(
        "File size:     {} bytes ({:.1} MB)",
        file_size,
        file_size as f64 / 1_000_000.0
    );
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

    println!(
        "Searching NCBI for: {} ({})",
        args.organism, args.molecule_type
    );

    let search_start = Instant::now();
    let result = match client.search(&format!("{}[Organism]", args.organism)).await {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Search failed: {}", e);
            std::process::exit(1);
        }
    };
    println!(
        "Search completed in {:.2}s",
        search_start.elapsed().as_secs_f64()
    );
    println!("Found {} results", result.count);

    if result.idlist.is_empty() {
        println!("No results found. Try a different organism name or molecule type.");
        return;
    }

    let display_count = result.idlist.len().min(args.max_results);
    println!("\nTop {} results:", display_count);
    println!(
        "{:<25} {:<50} {:<10} Type",
        "Accession", "Description", "Length"
    );
    println!("{:-<100}", "");

    if result.idlist.len() > display_count {
        println!(
            "  ... and {} more results (use -n to increase)",
            result.idlist.len() - display_count
        );
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
    let filtered: Vec<&bit_pop::ncbi::DocSum> = all_docsums
        .iter()
        .filter(|ds| {
            let is_refseq = ds
                .title
                .as_ref()
                .map(|t| t.contains("RefSeq"))
                .unwrap_or(false);
            let is_genomic = ds
                .nuc_genesim
                .as_ref()
                .map(|n| n.contains("Genomic DNA"))
                .unwrap_or(false);
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
            println!(
                "{:<25} {:<50} {:<10} RefSeq",
                accession, title_display, pavg
            );
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
            let base = if parts.len() >= 2 {
                parts[0]
            } else {
                accession
            };
            cache
                .cache_sequence(accession, version, base, &f)
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
                let file_size = std::fs::metadata(&args.output)
                    .map(|m| m.len())
                    .unwrap_or(0);
                for (name, _) in &genomes {
                    if let Some(_genome) = cache.manifest().get(name) {
                        let _ = cache.cache_index(name, &args.output, args.k);
                    }
                }

                println!("\nDone!");
                println!("  Genomes:    {}", genomes.len());
                println!(
                    "  Index size: {} bytes ({:.1} MB)",
                    file_size,
                    file_size as f64 / 1_000_000.0
                );
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

    println!(
        "Checking for updates in {} genome(s)...",
        cache.manifest().len()
    );

    if cache.manifest().is_empty() {
        println!("No genomes cached. Use 'fetch' to download genomes first.");
        return Ok(());
    }

    let mut updated = Vec::new();
    let mut already_current = Vec::new();

    let genomes_list: Vec<(String, String, String, String)> = cache
        .list_genomes()
        .iter()
        .map(|g| {
            (
                g.accession.clone(),
                g.version.clone(),
                g.base_accession.clone(),
                g.checksum.clone(),
            )
        })
        .collect();

    for (acc, version, base_accession, checksum) in genomes_list {
        print!("  Checking {}... ", acc);
        match client.fetch_by_accession_version(&acc).await {
            Ok(fasta) => {
                use sha2::{Digest, Sha256};
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
    use sha2::{Digest, Sha256};
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
        name = name
            .trim_end_matches(".fastq")
            .trim_end_matches(".fasta")
            .to_string();
    }
    parent.join(format!("{}.sam", name))
}

fn find_or_build_index(
    genome_paths: &[PathBuf],
    k: usize,
    auto_k: bool,
    force: bool,
) -> Result<BitPop, String> {
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
            let genome_mtime = std::fs::metadata(genome_path)
                .map_err(|e| e.to_string())?
                .modified()
                .map_err(|e| e.to_string())?;

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

    println!(
        "  Building index ({} genomes, k={})...",
        genome_paths.len(),
        k
    );
    let build_start = Instant::now();

    let mut bp = BitPop::new(k);
    bp.set_auto_k(auto_k);
    for path in genome_paths {
        let path_str = path.to_string_lossy();
        let ids = bp
            .load_genome_fasta(&path_str)
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
fn find_or_build_index_mmap(
    genome_paths: &[PathBuf],
    k: usize,
    auto_k: bool,
    force: bool,
) -> Result<BitPop, String> {
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
            let genome_mtime = std::fs::metadata(genome_path)
                .map_err(|e| e.to_string())?
                .modified()
                .map_err(|e| e.to_string())?;

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

    println!(
        "  Building index (mmap, {} genomes, k={})...",
        genome_paths.len(),
        k
    );
    let build_start = Instant::now();

    let mut bp = BitPop::new(k);
    bp.set_auto_k(auto_k);
    for path in genome_paths {
        let path_str = path.to_string_lossy();
        let ids = bp
            .load_genome_fasta_mmap(&path_str)
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

            let accessions =
                if genome.starts_with("NC_") || genome.starts_with("AC_") || genome.contains('.') {
                    vec![genome.clone()]
                } else {
                    let search_result = client
                        .search(&format!("{}[Organism]", genome))
                        .await
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
                    let f = client
                        .fetch_by_accession_version(acc)
                        .await
                        .map_err(|e| format!("Failed to fetch {}: {}", acc, e))?;
                    let parts: Vec<&str> = acc.split('.').collect();
                    let version = if parts.len() >= 2 { parts[1] } else { "1" };
                    let base = if parts.len() >= 2 { parts[0] } else { acc };
                    cache
                        .cache_sequence(acc, version, base, &f)
                        .map_err(|e| e.to_string())?;
                    println!(
                        "({} bases)",
                        f.lines()
                            .filter(|l| !l.starts_with('>'))
                            .map(|l| l.len())
                            .sum::<usize>()
                            / 2
                    );
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
                        p.extension()
                            .map(|e| e == "fna" || e == "fasta" || e == "fa")
                            .unwrap_or(false)
                    })
                    .collect();
                if entries.is_empty() {
                    return Err(format!("No .fna/.fasta files found in {}", path.display()));
                }
                println!(
                    "  Found {} genome file(s) in {}",
                    entries.len(),
                    path.display()
                );
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
        if let Some(pattern) = &args.spaced_seed_pattern {
            bp.set_spaced_seed_pattern(pattern);
        }
        println!(
            "  Spaced seed: enabled (pattern: {})",
            bp.spaced_seed_pattern()
        );
        bp.set_spaced_seed(true);
    }

    if args.method != "none" {
        let fuzzy_method = match args.method.as_str() {
            "fuzzy-kmer" => FuzzyMethod::FuzzyKmer,
            "fuzzy-seed" => FuzzyMethod::FuzzySeed,
            "neighborhood" => FuzzyMethod::Neighborhood,
            _ => FuzzyMethod::None,
        };
        bp.set_fuzzy_method(fuzzy_method);
        bp.set_fuzzy_mismatches(args.fuzzy_mismatches);
        println!(
            "  Fuzzy method: {} (mismatches: {})",
            args.method, args.fuzzy_mismatches
        );
    }

    if args.top_n > 1 {
        bp.set_top_n(args.top_n);
    }

    bp.set_read_type(&args.read_type);
    println!("  Read type: {}", args.read_type);

    bp.set_search_radius(args.search_radius);
    println!("  Search radius: {}bp", args.search_radius);

    if let Some(chunk_size) = args.chunk_size {
        bp.set_chunk_size(chunk_size);
        println!("  Chunk size: {}bp (PacBio mode)", chunk_size);
    }

    if let Some(chunk_pct) = args.chunk_pct {
        bp.set_chunk_pct(chunk_pct);
    }

    if let Some(chunk_min) = args.chunk_min {
        bp.set_chunk_min(chunk_min);
    }

    if let Some(chunk_max) = args.chunk_max {
        bp.set_chunk_max(chunk_max);
    }

    if bp.chunk_pct() > 0.0 {
        println!(
            "  Chunk pct: {:.2}% (dynamic, clamped {}-{}bp)",
            bp.chunk_pct() * 100.0,
            bp.chunk_min(),
            bp.chunk_max()
        );
    }

    if let Some(threshold) = args.chunk_vote_threshold {
        bp.set_chunk_vote_threshold(threshold);
        println!(
            "  Chunk vote threshold: {:.0}% (requires {:.0}% chunk agreement)",
            threshold * 100.0,
            threshold * 100.0
        );
    }

    if let Some(top_n) = args.chunk_top_n {
        bp.set_chunk_top_n(top_n);
        println!("  Chunk top-N: {} genomes per read", top_n);
    }

    let run_chunk_strategy = match args.chunk_strategy.as_str() {
        "golden" => bit_pop::ChunkAnchorStrategy::Golden,
        "spaced" => bit_pop::ChunkAnchorStrategy::Spaced,
        _ => bit_pop::ChunkAnchorStrategy::Rarest,
    };
    bp.set_chunk_anchor_strategy(run_chunk_strategy);
    if args.chunk_strategy != "rarest" {
        println!("  Chunk strategy: {}", args.chunk_strategy);
    }

    let score_mode = match args.score_mode.as_str() {
        "base" => bit_pop::ChunkScoreMode::Base,
        _ => bit_pop::ChunkScoreMode::Quality,
    };
    bp.set_chunk_score_mode(score_mode);
    if args.score_mode != "quality" {
        println!("  Score mode: {}", args.score_mode);
    }

    bp.set_chunk_anchor_min_score(args.anchor_min_score);
    if args.anchor_min_score != 0.5 {
        println!("  Anchor min score: {}", args.anchor_min_score);
    }

    bp.set_chunk_min_score(args.min_score);
    if args.min_score > 0.0 {
        println!("  Chunk min score: {}", args.min_score);
    }

    if args.anchor_filter {
        bp.set_chunk_use_anchor_filter(true);
        println!("  Anchor filter: enabled (legacy mode)");
    }

    if args.snp_detect {
        bp.set_snp_detect(true);
        bp.set_snp_min_support(args.snp_min_support);
        println!(
            "  SNP detection: enabled (min support: {}, penalty: {})",
            args.snp_min_support,
            bp.snp_penalty()
        );
    }

    let run_align_mode = match args.align_mode.as_str() {
        "sw" => AlignMode::Sw,
        "hybrid" => AlignMode::Hybrid,
        "softclip" => AlignMode::Softclip,
        "chain" => AlignMode::Chain,
        _ => AlignMode::Xor,
    };
    bp.set_align_mode(run_align_mode);
    println!("  Alignment mode: {}", run_align_mode);

    // Step 3 (or 2): Map reads
    let step_map = if use_index { 2 } else { 3 };
    println!("\n[{}] Mapping reads...", step_map);

    let _mapped_count = if let (Some(r1_path), Some(r2_path)) = (&args.reads_1, &args.reads_2) {
        // Paired-end mode
        println!("  Paired-end mode");
        println!("    R1: {}", r1_path.display());
        println!("    R2: {}", r2_path.display());

        let pairs = bit_pop::fastq::parse_paired_fastq(
            r1_path.to_str().unwrap(),
            r2_path.to_str().unwrap(),
        )
        .map_err(|e| format!("Failed to parse paired FASTQ: {}", e))?;
        println!("  Loaded {} read pairs", pairs.len());

        let output_path = if let Some(ref p) = args.output {
            p.clone()
        } else {
            default_output_path(r1_path)
        };

        let total_pairs = pairs.len();
        let pb = ProgressBar::new(total_pairs as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner} Mapping pairs: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap(),
        );

        let mapped = if args.min_quality > 0 {
            let result = bp
                .map_paired_reads_parallel_quality(
                    &pairs,
                    output_path.to_str().unwrap(),
                    args.min_quality,
                    50,
                    true,
                    5,
                    args.bam,
                )
                .map_err(|e| format!("Mapping failed: {}", e))?;

            pb.set_position(total_pairs as u64);
            pb.set_message(format!("{} pairs", total_pairs));
            report_atomic_progress(total_pairs as u64, total_pairs as u64);
            pb.finish_with_message("Mapping complete");
            result
        } else {
            let result = bp
                .map_paired_reads_parallel(
                    &pairs,
                    output_path.to_str().unwrap(),
                    50,
                    true,
                    5,
                    args.bam,
                )
                .map_err(|e| format!("Mapping failed: {}", e))?;

            pb.set_position(total_pairs as u64);
            pb.set_message(format!("{} pairs", total_pairs));
            report_atomic_progress(total_pairs as u64, total_pairs as u64);
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
        let reads_path = args
            .reads
            .as_ref()
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

        let filtered_reads_fasta: Vec<(String, String)> =
            if args.min_quality > 0 && reads_format.has_quality() {
                if let ReadsFormat::Fastq(reads) = &reads_format {
                    let passed = bit_pop::fastq::filter_by_quality(reads, args.min_quality);
                    println!(
                        "  Quality filter (min Q{}): {}/{} reads passed",
                        args.min_quality,
                        passed.len(),
                        reads.len()
                    );
                    passed
                        .iter()
                        .map(|&i| (reads[i].0.clone(), reads[i].1.clone()))
                        .collect()
                } else {
                    println!("  Warning: quality filtering ignored for FASTA input");
                    reads_format
                        .iter_fasta()
                        .map(|(n, s)| (n.to_string(), s.to_string()))
                        .collect()
                }
            } else {
                reads_format
                    .iter_fasta()
                    .map(|(n, s)| (n.to_string(), s.to_string()))
                    .collect()
            };

        let reads_refs: Vec<(&str, &str)> = filtered_reads_fasta
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.as_str()))
            .collect();

        let total_reads = reads_refs.len();
        let mapped_count =
            if args.chunk_size.is_some() || args.chunk_pct.is_some() {
                // Chunk-based mapping for PacBio long reads
                let pb = ProgressBar::new(total_reads as u64);
                pb.set_style(ProgressStyle::default_bar()
                .template("{spinner} Chunk mapping: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap());

                let result = if args.threads > 1 {
                    bp.map_reads_with_chunking_parallel(
                        &reads_refs,
                        output_path.to_str().unwrap(),
                        50,
                    )
                    .map_err(|e| format!("Mapping failed: {}", e))?;
                    total_reads
                } else {
                    let pb_inner = pb.clone();
                    let result = bp
                        .map_reads_to_output_with_progress(
                            &reads_refs,
                            output_path.to_str().unwrap(),
                            50,
                            if total_reads > 1000 { 100 } else { 10 },
                            args.bam,
                            move |completed, total| {
                                pb_inner.set_position(completed as u64);
                                pb_inner.set_message(format!("{}/{} reads", completed, total));
                            },
                        )
                        .map_err(|e| format!("Mapping failed: {}", e))?;
                    result
                };

                pb.finish_with_message("Chunk mapping complete");
                result
            } else if args.threads > 1 {
                let pb = ProgressBar::new(total_reads as u64);
                pb.set_style(ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap());
                let pb_clone = pb.clone();

                let result = bp
                    .map_reads_parallel_with_progress(
                        &reads_refs,
                        output_path.to_str().unwrap(),
                        50,
                        if total_reads > 1000 { 100 } else { 10 },
                        move |completed, total| {
                            pb_clone.set_position(completed as u64);
                            pb_clone.set_message(format!("{}/{} reads", completed, total));
                            report_atomic_progress(completed as u64, total as u64);
                        },
                    )
                    .map_err(|e| format!("Mapping failed: {}", e))?;

                pb.finish_with_message("Mapping complete");
                result
            } else {
                let pb = ProgressBar::new(total_reads as u64);
                pb.set_style(ProgressStyle::default_bar()
                .template("{spinner} Mapping reads: [{elapsed_precise} {bar:40} {pos}/{len}] {msg}")
                .unwrap());
                let pb_clone = pb.clone();

                let result = bp
                    .map_reads_to_output_with_progress(
                        &reads_refs,
                        output_path.to_str().unwrap(),
                        50,
                        if total_reads > 1000 { 100 } else { 10 },
                        args.bam,
                        move |completed, total| {
                            pb_clone.set_position(completed as u64);
                            pb_clone.set_message(format!("{}/{} reads", completed, total));
                            report_atomic_progress(completed as u64, total as u64);
                        },
                    )
                    .map_err(|e| format!("Mapping failed: {}", e))?;

                pb.finish_with_message("Mapping complete");
                result
            };

        let elapsed = start.elapsed();

        // Parse SAM and show results
        println!("\n═══════════");
        println!("Done!");
        println!(
            "  Mapped:     {}/{} reads",
            mapped_count,
            filtered_reads_fasta.len()
        );
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

fn cmd_em(args: &EmArgs) {
    use bit_pop::em::{EMClassifier, EMConfig, ReadMappings};
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};
    use std::time::Instant;

    fn strip_read_suffix(name: &str) -> &str {
        // Match Python: rstrip("/1").rstrip("/2")
        // Python rstrip strips CHARACTERS not suffix, so:
        // "READ/1".rstrip("/1") → "READ" (strips "1" and "/")
        // "READ/2".rstrip("/1") → "READ/2" (no change)
        // "READ/2".rstrip("/2") → "READ"
        let name = name.strip_suffix("/1").unwrap_or(name);
        name.strip_suffix("/2").unwrap_or(name)
    }

    fn extract_nm_tag(parts: &[&str]) -> u32 {
        for field in parts.iter().skip(11) {
            if let Some(nm_str) = field.strip_prefix("NM:i:") {
                return nm_str.parse().unwrap_or(0);
            }
        }
        0
    }

    let start = Instant::now();

    println!("EM Classifier - Soft Assignment for Strain Resolution");
    println!("======================================================");
    println!();
    println!("Loading SAM: {}", args.input.to_string_lossy());

    let file = std::fs::File::open(&args.input).unwrap_or_else(|e| {
        eprintln!("Error opening SAM file: {}", e);
        std::process::exit(1);
    });
    let reader = BufReader::new(file);

    // Collect header lines and mappings
    let mut header_lines = Vec::new();
    let mut read_genomes: HashMap<String, Vec<(String, f64)>> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap_or_else(|e| {
            eprintln!("Error reading line: {}", e);
            std::process::exit(1);
        });

        if line.starts_with('@') {
            header_lines.push(line);
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 11 {
            continue;
        }

        let qname = strip_read_suffix(parts[0]).to_string();
        let flag: u16 = parts[1].parse().unwrap_or(0);
        let rname = parts[2];
        let mapq: u16 = if parts[4] != "0" {
            parts[4].parse().unwrap_or(0)
        } else {
            0
        };

        let is_supplementary = (flag & 0x800) != 0;
        let mut score = mapq as f64 / 60.0;

        if is_supplementary {
            score *= 0.5;
        }

        let nm = extract_nm_tag(&parts);
        let seq = parts.get(9).unwrap_or(&"");
        let read_len = seq.len() as f64;
        if read_len > 0.0 && nm > 0 {
            let mismatch_rate = nm as f64 / read_len;
            let nm_score = 1.0 - mismatch_rate;
            score = score * 0.7 + nm_score * 0.3;
        }

        if rname != "*" {
            read_genomes
                .entry(qname)
                .or_default()
                .push((rname.to_string(), score));
        }
    }

    let total_reads = read_genomes.len();
    let total_alignments: usize = read_genomes.values().map(|v| v.len()).sum();

    println!("  Loaded {} reads with mappings", total_reads);
    println!("  Total alignments: {}", total_alignments);

    // Build ReadMappings for EM
    let em_mappings: ReadMappings = read_genomes
        .iter()
        .flat_map(|(read_name, genomes)| {
            genomes
                .iter()
                .map(|(genome, score)| (read_name.clone(), genome.clone(), *score))
                .collect::<Vec<_>>()
        })
        .collect();

    // Run EM
    println!();
    println!("Running EM algorithm...");
    println!("  Convergence threshold: {}", args.convergence);
    println!("  Max iterations: {}", args.max_iterations);
    println!("  Temperature: {}", args.temperature);
    println!("  Top-K: {}", args.top_k);
    println!("  Confidence threshold: {}", args.confidence_threshold);
    println!();

    let em_start = Instant::now();

    let mut em = EMClassifier::new(EMConfig {
        convergence_threshold: args.convergence,
        max_iterations: args.max_iterations,
        temperature: args.temperature,
        top_k: args.top_k,
        confidence_threshold: args.confidence_threshold,
        ..EMConfig::default()
    });

    let hard_assignments = em.classify(&em_mappings);

    let em_time = em_start.elapsed();
    println!("EM completed in {:.2}s", em_time.as_secs_f64());
    println!("  Iterations: {}", em.iterations_run);
    println!("  Final KL divergence: {:.6}", em.final_kl);
    println!();

    // Create hard assignment lookup
    let hard_map: HashMap<String, Option<String>> = hard_assignments.into_iter().collect();

    // Write output SAM with EM improved classifications
    println!("Writing EM-improved SAM: {}", args.output.to_string_lossy());

    let mut out_file = std::fs::File::create(&args.output).unwrap_or_else(|e| {
        eprintln!("Error creating output file: {}", e);
        std::process::exit(1);
    });

    // Write headers
    for header in &header_lines {
        writeln!(out_file, "{}", header).unwrap();
    }

    // Re-read SAM and write with EM improvements
    let file = std::fs::File::open(&args.input).unwrap();
    let reader = BufReader::new(file);

    let mut changed = 0;
    let mut total = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('@') {
            writeln!(out_file, "{}", line).unwrap();
            continue;
        }

        let mut parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 11 {
            writeln!(out_file, "{}", line).unwrap();
            continue;
        }

        let qname = strip_read_suffix(parts[0]).to_string();
        let flag: u16 = parts[1].parse().unwrap_or(0);
        let current_rname = parts[2];
        let is_primary = (flag & 0x900) == 0;

        total += 1;

        // Check if this read has an EM assignment
        if let Some(Some(em_genome)) = hard_map.get(&qname) {
            if is_primary && current_rname != em_genome.as_str() {
                parts[2] = em_genome.as_str();
                parts[4] = "40";

                let new_line = parts.join("\t");
                writeln!(out_file, "{}", new_line).unwrap();
                changed += 1;
                continue;
            }
        }

        writeln!(out_file, "{}", line).unwrap();
    }

    let elapsed = start.elapsed();

    println!("  Total lines processed: {}", total);
    println!("  Classifications changed: {}", changed);
    println!();
    println!("EM completed in {:.2}s total", elapsed.as_secs_f64());
}

fn cmd_consensus(args: &ConsensusArgs) {
    use bit_pop::consensus::ConsensusStrategy;
    use std::time::Instant;

    let start = Instant::now();
    println!("Bit-Pop multi-k consensus");
    println!("==========================");
    println!();

    // Parse index paths (k-value will be read from index file)
    let mut index_paths: Vec<PathBuf> = Vec::new();
    for idx_arg in &args.indexes {
        let path = PathBuf::from(idx_arg);
        if !path.exists() {
            eprintln!("Error: index file not found: {}", path.display());
            std::process::exit(1);
        }
        index_paths.push(path);
    }

    println!("Indexes:");
    for path in &index_paths {
        println!("  {}", path.display());
    }
    println!();

    // Parse strategy
    let strategy = match args.strategy.as_str() {
        "majority" => ConsensusStrategy::Majority,
        "best_score" => ConsensusStrategy::BestScore,
        "base_score" => ConsensusStrategy::BaseScore,
        _ => ConsensusStrategy::WeightedScore,
    };

    println!("[1/3] Loading indexes...");
    let mut consensus =
        MultiKConsensus::from_paths(&index_paths, args.min_score).unwrap_or_else(|e| {
            eprintln!("Error loading indexes: {}", e);
            std::process::exit(1);
        });

    consensus.strategy = strategy;
    // Set top_n on each BitPop index (controls rare k-mer anchors)
    for bp in consensus.indexes.values_mut() {
        bp.set_top_n(args.top_n);
        bp.set_chunk_anchor_min_score(args.anchor_min_score);
        bp.set_chunk_min_score(args.min_score);
        bp.set_chunk_use_anchor_filter(args.anchor_filter);
    }
    consensus.chunk_size = args.chunk_size;
    consensus.chunk_pct = args.chunk_pct;
    consensus.chunk_min = args.chunk_min;
    consensus.chunk_max = args.chunk_max;
    consensus.enable_snp_detect = args.snp_detect;
    consensus.snp_min_support = args.snp_min_support;
    consensus.snp_penalty = args.snp_penalty;
    consensus.min_k_mappings = args.min_k_mappings;
    consensus.top_n = args.top_n;

    println!();
    println!("[2/3] Mapping reads...");
    println!("  Reads: {}", args.reads.display());
    println!("  Output: {}", args.output.display());
    println!(
        "  Strategy: {}",
        match strategy {
            ConsensusStrategy::Majority => "majority",
            ConsensusStrategy::WeightedScore => "weighted_score",
            ConsensusStrategy::BestScore => "best_score",
            ConsensusStrategy::BaseScore => "base_score",
        }
    );
    println!("  Threads: {}", args.threads);
    println!(
        "  Chunk size: {} (fixed), {}% (dynamic)",
        args.chunk_size, args.chunk_pct
    );

    // Two-pass mode: map each k separately, then combine (like Python script)
    if args.two_pass {
        println!("  Two-pass: enabled (map each k separately, then combine)");
        match consensus.map_reads_to_sam_two_pass(&args.reads, &args.output, args.threads, false) {
            Ok((mapped, total)) => {
                println!();
                println!("==========================");
                println!("Done!");
                println!("  Mapped: {} / {} reads", mapped, total);
                println!("  Output: {}", args.output.display());
                println!();
                println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
            }
            Err(e) => {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    // Streaming mode
    if args.stream {
        let chunk_size = parse_stream_chunk_size(&args.max_ram);
        println!("  Streaming: enabled (chunk={} reads)", chunk_size);
        match consensus.map_reads_to_sam_stream(&args.reads, &args.output, args.threads, chunk_size)
        {
            Ok((mapped, total)) => {
                println!();
                println!("==========================");
                println!("Done!");
                println!("  Mapped: {} / {} reads", mapped, total);
                println!("  Output: {}", args.output.display());
                println!();
                println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
            }
            Err(e) => {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    match consensus.map_reads_to_sam(&args.reads, &args.output, args.threads) {
        Ok((mapped, total)) => {
            println!();
            println!("==========================");
            println!("Done!");
            println!("  Mapped: {} / {} reads", mapped, total);
            println!("  Output: {}", args.output.display());
            println!();
            println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_concon(args: &ConConArgs) {
    use std::time::Instant;

    let start = Instant::now();
    println!("Bit-Pop consensus (subprocess)");
    println!("=====================================");
    println!();

    // Parse index paths (k-value will be read from index file)
    let mut index_paths: Vec<PathBuf> = Vec::new();
    for idx_arg in &args.indexes {
        let path = PathBuf::from(idx_arg);
        if !path.exists() {
            eprintln!("Error: index file not found: {}", path.display());
            std::process::exit(1);
        }
        index_paths.push(path);
    }

    println!("Indexes:");
    for path in &index_paths {
        println!("  {}", path.display());
    }
    println!();

    // Auto-detect bit-pop executable
    let bit_pop_exe = if let Some(ref p) = args.bit_pop {
        p.clone()
    } else {
        // Try to find the current executable's directory
        let current_exe = std::env::current_exe().unwrap_or_default();
        let parent = current_exe.parent().unwrap_or(Path::new("."));
        let candidate = parent.join("bit-pop.exe");
        if candidate.exists() {
            candidate
        } else {
            // Fallback: look in target/release
            let script_dir = std::path::PathBuf::from("target/release");
            let candidate = script_dir.join("bit-pop.exe");
            if candidate.exists() {
                candidate
            } else {
                eprintln!("Error: bit-pop.exe not found. Use --bit-pop to specify path.");
                std::process::exit(1);
            }
        }
    };

    println!("bit-pop executable: {}", bit_pop_exe.display());
    println!();

    // Parse strategy
    let strategy = match args.strategy.as_str() {
        "majority" => bit_pop::concon::ConsensusStrategy::Majority,
        "best_score" => bit_pop::concon::ConsensusStrategy::BestScore,
        "base_score" => bit_pop::concon::ConsensusStrategy::BaseScore,
        _ => bit_pop::concon::ConsensusStrategy::WeightedScore,
    };

    let mut concon = ConCon::new(
        index_paths,
        bit_pop_exe,
        args.min_score,
        args.context_window,
        args.anchor_min_score,
        args.anchor_filter,
    )
    .unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });

    concon.strategy = strategy;
    concon.min_k_mappings = args.min_k_mappings;
    concon.map_top_n = args.top_n;
    concon.top_n = args.consensus_top_n;
    concon.chunk_size = args.chunk_size;
    concon.chunk_pct = args.chunk_pct;
    concon.chunk_min = args.chunk_min;
    concon.chunk_max = args.chunk_max;

    println!("[1/3] Configuration");
    println!("  Reads: {}", args.reads.display());
    println!("  Output: {}", args.output.display());
    println!(
        "  Strategy: {}",
        match strategy {
            bit_pop::concon::ConsensusStrategy::Majority => "majority",
            bit_pop::concon::ConsensusStrategy::WeightedScore => "weighted_score",
            bit_pop::concon::ConsensusStrategy::BestScore => "best_score",
            bit_pop::concon::ConsensusStrategy::BaseScore => "base_score",
        }
    );
    println!("  Threads per map: {}", args.threads);
    if args.chunk_pct > 0.0 {
        println!(
            "  Chunk pct: {:.2}% (dynamic, clamped {}-{}bp)",
            args.chunk_pct * 100.0,
            args.chunk_min,
            args.chunk_max
        );
    } else if args.chunk_size > 0 {
        println!("  Chunk size: {}bp (fixed)", args.chunk_size);
    }
    println!();

    println!("[2/3] Running consensus...");
    match concon.run(&args.reads, &args.output, args.threads) {
        Ok((mapped, total)) => {
            println!();
            println!("=====================================");
            println!("Done!");
            println!("  Mapped: {} / {} reads", mapped, total);
            println!("  Output: {}", args.output.display());
            println!();
            println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_chunk_consensus(args: &ChunkConsensusArgs) {
    use bit_pop::chunk_consensus::ConsensusStrategy;
    use std::time::Instant;

    let start = Instant::now();
    println!("Bit-Pop chunk-consensus");
    println!("========================");
    println!();

    // Parse chunk percentages
    let chunk_pcts: Vec<f64> = args
        .chunk_pcts
        .split(',')
        .map(|s| {
            s.trim().parse::<f64>().unwrap_or_else(|e| {
                eprintln!("Error: invalid chunk_pct '{}': {}", s.trim(), e);
                std::process::exit(1);
            })
        })
        .collect();

    if chunk_pcts.is_empty() {
        eprintln!("Error: at least one chunk_pct is required");
        std::process::exit(1);
    }

    let index_path = args.index.to_str().unwrap();
    if !args.index.exists() {
        eprintln!("Error: index file not found: {}", args.index.display());
        std::process::exit(1);
    }

    println!("Index: {}", args.index.display());
    println!(
        "Chunk configs: {:?}",
        chunk_pcts
            .iter()
            .map(|p| format!("{:.0}%", p * 100.0))
            .collect::<Vec<_>>()
    );
    println!();

    let strategy = match args.strategy.as_str() {
        "majority" => ConsensusStrategy::Majority,
        "base_score" => ConsensusStrategy::BaseScore,
        _ => ConsensusStrategy::WeightedScore,
    };

    println!("[1/3] Loading index ({} configs)...", chunk_pcts.len());
    let mut consensus = MultiChunkConsensus::from_path(
        index_path,
        &chunk_pcts,
        args.chunk_min,
        args.chunk_max,
        args.min_score,
        args.anchor_min_score,
    )
    .unwrap_or_else(|e| {
        eprintln!("Error loading index: {}", e);
        std::process::exit(1);
    });

    consensus.strategy = strategy;
    if let Some(min_agree) = args.min_agreement {
        consensus.min_agreement = min_agree;
    }
    consensus.top_n = args.top_n;

    println!();
    println!("[2/3] Mapping reads...");
    println!("  Reads: {}", args.reads.display());
    println!("  Output: {}", args.output.display());
    println!(
        "  Strategy: {}",
        if strategy == ConsensusStrategy::Majority {
            "majority"
        } else {
            "weighted_score"
        }
    );
    println!(
        "  Min agreement: {}/{}",
        consensus.min_agreement,
        chunk_pcts.len()
    );
    println!("  Threads: {}", args.threads);

    let result = if args.stream {
        let chunk_size = if let Some(ref max_ram) = args.max_ram {
            parse_stream_chunk_size(&Some(max_ram.clone()))
        } else {
            20_000_000
        };
        println!("  Streaming: chunk size = {} reads", chunk_size);
        consensus.map_reads_to_sam_stream(&args.reads, &args.output, args.threads, chunk_size)
    } else {
        consensus.map_reads_to_sam(&args.reads, &args.output, args.threads)
    };

    match result {
        Ok((mapped, total)) => {
            println!();
            println!("========================");
            println!("Done!");
            println!(
                "  Mapped: {} / {} reads ({:.1}%)",
                mapped,
                total,
                mapped as f64 / total as f64 * 100.0
            );
            println!("  Output: {}", args.output.display());
            println!();
            println!("Total time: {:.2}s", start.elapsed().as_secs_f64());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_tax(args: &TaxArgs) {
    use bit_pop::taxonomy::{compute_tax_report, format_tax_report, TaxonomyTree};
    use std::collections::HashMap;
    use std::io::BufRead;
    use std::time::Instant;

    let start = Instant::now();

    println!("Bit-Pop Taxonomic Classification Report");
    println!("========================================");
    println!();

    // Load taxonomy tree
    println!("Loading taxonomy tree...");
    let tax_start = Instant::now();
    let mut tree = match TaxonomyTree::load(
        args.nodes_dmp.to_str().unwrap(),
        args.names_dmp.to_str().unwrap(),
    ) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error loading taxonomy: {}", e);
            eprintln!();
            eprintln!("Download from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz");
            eprintln!(
                "Extract nodes.dmp and names.dmp, then pass with --nodes-dmp and --names-dmp"
            );
            std::process::exit(1);
        }
    };
    let tax_time = tax_start.elapsed();
    println!(
        "  Loaded {} taxonomy nodes in {:.2}s",
        tree.node_count(),
        tax_time.as_secs_f64()
    );

    // Parse SAM file
    println!("Loading SAM: {}", args.input.to_string_lossy());
    let file = match std::fs::File::open(&args.input) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening SAM file: {}", e);
            std::process::exit(1);
        }
    };
    let reader = std::io::BufReader::new(file);

    let mut genome_counts: HashMap<String, usize> = HashMap::new();
    let mut seen_reads: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut unmapped: usize = 0;

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };

        if line.starts_with('@') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let read_name = fields[0];
        let ref_name = fields[2];

        if ref_name == "*" {
            unmapped += 1;
            continue;
        }

        // Count each read only once (primary mapping)
        if seen_reads.insert(read_name.to_string()) {
            *genome_counts.entry(ref_name.to_string()).or_insert(0) += 1;
        }
    }

    let total_mapped: usize = genome_counts.values().sum();
    println!(
        "  Loaded {} mapped reads (+ {} unmapped)",
        total_mapped, unmapped
    );
    println!("  Found {} unique genomes", genome_counts.len());

    // Map genome names to taxonomy IDs
    println!("Mapping genomes to taxonomy...");
    let mut mapped = 0usize;
    let mut unmapped_genomes = Vec::new();

    for genome_name in genome_counts.keys() {
        if tree.map_genome_name(genome_name).is_some() {
            mapped += 1;
        } else {
            unmapped_genomes.push(genome_name.clone());
        }
    }

    println!(
        "  Mapped {}/{} genomes to taxonomy",
        mapped,
        genome_counts.len()
    );

    if !unmapped_genomes.is_empty() {
        println!(
            "  Warning: {} genomes not found in taxonomy:",
            unmapped_genomes.len()
        );
        for name in &unmapped_genomes[..unmapped_genomes.len().min(5)] {
            println!("    - {}", name);
        }
        if unmapped_genomes.len() > 5 {
            println!("    ... and {} more", unmapped_genomes.len() - 5);
        }
    }

    // Compute taxonomic report
    println!();
    let report = compute_tax_report(&tree, &genome_counts);

    // Output report
    let report_text = format_tax_report(&report, args.top_n);

    if args.format == "json" {
        let mut json_entries: Vec<serde_json::Value> = Vec::new();
        for (rank, entries) in &report.ranks {
            if rank == "no rank" {
                continue;
            }
            for (name, count, pct) in entries {
                json_entries.push(serde_json::json!({
                    "rank": rank,
                    "name": name,
                    "count": count,
                    "percentage": (pct * 100.0).round() / 100.0,
                }));
            }
        }
        let json_output = serde_json::json!({
            "total_reads": report.total_reads,
            "taxonomic_breakdown": json_entries,
        });
        print!("{}", serde_json::to_string_pretty(&json_output).unwrap());
    } else {
        print!("{}", report_text);
    }

    // Write to file if specified
    if let Some(ref output_path) = args.output {
        let content = if args.format == "json" {
            let mut json_entries: Vec<serde_json::Value> = Vec::new();
            for (rank, entries) in &report.ranks {
                if rank == "no rank" {
                    continue;
                }
                for (name, count, pct) in entries {
                    json_entries.push(serde_json::json!({
                        "rank": rank,
                        "name": name,
                        "count": count,
                        "percentage": (pct * 100.0).round() / 100.0,
                    }));
                }
            }
            let json_output = serde_json::json!({
                "total_reads": report.total_reads,
                "taxonomic_breakdown": json_entries,
            });
            serde_json::to_string_pretty(&json_output).unwrap()
        } else {
            report_text.clone()
        };
        match std::fs::write(output_path, content) {
            Ok(_) => println!("\nReport saved to: {}", output_path.display()),
            Err(e) => eprintln!("Error writing report: {}", e),
        }
    }

    let elapsed = start.elapsed();
    println!("\nTotal time: {:.2}s", elapsed.as_secs_f64());
}

/// Parse max-ram string and calculate optimal chunk size.
/// Returns chunk size in number of reads.
fn parse_stream_chunk_size(max_ram: &Option<String>) -> usize {
    let total_ram_bytes = if let Some(ref ram_str) = max_ram {
        parse_ram_bytes(ram_str)
    } else {
        // Default: use available system memory (estimate 64GB if unknown)
        64 * 1024 * 1024 * 1024
    };

    // Reserve 10% for OS and indexes
    let available = total_ram_bytes as f64 * 0.9;

    // Estimate ~500 bytes per read in memory (name + seq + qual + overhead)
    let bytes_per_read = 500.0;
    let chunk_size = (available / bytes_per_read) as usize;

    // Clamp to reasonable range
    chunk_size.clamp(500_000, 20_000_000)
}

/// Parse RAM string like "32G", "16GB", "8G" to bytes.
fn parse_ram_bytes(s: &str) -> u64 {
    let s = s.trim().to_uppercase();
    let (num_str, unit) = if s.ends_with("GB") {
        (&s[..s.len() - 2], "GB")
    } else if s.ends_with('G') {
        (&s[..s.len() - 1], "G")
    } else if s.ends_with("MB") {
        (&s[..s.len() - 2], "MB")
    } else if s.ends_with('M') {
        (&s[..s.len() - 1], "M")
    } else {
        (s.as_str(), "B")
    };

    let num: f64 = num_str.trim().parse().unwrap_or(64.0);
    match unit {
        "GB" | "G" => (num * 1024.0 * 1024.0 * 1024.0) as u64,
        "MB" | "M" => (num * 1024.0 * 1024.0) as u64,
        _ => num as u64,
    }
}
