use bit_pop::persisted::load_bitpop;
use bit_pop::decode_base;
use std::io::Write;

/// Extract genome sequences from a bitpop index file and write to FASTA
fn main() {
    let args: Vec<String> = std::env::args().collect();
    
    if args.len() < 3 {
        eprintln!("Usage: extract_seqs <index.bitpop> <output.fasta>");
        std::process::exit(1);
    }
    
    let index_path = &args[1];
    let output_path = &args[2];
    
    println!("Loading index from: {}", index_path);
    let bp = load_bitpop(index_path).expect("Failed to load index");
    
    println!("Index loaded: {} genomes, k={}", bp.genome_count(), bp.k());
    
    let names = bp.genome_names_ordered();
    println!("Writing {} sequences to: {}", names.len(), output_path);
    
    let file = std::fs::File::create(output_path).expect("Failed to create output file");
    let mut writer = std::io::BufWriter::new(file);
    
    for (gid, name) in names.iter().enumerate() {
        if let Some(seq_bytes) = bp.get_genome_seq(gid as u32) {
            let seq: String = seq_bytes.iter().map(|&b| decode_base(b)).collect();
            writeln!(writer, ">{}", name).expect("Failed to write header");
            for i in (0..seq.len()).step_by(80) {
                let end = (i + 80).min(seq.len());
                writeln!(writer, "{}", &seq[i..end]).expect("Failed to write");
            }
            println!("  Written {}: {} bases", name, seq.len());
        }
    }
    
    writer.flush().expect("Failed to flush");
    println!("Done!");
}
