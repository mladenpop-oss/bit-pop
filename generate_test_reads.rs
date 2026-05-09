use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        SimpleRng { state: seed }
    }
    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.state
    }
    fn next_range(&mut self, max: u64) -> u64 {
        self.next() % max
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <fasta_path> <fastq_output> [num_reads_per_genome]", args[0]);
        std::process::exit(1);
    }

    let fasta_path = Path::new(&args[1]);
    let fastq_out = Path::new(&args[2]);
    let num_reads_per_genome: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    
    let file = File::open(fasta_path).expect("Cannot open FASTA");
    let reader = BufReader::new(file);
    
    let mut headers: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();
    
    for line in reader.lines() {
        let line = line.expect("Cannot read line");
        let line = line.trim().replace("\r", "");
        if line.starts_with('>') {
            if !current_header.is_empty() {
                headers.push(current_header.clone());
                sequences.push(current_seq.clone());
            }
            current_header = line[1..].to_string();
            current_seq = String::new();
        } else if !line.is_empty() {
            current_seq.push_str(&line);
        }
    }
    if !current_header.is_empty() {
        headers.push(current_header);
        sequences.push(current_seq);
    }
    
    println!("Loaded {} genomes, {} total bases", headers.len(), sequences.iter().map(|s| s.len()).sum::<usize>());
    
    let mut out = File::create(fastq_out).expect("Cannot create FASTQ");
    let mut rng = SimpleRng::new(42);
    let mut read_id = 0u32;
    let bases = [b'A', b'C', b'G', b'T'];
    
    for (i, seq) in sequences.iter().enumerate() {
        if seq.len() < 100 { continue; }
        
        let num_reads = if num_reads_per_genome > 0 {
            num_reads_per_genome
        } else {
            3 + (rng.next_range(3) as usize)
        };
        for _ in 0..num_reads {
            if seq.len() < 75 { continue; }
            let max_start = seq.len() - 75 + 1;
            let start = rng.next_range(max_start as u64) as usize;
            let read_len = 75 + rng.next_range(51) as usize;
            let end = (start + read_len).min(seq.len());
            let read = &seq.as_bytes()[start..end];
            
            let mut error_read = Vec::new();
            for (j, &c) in read.iter().enumerate() {
                if rng.next_range(50) == 0 && j > 2 && j < read.len() - 2 {
                    error_read.push(bases[rng.next_range(4) as usize]);
                } else {
                    error_read.push(c);
                }
            }
            
            let quality: String = (0..error_read.len()).map(|_| {
                if rng.next_range(10) == 0 { '!' } else { 'I' }
            }).collect();
            
            let header = format!("@read_{}_genome_{} pos:{} len:{}", read_id, i, start, error_read.len());
            let fastq = format!("{}\n{}\n+\n{}\n", header, String::from_utf8_lossy(&error_read), quality);
            out.write_all(fastq.as_bytes()).expect("Write failed");
            read_id += 1;
        }
    }
    
    println!("Generated {} reads -> {}", read_id, fastq_out.display());
}
