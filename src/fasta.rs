use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Reads FASTA files and yields (header, sequence) pairs.
/// Streaming: does not load entire file into memory at once.
pub struct FastaReader {
    reader: BufReader<File>,
    buffer: String,
    /// Buffer for a lookahead line (next header found while reading sequences)
    lookahead: Option<String>,
}

impl FastaReader {
    /// Create a new FastaReader from a file path.
    pub fn new(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
            buffer: String::new(),
            lookahead: None,
        })
    }

    /// Read the next (header, sequence) pair from the FASTA file.
    /// Header is everything after '>' on the header line (trimmed).
    /// Sequence is all subsequent non-header lines joined together (uppercased).
    pub fn next(&mut self) -> Option<io::Result<(String, String)>> {
        // Use lookahead if available
        let first_line = if let Some(line) = self.lookahead.take() {
            line
        } else {
            loop {
                self.buffer.clear();
                match self.reader.read_line(&mut self.buffer) {
                    Ok(0) => return None,
                    Ok(_) => {},
                    Err(e) => return Some(Err(e)),
                }

                let line = self.buffer.trim().to_string();
                if line.is_empty() {
                    continue;
                }

                if !line.starts_with('>') {
                    continue;
                }

                break line;
            }
        };

        let header = first_line[1..].trim().to_string();
        let mut sequence = String::new();

        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => break,
                Ok(_) => {},
                Err(e) => return Some(Err(e)),
            }

            let seq_line = self.buffer.trim();
            if seq_line.is_empty() {
                continue;
            }

            if seq_line.starts_with('>') {
                self.lookahead = Some(seq_line.to_string());
                break;
            }

            sequence.push_str(seq_line.trim().to_uppercase().as_str());
        }

        Some(Ok((header, sequence)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    struct TempFile {
        path: std::path::PathBuf,
    }

    impl TempFile {
        fn new(content: &str) -> Self {
            let dir = std::env::temp_dir();
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let path = dir.join(format!("bitpop_test_{}_{}.fasta", std::process::id(), id));
            let mut f = File::create(&path).expect("Create temp file");
            f.write_all(content.as_bytes()).expect("Write to temp file");
            Self { path }
        }
    }

    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.path);
        }
    }

    #[test]
    fn test_single_fasta() {
        let file = TempFile::new(">seq1 test sequence\nACGTACGT\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let result = reader.next().unwrap().unwrap();
        assert_eq!(result.0, "seq1 test sequence");
        assert_eq!(result.1, "ACGTACGT");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_multi_fasta() {
        let file = TempFile::new(">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (h1, s1) = reader.next().unwrap().unwrap();
        assert_eq!(h1, "chr1");
        assert_eq!(s1, "ACGTACGT");

        let (h2, s2) = reader.next().unwrap().unwrap();
        assert_eq!(h2, "chr2");
        assert_eq!(s2, "TTTTGGGG");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_multiline_sequence() {
        let file = TempFile::new(">seq1\nACGT\nACGT\nACGT\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (_, seq) = reader.next().unwrap().unwrap();
        assert_eq!(seq, "ACGTACGTACGT");
    }

    #[test]
    fn test_lowercase_converted() {
        let file = TempFile::new(">seq1\nacgtacgt\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (_, seq) = reader.next().unwrap().unwrap();
        assert_eq!(seq, "ACGTACGT");
    }

    #[test]
    fn test_empty_file() {
        let file = TempFile::new("");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_header_only() {
        let file = TempFile::new(">seq1\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (header, seq) = reader.next().unwrap().unwrap();
        assert_eq!(header, "seq1");
        assert_eq!(seq, "");
    }

    #[test]
    fn test_whitespace_in_header() {
        let file = TempFile::new(">  seq1  description  \nACGT\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (header, _) = reader.next().unwrap().unwrap();
        assert_eq!(header, "seq1  description");
    }

    #[test]
    fn test_nonexistent_file() {
        let result = FastaReader::new("/nonexistent/path/file.fasta");
        assert!(result.is_err());
    }

    #[test]
    fn test_blank_lines_skipped() {
        let file = TempFile::new("\n\n>seq1\n\nACGT\n\nACGT\n\n");
        let mut reader = FastaReader::new(file.path.to_str().unwrap()).unwrap();

        let (_, seq) = reader.next().unwrap().unwrap();
        assert_eq!(seq, "ACGTACGT");
    }
}
