use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const CACHE_DIR_NAME: &str = ".bitpop";
const MANIFEST_FILE: &str = "manifest.json";
const SEQUENCES_DIR: &str = "sequences";
const INDEXES_DIR: &str = "indexes";

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CachedGenome {
    pub accession: String,
    pub version: String,
    pub base_accession: String,
    pub download_date: String,
    pub checksum: String,
    pub fasta_path: String,
    pub index_path: Option<String>,
    pub genome_size: usize,
    pub kmer_size: Option<usize>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CacheManifest {
    pub version: u32,
    pub genomes: HashMap<String, CachedGenome>,
}

impl Default for CacheManifest {
    fn default() -> Self {
        Self::new()
    }
}

impl CacheManifest {
    pub fn new() -> Self {
        Self {
            version: 1,
            genomes: HashMap::new(),
        }
    }

    pub fn from_file(path: &Path) -> io::Result<Self> {
        let content = fs::read_to_string(path)?;
        let manifest: CacheManifest = serde_json::from_str(&content)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(manifest)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let content = serde_json::to_string_pretty(self)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        fs::write(path, content)
    }

    pub fn get(&self, accession: &str) -> Option<&CachedGenome> {
        self.genomes.get(accession)
    }

    pub fn insert(&mut self, accession: String, genome: CachedGenome) {
        self.genomes.insert(accession, genome);
    }

    pub fn remove(&mut self, accession: &str) -> bool {
        self.genomes.remove(accession).is_some()
    }

    pub fn len(&self) -> usize {
        self.genomes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.genomes.is_empty()
    }

    pub fn compute_checksum(path: &Path) -> io::Result<String> {
        let content = fs::read(path)?;
        let mut hasher = Sha256::new();
        hasher.update(&content);
        Ok(format!("{:x}", hasher.finalize()))
    }

    pub fn needs_update(&self, accession: &str, current_checksum: &str) -> bool {
        match self.genomes.get(accession) {
            Some(cached) => cached.checksum != current_checksum,
            None => true,
        }
    }
}

pub struct CacheManager {
    cache_dir: PathBuf,
    manifest_path: PathBuf,
    manifest: CacheManifest,
}

impl CacheManager {
    pub fn new(custom_dir: Option<PathBuf>) -> io::Result<Self> {
        let cache_dir = match custom_dir {
            Some(dir) => dir,
            None => {
                let home = dirs::cache_dir()
                    .or_else(|| std::env::current_dir().ok())
                    .unwrap_or_else(|| PathBuf::from("."));
                home.join(CACHE_DIR_NAME)
            }
        };

        let sequences_dir = cache_dir.join(SEQUENCES_DIR);
        let indexes_dir = cache_dir.join(INDEXES_DIR);

        fs::create_dir_all(&sequences_dir)?;
        fs::create_dir_all(&indexes_dir)?;

        let manifest_path = cache_dir.join(MANIFEST_FILE);
        let manifest = if manifest_path.exists() {
            CacheManifest::from_file(&manifest_path).unwrap_or_else(|_| CacheManifest::new())
        } else {
            CacheManifest::new()
        };

        Ok(Self {
            cache_dir,
            manifest_path,
            manifest,
        })
    }

    pub fn cache_dir(&self) -> &Path {
        &self.cache_dir
    }

    pub fn sequences_dir(&self) -> PathBuf {
        self.cache_dir.join(SEQUENCES_DIR)
    }

    pub fn indexes_dir(&self) -> PathBuf {
        self.cache_dir.join(INDEXES_DIR)
    }

    pub fn manifest(&self) -> &CacheManifest {
        &self.manifest
    }

    pub fn manifest_mut(&mut self) -> &mut CacheManifest {
        &mut self.manifest
    }

    pub fn save_manifest(&self) -> io::Result<()> {
        self.manifest.save(&self.manifest_path)
    }

    pub fn get_fasta_path(&self, accession: &str) -> PathBuf {
        let filename = format!("{}.fasta", accession);
        self.cache_dir.join(SEQUENCES_DIR).join(filename)
    }

    pub fn get_index_path(&self, accession: &str, k: usize) -> PathBuf {
        let filename = format!("{}_k{}.bitpop", accession, k);
        self.cache_dir.join(INDEXES_DIR).join(filename)
    }

    pub fn has_sequence(&self, accession: &str) -> bool {
        let genome = self.manifest.get(accession);
        match genome {
            Some(g) => {
                let path = self.get_fasta_path(&g.accession);
                path.exists()
            }
            None => false,
        }
    }

    pub fn has_index(&self, accession: &str, k: usize) -> bool {
        let genome = self.manifest.get(accession);
        match genome {
            Some(g) => match &g.index_path {
                Some(path) => Path::new(path).exists(),
                None => {
                    let path = self.get_index_path(accession, k);
                    path.exists()
                }
            },
            None => false,
        }
    }

    pub fn cache_sequence(
        &mut self,
        accession: &str,
        version: &str,
        base_accession: &str,
        fasta_content: &str,
    ) -> io::Result<()> {
        let fasta_path = self.get_fasta_path(accession);
        fs::write(&fasta_path, fasta_content)?;

        let checksum = CacheManifest::compute_checksum(&fasta_path)?;
        let genome_size = fasta_content.len();

        let genome = CachedGenome {
            accession: accession.to_string(),
            version: version.to_string(),
            base_accession: base_accession.to_string(),
            download_date: chrono::Utc::now().to_rfc3339(),
            checksum,
            fasta_path: fasta_path.to_string_lossy().to_string(),
            index_path: None,
            genome_size,
            kmer_size: None,
        };

        self.manifest.insert(accession.to_string(), genome);
        self.save_manifest()?;

        Ok(())
    }

    pub fn cache_index(&mut self, accession: &str, index_path: &Path, k: usize) -> io::Result<()> {
        let dest = self.get_index_path(accession, k);
        fs::copy(index_path, &dest)?;

        if let Some(genome) = self.manifest.genomes.get_mut(accession) {
            genome.index_path = Some(dest.to_string_lossy().to_string());
            genome.kmer_size = Some(k);
        }

        self.save_manifest()?;
        Ok(())
    }

    pub fn remove_genome(&mut self, accession: &str) -> io::Result<()> {
        if let Some(genome) = self.manifest.genomes.get(accession) {
            let fasta = Path::new(&genome.fasta_path);
            if fasta.exists() {
                fs::remove_file(fasta)?;
            }
            if let Some(ref idx) = genome.index_path {
                let index = Path::new(idx);
                if index.exists() {
                    fs::remove_file(index)?;
                }
            }
        }
        self.manifest.remove(accession);
        self.save_manifest()?;
        Ok(())
    }

    pub fn list_genomes(&self) -> Vec<&CachedGenome> {
        self.manifest.genomes.values().collect()
    }

    pub fn get_or_download<F>(
        &mut self,
        accession: &str,
        _k: usize,
        fetch_fn: F,
    ) -> io::Result<Option<PathBuf>>
    where
        F: FnOnce(&str) -> Result<String, String>,
    {
        if self.has_sequence(accession) {
            let genome = self.manifest.get(accession).unwrap();
            return Ok(Some(Path::new(&genome.fasta_path).to_path_buf()));
        }

        match fetch_fn(accession) {
            Ok(fasta) => {
                let parts: Vec<&str> = accession.split('.').collect();
                let version = if parts.len() >= 2 { parts[1] } else { "1" };
                let base = if parts.len() >= 2 {
                    parts[0]
                } else {
                    accession
                };
                self.cache_sequence(accession, version, base, &fasta)?;
                let genome = self.manifest.get(accession).unwrap();
                Ok(Some(Path::new(&genome.fasta_path).to_path_buf()))
            }
            Err(e) => {
                eprintln!("  Error downloading {}: {}", accession, e);
                Ok(None)
            }
        }
    }
}
