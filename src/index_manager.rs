use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use crate::cache::{CacheManager, CacheManifest};
use crate::ncbi::NcbiClient;
use crate::BitPop;

#[derive(Debug, Serialize, Deserialize)]
pub struct IndexMetadata {
    pub version: u32,
    pub created_at: String,
    pub updated_at: String,
    pub kmer_size: usize,
    pub genomes: HashMap<u32, GenomeEntry>,
    pub genome_count: u32,
    pub total_bases: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct GenomeEntry {
    pub name: String,
    pub accession: Option<String>,
    pub local_path: Option<String>,
    pub sequence_length: usize,
    pub added_at: String,
    pub checksum: String,
}

impl IndexMetadata {
    pub fn new(kmer_size: usize) -> Self {
        Self {
            version: 1,
            created_at: chrono::Utc::now().to_rfc3339(),
            updated_at: chrono::Utc::now().to_rfc3339(),
            kmer_size,
            genomes: HashMap::new(),
            genome_count: 0,
            total_bases: 0,
        }
    }

    pub fn from_file(path: &Path) -> io::Result<Self> {
        let content = fs::read_to_string(path)?;
        let meta: IndexMetadata = serde_json::from_str(&content)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(meta)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let content = serde_json::to_string_pretty(self)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        fs::write(path, content)
    }

    pub fn add_genome(&mut self, name: &str, accession: Option<String>, local_path: Option<String>, sequence_length: usize, checksum: String) {
        let gid = self.genome_count;
        self.genomes.insert(gid, GenomeEntry {
            name: name.to_string(),
            accession,
            local_path,
            sequence_length,
            added_at: chrono::Utc::now().to_rfc3339(),
            checksum,
        });
        self.genome_count += 1;
        self.total_bases += sequence_length;
        self.updated_at = chrono::Utc::now().to_rfc3339();
    }

    pub fn get_genome(&self, gid: u32) -> Option<&GenomeEntry> {
        self.genomes.get(&gid)
    }
}

pub struct DynamicIndexManager {
    metadata: IndexMetadata,
    cache: CacheManager,
    index_path: PathBuf,
    metadata_path: PathBuf,
}

impl DynamicIndexManager {
    pub fn new(index_path: PathBuf, custom_cache: Option<PathBuf>) -> io::Result<Self> {
        let cache = CacheManager::new(custom_cache)?;
        let metadata_path = index_path.with_extension("bitpop.meta");

        let metadata = if metadata_path.exists() {
            IndexMetadata::from_file(&metadata_path).unwrap_or_else(|_| {
                IndexMetadata::new(8)
            })
        } else {
            IndexMetadata::new(8)
        };

        Ok(Self {
            metadata,
            cache,
            index_path,
            metadata_path,
        })
    }

    pub fn from_existing(index_path: PathBuf) -> io::Result<Self> {
        let meta_path = index_path.with_extension("bitpop.meta");
        if meta_path.exists() {
            let metadata = IndexMetadata::from_file(&meta_path)?;
            let cache = CacheManager::new(None)?;
            Ok(Self {
                metadata,
                cache,
                index_path,
                metadata_path: meta_path,
            })
        } else {
            Self::new(index_path, None)
        }
    }

    pub fn metadata(&self) -> &IndexMetadata {
        &self.metadata
    }

    pub fn metadata_mut(&mut self) -> &mut IndexMetadata {
        &mut self.metadata
    }

    pub fn save_metadata(&self) -> io::Result<()> {
        self.metadata.save(&self.metadata_path)
    }

    pub fn cache(&self) -> &CacheManager {
        &self.cache
    }

    pub fn genome_count(&self) -> u32 {
        self.metadata.genome_count
    }

    pub fn total_bases(&self) -> usize {
        self.metadata.total_bases
    }

    pub async fn fetch_from_ncbi(
        &mut self,
        client: &mut NcbiClient,
        accessions: &[&str],
        _k: usize,
    ) -> io::Result<Vec<PathBuf>> {
        let mut downloaded = Vec::new();

        for accession in accessions {
            println!("  Fetching {} from NCBI...", accession);

            let fasta = match client.fetch_by_accession_version(accession).await {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("  Failed to fetch {}: {}", accession, e);
                    continue;
                }
            };

            let result = self.cache.get_or_download(accession, _k, |_acc| Ok(fasta.clone()))?;

            if let Some(path) = result {
                if let Some(genome) = self.cache.manifest().get(accession) {
                    let checksum = CacheManifest::compute_checksum(&path)?;
                    self.metadata.add_genome(
                        &genome.accession,
                        Some(genome.accession.clone()),
                        Some(path.to_string_lossy().to_string()),
                        genome.genome_size,
                        checksum,
                    );
                    downloaded.push(path.to_path_buf());
                }
            }
        }

        self.save_metadata()?;
        Ok(downloaded)
    }

    pub fn add_local_genome(
        &mut self,
        fasta_path: &Path,
        name: Option<&str>,
        accession: Option<&str>,
    ) -> io::Result<()> {
        let content = fs::read_to_string(fasta_path)?;
        let checksum = CacheManifest::compute_checksum(fasta_path)?;

        let fasta_content = content.as_str();
        let mut lines = fasta_content.lines();
        let header = lines.next().unwrap_or("");
        let genome_name = name.unwrap_or_else(|| {
            parse_header_name(header)
        });

        let sequence_length = fasta_content.len();

        self.metadata.add_genome(
            genome_name,
            accession.map(|s| s.to_string()),
            Some(fasta_path.to_string_lossy().to_string()),
            sequence_length,
            checksum,
        );

        // Cache the sequence
        let acc = accession.unwrap_or("local");
        self.cache.cache_sequence(acc, "1", acc, &content)?;

        self.save_metadata()?;
        Ok(())
    }

    pub fn build_or_update(&self, genomes: Vec<(&str, &str)>, k: usize) -> io::Result<BitPop> {
        let mut bp = BitPop::new(k);

        for (name, seq) in &genomes {
            bp.add_genome(name, seq);
        }

        bp.build();
        Ok(bp)
    }

    pub fn save_index(&self, bp: &BitPop) -> io::Result<()> {
        bp.serialize_to_file(self.index_path.to_str().unwrap())
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))
    }

    pub fn list_genomes(&self) -> Vec<(&u32, &GenomeEntry)> {
        let mut entries: Vec<(&u32, &GenomeEntry)> = self.metadata.genomes.iter().collect();
        entries.sort_by_key(|(k, _)| *k);
        entries
    }

    pub async fn check_for_updates(&mut self, client: &mut NcbiClient) -> Vec<String> {
        let mut updated = Vec::new();

        for (gid, genome) in self.list_genomes() {
            if let Some(ref acc) = genome.accession {
                match client.fetch_by_accession_version(acc).await {
                    Ok(fasta) => {
                        let current_checksum = CacheManifest::compute_checksum(
                            Path::new(genome.local_path.as_deref().unwrap_or(""))
                        ).unwrap_or_default();

                        let new_checksum = {
                            let mut hasher = Sha256::new();
                            hasher.update(fasta.as_bytes());
                            format!("{:x}", hasher.finalize())
                        };

                        if current_checksum != new_checksum {
                            updated.push(acc.clone());
                        }
                    }
                    Err(_) => {}
                }
            }
        }

        updated
    }
}

fn parse_header_name(header: &str) -> &str {
    let trimmed = header.trim_start_matches('>');
    let parts: Vec<&str> = trimmed.split_whitespace().collect();
    if parts.is_empty() {
        return trimmed;
    }
    // Try to extract accession from header
    for part in &parts {
        if part.contains('_') || part.starts_with('N') && part.len() > 5 {
            return part;
        }
    }
    parts[0]
}


