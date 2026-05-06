use reqwest::Client;
use serde::Deserialize;
use std::collections::HashMap;
use std::time::Duration;

const NCBI_BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
const DEFAULT_RATE_LIMIT: f64 = 3.0;

#[derive(Debug, Clone)]
pub struct NcbiConfig {
    api_key: Option<String>,
    rate_limit: f64,
    email: Option<String>,
    timeout_secs: u64,
}

impl Default for NcbiConfig {
    fn default() -> Self {
        Self {
            api_key: None,
            rate_limit: DEFAULT_RATE_LIMIT,
            email: None,
            timeout_secs: 30,
        }
    }
}

impl NcbiConfig {
    pub fn new() -> Self {
        Self {
            api_key: None,
            rate_limit: DEFAULT_RATE_LIMIT,
            email: Some("bit-pop-user@example.com".to_string()),
            timeout_secs: 30,
        }
    }

    pub fn with_api_key(mut self, key: String) -> Self {
        self.api_key = Some(key.clone());
        self.rate_limit = if key.is_empty() { 3.0 } else { 10.0 };
        self
    }

    pub fn with_email(mut self, email: String) -> Self {
        self.email = Some(email);
        self
    }

    pub fn with_rate_limit(mut self, rate: f64) -> Self {
        self.rate_limit = rate;
        self
    }

    pub fn with_timeout(mut self, secs: u64) -> Self {
        self.timeout_secs = secs;
        self
    }

    fn build_query_params(&self) -> HashMap<String, String> {
        let mut params: HashMap<String, String> = HashMap::new();
        if let Some(ref key) = self.api_key {
            params.insert("api_key".to_string(), key.clone());
        }
        if let Some(ref email) = self.email {
            params.insert("email".to_string(), email.clone());
        }
        params
    }
}

#[derive(Debug, Deserialize)]
pub struct ESearchResult {
    #[serde(default)]
    pub count: String,
    #[serde(default)]
    pub retmax: String,
    #[serde(default)]
    pub retstart: String,
    #[serde(default)]
    pub idlist: Vec<String>,
}

impl Default for ESearchResult {
    fn default() -> Self {
        Self {
            count: "0".to_string(),
            retmax: "0".to_string(),
            retstart: "0".to_string(),
            idlist: Vec::new(),
        }
    }
}

#[derive(Debug, Deserialize)]
struct ESearchEnvelope {
    esearchresult: Option<ESearchResult>,
}

impl ESearchResult {
    fn from_envelope(env: ESearchEnvelope) -> Self {
        env.esearchresult.unwrap_or_default()
    }
}

#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct DocSum {
    #[serde(rename = "Id")]
    pub id: String,
    #[serde(rename = "Title")]
    pub title: Option<String>,
    #[serde(rename = "Organism")]
    pub organism: Option<String>,
    #[serde(rename = "NucGeneSim")]
    pub nuc_genesim: Option<String>,
    #[serde(rename = "Pavg")]
    pub pavg: Option<String>,
}

pub struct NcbiClient {
    client: Client,
    config: NcbiConfig,
    last_request_time: Option<std::time::Instant>,
}

impl NcbiClient {
    pub fn new(config: NcbiConfig) -> Self {
        let client = Client::builder()
            .timeout(Duration::from_secs(config.timeout_secs))
            .build()
            .expect("Failed to create HTTP client");

        Self {
            client,
            config,
            last_request_time: None,
        }
    }

    async fn enforce_rate_limit(&mut self) -> Result<(), NcbiError> {
        let now = std::time::Instant::now();
        if let Some(last) = self.last_request_time {
            let min_interval = Duration::from_secs_f64(1.0 / self.config.rate_limit);
            let elapsed = now - last;
            if elapsed < min_interval {
                let sleep_duration = min_interval - elapsed;
                tokio::time::sleep(sleep_duration).await;
            }
        }
        self.last_request_time = Some(std::time::Instant::now());
        Ok(())
    }

    pub async fn search(&mut self, term: &str) -> Result<ESearchResult, NcbiError> {
        self.enforce_rate_limit().await?;

        let mut params = self.config.build_query_params();
        params.insert("db".to_string(), "nucleotide".to_string());
        params.insert("term".to_string(), term.to_string());
        params.insert("retmode".to_string(), "json".to_string());
        params.insert("rettype".to_string(), "summary".to_string());

        let url = format!("{}/esearch.fcgi", NCBI_BASE);
        let response = self.client.get(&url).query(&params).send().await?;

        if !response.status().is_success() {
            let status = response.status();
            let body = response.text().await.unwrap_or_default();
            if status == 429 {
                return Err(NcbiError::RateLimited);
            }
            return Err(NcbiError::HttpError(format!("{}: {}", status, body)));
        }

        let text = response.text().await?;
        let envelope: ESearchEnvelope = serde_json::from_str(&text)
            .map_err(|e| NcbiError::HttpError(format!("JSON: {}: {}", e, text)))?;
        let result = ESearchResult::from_envelope(envelope);
        Ok(result)
    }

    pub async fn fetch_fasta(&mut self, accession: &str) -> Result<String, NcbiError> {
        self.enforce_rate_limit().await?;

        let mut params = self.config.build_query_params();
        params.insert("db".to_string(), "nucleotide".to_string());
        params.insert("id".to_string(), accession.to_string());
        params.insert("rettype".to_string(), "fasta".to_string());
        params.insert("retmode".to_string(), "text".to_string());

        let url = format!("{}/efetch.fcgi", NCBI_BASE);
        let response = self.client.get(&url).query(&params).send().await?;

        if !response.status().is_success() {
            let status = response.status();
            let body = response.text().await.unwrap_or_default();
            if status == 429 {
                return Err(NcbiError::RateLimited);
            }
            return Err(NcbiError::HttpError(format!("{}: {}", status, body)));
        }

        let text = response.text().await?;
        Ok(text)
    }

    pub async fn fetch_fasta_batch(&mut self, accessions: &[&str]) -> Result<Vec<(String, String)>, NcbiError> {
        let mut results = Vec::new();
        for acc in accessions {
            match self.fetch_fasta(acc).await {
                Ok(seq) => results.push((acc.to_string(), seq)),
                Err(e) => {
                    eprintln!("  Warning: failed to fetch {}: {}", acc, e);
                }
            }
        }
        Ok(results)
    }

    pub async fn summary(&mut self, accessions: &[&str]) -> Result<Vec<DocSum>, NcbiError> {
        self.enforce_rate_limit().await?;

        let mut params = self.config.build_query_params();
        params.insert("db".to_string(), "nucleotide".to_string());
        params.insert("id".to_string(), accessions.join(","));
        params.insert("retmode".to_string(), "json".to_string());

        let url = format!("{}/esummary.fcgi", NCBI_BASE);
        let response = self.client.get(&url).query(&params).send().await?;

        if !response.status().is_success() {
            let status = response.status();
            let body = response.text().await.unwrap_or_default();
            if status == 429 {
                return Err(NcbiError::RateLimited);
            }
            return Err(NcbiError::HttpError(format!("{}: {}", status, body)));
        }

        #[derive(Deserialize)]
        struct EsummaryResult {
            docsums: Vec<DocSum>,
        }

        let result: EsummaryResult = response.json().await?;
        Ok(result.docsums)
    }

    pub async fn fetch_by_accession_version(&mut self, accession_version: &str) -> Result<String, NcbiError> {
        let parts: Vec<&str> = accession_version.split('.').collect();
        if parts.len() >= 2 {
            let accession = parts[0];
            let fasta = self.fetch_fasta(accession).await?;
            let expected_header = format!(">{}", accession_version);
            if fasta.starts_with(&expected_header) {
                Ok(fasta)
            } else {
                let fasta2 = self.fetch_fasta(accession_version).await?;
                Ok(fasta2)
            }
        } else {
            self.fetch_fasta(accession_version).await
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum NcbiError {
    #[error("HTTP error: {0}")]
    HttpError(String),
    #[error("Rate limit exceeded by NCBI")]
    RateLimited,
    #[error("Network/error: {0}")]
    Network(#[from] reqwest::Error),
    #[error("No results found for query")]
    NoResults,
}

pub fn build_organism_term(genus: &str, species: &str, strain: Option<&str>) -> String {
    let mut term = format!("\"{} {}\"", genus, species);
    if let Some(s) = strain {
        term.push_str(&format!(" \"{}\"", s));
    }
    term
}

pub fn parse_fasta_header(header: &str) -> &str {
    header
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("")
}

pub fn parse_fasta_sequence(fasta: &str) -> (String, String) {
    let mut lines = fasta.lines();
    let header = lines.next().unwrap_or("");
    let sequence: String = lines.collect();
    (header.trim().to_string(), sequence.trim().to_string())
}
