# Bit-Pop Usage Guide

## Quick Start

```bash
# One-command workflow: build index + map reads
bit-pop run data/genomes/ -r reads.fastq

# Paired-end
bit-pop run data/genomes/ -1 R1.fastq -2 R2.fastq

# Download from NCBI and map
bit-pop run --ncbi "Escherichia coli" reads.fastq
```

## Commands

### `run` — One-Command Workflow

Builds an index (if needed) and maps reads in a single command.

```bash
bit-pop run <GENOME> [OPTIONS]
```

**Arguments:**
- `<GENOME>` — FASTA file, folder of FASTA files, or organism name with `--ncbi`

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --index` | Use existing `.bitpop` index | build from genomes |
| `-r, --reads` | Reads file (FASTQ/FASTA) | required |
| `-1, --reads-1` | R1 FASTQ (paired-end) | — |
| `-2, --reads-2` | R2 FASTQ (paired-end) | — |
| `-n, --ncbi` | Fetch genome from NCBI | — |
| `-o, --output` | Output SAM file | `<reads_name>.sam` |
| `-k, --k` | K-mer size | 10 |
| `--auto-k` | Auto-scale k by genome size | — |
| `-a, --align-mode` | xor, sw, hybrid | hybrid |
| `-m, --min-score` | Minimum alignment score (0-1) | 0.7 |
| `-t, --threads` | Number of threads | 1 |
| `--top-n` | Top rarest k-mer anchors | 1 |
| `--em` | Apply EM post-processing | — |

**Examples:**
```bash
# Simple mapping
bit-pop run data/genomes/Ecoli_K12.fna -r reads.fastq

# With EM post-processing
bit-pop run data/genomes/ -r reads.fastq --em

# NCBI download + map
bit-pop run --ncbi "Escherichia coli" -r reads.fastq --email "user@example.com"
```

---

### `build` — Build FM-Index

```bash
bit-pop build -f <FASTA> -o <OUTPUT> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-f, --fasta` | Input FASTA file(s) | required |
| `-o, --output` | Output index path | required |
| `-k, --k` | K-mer size | 8 |
| `--auto-k` | Auto-scale k by genome size | — |
| `--read-type` | short (Illumina) or long (Nanopore/PacBio) | short |
| `-t, --threads` | Number of threads | 1 |
| `--cami` | CAMI mode: extract genome name from filename | — |
| `--pacbio` | PacBio mode: extract genome name from filename | — |
| `--hf` | Homopolymer fingerprint scoring | — |
| `--spaced-seed` | Spaced seed pattern | — |

**Examples:**
```bash
# Single genome
bit-pop build -f genome.fna -o index.bitpop -k 21

# Multiple genomes
bit-pop build -f genome1.fna genome2.fna genome3.fna -o index.bitpop -k 21

# CAMI mode
bit-pop build -f *.fasta -o index.bitpop --cami --auto-k
```

---

### `map` — Map Reads

```bash
bit-pop map -i <INDEX> -r <READS> -o <OUTPUT> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --index` | Input index path | required |
| `-r, --reads` | Input reads (FASTA/FASTQ) | required |
| `-o, --output` | Output SAM file | required |
| `-m, --min-score` | Minimum alignment score (0-1) | 0.7 |
| `-a, --align-mode` | xor, sw, hybrid | xor |
| `-q, --min-quality` | Minimum average quality (0 = no filter) | 0 |
| `-t, --reads-threads` | Number of threads | 1 |
| `--top-n` | Top rarest k-mer anchors | 1 |
| `--em` | Apply EM post-processing | — |

**Chunking Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `--chunk-size` | Fixed chunk size (0 = auto) | 0 |
| `--chunk-pct` | Chunk as % of read length (0 = disabled) | 0 |
| `--chunk-min` | Minimum chunk size clamp | 20 |
| `--chunk-max` | Maximum chunk size clamp | 500 |
| `--chunk-vote-threshold` | Min fraction of chunks to agree (0 = none) | 0 |
| `--chunk-top-n` | Top genomes per read | 1 |
| `--chunk-strategy` | rarest, golden, spaced | rarest |
| `--score-mode` | quality (score\*score), base (raw sum) | quality |
| `--anchor-min-score` | Min anchor score threshold | 0.5 |
| `--anchor-filter` | Use legacy anchor_filter | — |

**Advanced Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `--bam` | Output BAM format | — |
| `--stream` | Stream reads (limits memory) | — |
| `--max-ram` | Max RAM for streaming (e.g., "32G") | — |
| `--two-pass` | Re-map unmapped with lower threshold | — |
| `--second-pass-score` | Min score for second pass | 0.5 |
| `--diagnose-unmapped` | Diagnose why reads failed | — |
| `--snp-detect` | SNP-aware scoring | — |
| `--hf` | Homopolymer fingerprint | — |

**Examples:**
```bash
# Basic mapping
bit-pop map -i index.bitpop -r reads.fastq -o output.sam -t 16

# Clinical sample (optimal settings)
bit-pop map -i index_k21.bitpop -r clinical.fastq -o output.sam \
  --top-n 4 -t 16 --chunk-pct 0.03 --chunk-min 125 --chunk-max 130

# Long reads with chunking
bit-pop map -i index.bitpop -r nanopore.fastq -o output.sam \
  --chunk-pct 0.02 --chunk-min 20 --chunk-max 500

# Two-pass mapping
bit-pop map -i index.bitpop -r reads.fastq -o output.sam --two-pass
```

---

### `em` — EM Post-Processing

Improves strain resolution through soft-assignment classification.

```bash
bit-pop em -i <INPUT> -o <OUTPUT> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input SAM file | required |
| `-o, --output` | Output SAM file | required |
| `--convergence` | KL divergence threshold | 0.001 |
| `--max-iterations` | Maximum EM iterations | 50 |
| `--temperature` | Softmax temperature (lower = sharper) | 0.1 |
| `--top-k` | Top-K genomes per read | 10 |
| `--confidence-threshold` | Min probability to apply (0 = always) | 0.0 |

**Examples:**
```bash
# Basic EM
bit-pop em -i mapped.sam -o mapped_em.sam

# Recommended for clinical samples
bit-pop em -i mapped.sam -o mapped_em.sam --top-k 10 --temperature 0.1

# Conservative (only high-confidence)
bit-pop em -i mapped.sam -o mapped_em.sam --confidence-threshold 0.75
```

---

### `consensus` — Multi-K Consensus

Map reads against multiple k-indexes with voting.

```bash
bit-pop consensus -i <INDEXES>... -r <READS> -o <OUTPUT> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --indexes` | List of index files | required |
| `-r, --reads` | Reads file (FASTQ) | required |
| `-o, --output` | Output SAM file | required |
| `--strategy` | weighted_score, majority, best_score, base_score | weighted_score |
| `--min-score` | Min alignment score (0 = no filter) | 0.0 |
| `--min-k-mappings` | Min k-values that must map (0 = any) | 1 |
| `--top-n` | Top candidates per read (0 = winner only) | 1 |
| `--two-pass` | Map each k separately (faster) | — |

**Examples:**
```bash
# Multi-k consensus
bit-pop consensus -i index_k13.bitpop index_k15.bitpop \
  -r reads.fastq -o output.sam --strategy weighted_score -t 16

# With chunking
bit-pop consensus -i index_k13.bitpop index_k15.bitpop \
  -r reads.fastq -o output.sam --chunk-pct 0.02 -t 16
```

---

### `concon` — Consensus

Runs `bit-pop map` for each index, then combines results (subprocess-based, like Python script).

```bash
bit-pop concon -i <INDEXES>... -r <READS> -o <OUTPUT> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --indexes` | List of index files | required |
| `-r, --reads` | Reads file (FASTQ) | required |
| `-o, --output` | Output SAM file | required |
| `--strategy` | weighted_score, majority, best_score, base_score | weighted_score |
| `--min-score` | Min alignment score (0 = no filter) | 0.0 |
| `-t, --threads` | Threads per map | 1 |
| `--top-n` | Top-N rarest k-mer anchors per map | 1 |
| `--consensus-top-n` | Top-N consensus candidates per read | 1 |
| `--chunk-pct` | Chunk as % of read length | 0 |
| `--chunk-min` | Min chunk size | 20 |
| `--chunk-max` | Max chunk size | 500 |

**Examples:**
```bash
# Fast consensus with two indexes
bit-pop concon -i index_k17.bitpop index_k15.bitpop \
  -r reads.fastq -o output.sam --strategy weighted_score \
  --top-n 4 -t 16 --consensus-top-n 2

# With chunking
bit-pop concon -i index_k17.bitpop index_k15.bitpop \
  -r reads.fastq -o output.sam --chunk-pct 0.03 \
  --chunk-min 125 --chunk-max 200 -t 16
```

---

### `load` — Incremental Index Update

Add genomes to an existing index.

```bash
bit-pop load -i <INDEX> -f <FASTA> -o <OUTPUT>
```

**Examples:**
```bash
# Add new genome to existing index
bit-pop load -i index.bitpop -f new_genome.fna -o updated.bitpop
```

---

### `stats` — Index Statistics

Show information about an index.

```bash
bit-pop stats -i <INDEX>
```

**Examples:**
```bash
# Show index statistics
bit-pop stats -i index.bitpop
```

---

### `search` — NCBI Search

Search for genome accessions by organism name.

```bash
bit-pop search <ORGANISM>
```

**Examples:**
```bash
# Search for E. coli genomes
bit-pop search "Escherichia coli"
```

---

### `fetch` — NCBI Fetch + Build

Fetch genome sequences from NCBI and build index.

```bash
bit-pop fetch <ORGANISM> -o <OUTPUT> [OPTIONS]
```

**Examples:**
```bash
# Fetch and build index
bit-pop fetch "Escherichia coli" -o index.bitpop --email "user@example.com"
```

---

### `update` — Update Cached Genomes

Update cached genomes with latest versions from NCBI.

```bash
bit-pop update [OPTIONS]
```

**Examples:**
```bash
# Update all cached genomes
bit-pop update --email "user@example.com"
```

---

### `tax` — Taxonomic Classification

Generate taxonomic classification report from SAM output.

```bash
bit-pop tax -i <INPUT> --nodes-dmp <NODES> --names-dmp <NAMES> [OPTIONS]
```

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input SAM file | required |
| `--nodes-dmp` | NCBI nodes.dmp file | required |
| `--names-dmp` | NCBI names.dmp file | required |
| `-o, --output` | Output file (default: stdout) | — |
| `--top-n` | Top entries per rank | 10 |
| `--format` | text, json | text |

**Examples:**
```bash
# Taxonomic report
bit-pop tax -i mapped.sam --nodes-dmp nodes.dmp --names-dmp names.dmp

# JSON output
bit-pop tax -i mapped.sam --nodes-dmp nodes.dmp --names-dmp names.dmp \
  --format json -o report.json
```

---

## Recommended Configurations

### Clinical Sample Classification (Best Accuracy)

```bash
bit-pop map -i index_k21.bitpop -r clinical.fastq -o output.sam \
  --top-n 4 -t 16 --chunk-pct 0.03 --chunk-min 125 --chunk-max 130

bit-pop em -i output.sam -o final.sam --top-k 10 --temperature 0.1
```

**Results:** 99.3% overall accuracy, 99.5% Ebola, 67 FP (11,685 reads)

### Ebola Strain Identification

```bash
bit-pop concon -i ebola_k13.bitpop ebola_k15.bitpop \
  -r reads.fastq -o output.sam \
  --strategy weighted_score --top-n 4 \
  --chunk-pct 0.02 --consensus-top-n 2 -t 16
```

**Results:** 99.98% strain accuracy, 81% mapping rate

### CAMI Metagenomic Classification

**Precision mode:**
```bash
bit-pop map -i index_k13.bitpop -r reads.fastq -o output.sam -t 16
bit-pop em -i output.sam -o final.sam --temperature 1.0
```

**Balanced mode:**
```bash
bit-pop consensus -i index_k12.bitpop index_k13.bitpop index_k14.bitpop index_k15.bitpop \
  -r reads.fastq -o output.sam --top-n 2 -t 16
bit-pop em -i output.sam -o final.sam --temperature 1.0
```

**Coverage mode:**
```bash
bit-pop consensus -i index_k13.bitpop index_k22.bitpop \
  -r reads.fastq -o output.sam --top-n 2 -t 16
bit-pop em -i output.sam -o final.sam --temperature 1.0
```

### PacBio HiFi

```bash
bit-pop map -i index_k70.bitpop -r hifi.fastq -o output.sam -t 16
```

No chunking needed for reads >1,000 bp.

---

## Output Format

Bit-Pop outputs SAM format with custom tags:

| Tag | Description |
|-----|-------------|
| `AS:i` | Alignment score |
| `NM:i` | Number of mismatches |
| `XS:i` | Strand (1 = forward, -1 = reverse) |
| `YT:A` | Mapping type (U = unmapped, M = mapped) |

Genome name is encoded in the reference name field:
```
NC_014373.1 Bundibugyo ebolavirus, complete genome:12345
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ^^^^^
genome name                                             position
```

---

## Performance Tips

1. **Use k21 for clinical samples** — optimal accuracy/FP trade-off
2. **Narrow chunk range (125-130bp)** — best for mixed samples
3. **Always use EM post-processing** — +1.5-3.8% accuracy, -52% FP
4. **Single index for outbreak detection** — stronger unmapped signal than consensus
5. **Multi-k consensus for strain resolution** — k13+k15 or k13+k22
6. **Use `--stream` for large FASTQ** — limits memory usage
7. **`--two-pass` for difficult samples** — re-maps unmapped reads with lower threshold
