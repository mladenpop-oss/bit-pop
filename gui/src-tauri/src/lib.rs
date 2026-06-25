use std::path::PathBuf;
use std::io::BufRead;
use tauri::Emitter;
use tokio::process::Command as TokioCommand;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::task;

fn get_bitpop_path() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let repo_root = PathBuf::from(manifest_dir)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();

    if cfg!(windows) {
        repo_root.join("target").join("release").join("bit-pop.exe")
    } else {
        repo_root.join("target").join("release").join("bit-pop")
    }
}

fn get_repo_root() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest_dir)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[tauri::command]
async fn build_index(
    app_handle: tauri::AppHandle,
    genomes: String,
    output: String,
    _k: u32,
    threads: u32,
    cami: bool,
) -> Result<String, String> {
    let bit_pop_path = get_bitpop_path();
    let repo_root = get_repo_root();

    if !bit_pop_path.exists() {
        return Err(format!("bit-pop not found at {:?}", bit_pop_path));
    }

    app_handle.emit("run-started", "Building index...").unwrap();

    let mut args = vec![
        "build".to_string(),
        "-f".to_string(),
        genomes.clone(),
        "-o".to_string(),
        output.clone(),
        "-k".to_string(),
        _k.to_string(),
        "-t".to_string(),
        threads.to_string(),
    ];

    if cami {
        args.push("--cami".to_string());
    }

    let mut child = TokioCommand::new(&bit_pop_path)
        .current_dir(&repo_root)
        .args(&args)
        .env("BITPOP_PROGRESS", "atomic")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to spawn: {}", e))?;

    let stdout = child.stdout.take().ok_or("Failed to take stdout")?;
    let stderr = child.stderr.take().ok_or("Failed to take stderr")?;

    let stdout_handle = app_handle.clone();
    let stdout_task = tokio::spawn(async move { stream_output(stdout_handle, stdout).await });

    let stderr_handle = app_handle.clone();
    let stderr_task = tokio::spawn(async move { stream_stderr(stderr_handle, stderr).await });

    let status = child.wait().await.map_err(|e| format!("Failed to wait: {}", e))?;

    let _ = stdout_task.await;
    let _ = stderr_task.await;

    if status.success() {
        app_handle.emit("run-finished", "Index built successfully!").unwrap();
        Ok("Done!".to_string())
    } else {
        let code = status.code().map(|c| c.to_string()).unwrap_or_default();
        app_handle.emit("run-error", &format!("Build failed with code: {}", code)).unwrap();
        Err(format!("Build failed with code: {}", code))
    }
}

async fn stream_output(
    app_handle: tauri::AppHandle,
    mut stream: tokio::process::ChildStdout,
) {
    let reader = BufReader::new(&mut stream);
    let mut lines = reader.lines();
    let mut last_progress: f64 = 0.0;
    while let Ok(Some(mut line)) = lines.next_line().await {
        // Handle \r (carriage return) in progress bars - take the last segment
        if line.contains('\r') {
            let segments: Vec<&str> = line.split('\r').collect();
            line = segments.last().map(|s| *s).unwrap_or(&line).to_string();
        }
        // Clean ANSI escape codes
        let clean = strip_ansi(&line);
        let _ = app_handle.emit("run-log", &clean);
        if let Some(p) = parse_progress(&clean) {
            if (p - last_progress).abs() > 0.01 {
                last_progress = p;
                let _ = app_handle.emit("run-progress", p);
            }
        }
    }
}

async fn stream_stderr(
    app_handle: tauri::AppHandle,
    mut stream: tokio::process::ChildStderr,
) {
    let reader = BufReader::new(&mut stream);
    let mut lines = reader.lines();
    while let Ok(Some(line)) = lines.next_line().await {
        let clean = strip_ansi(&line);
        let _ = app_handle.emit("run-log", &clean);
        if let Some(p) = parse_progress(&clean) {
            let _ = app_handle.emit("run-progress", p);
        }
    }
}

fn strip_ansi(line: &str) -> String {
    let mut result = String::new();
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\x1b' {
            if chars.peek() == Some(&'[') {
                chars.next();
                loop {
                    match chars.next() {
                        Some(nc) if nc.is_ascii_alphabetic() => break,
                        None => break,
                        _ => {}
                    }
                }
            }
        } else {
            result.push(c);
        }
    }
    result
}

#[tauri::command]
async fn run_map(
    app_handle: tauri::AppHandle,
    index: String,
    reads: String,
    output: String,
    align_mode: String,
    threads: u32,
    use_top_n: bool,
    top_n: u32,
    use_chunk_pct: bool,
    chunk_pct: f64,
) -> Result<String, String> {
    let bit_pop_path = get_bitpop_path();
    let repo_root = get_repo_root();

    if !bit_pop_path.exists() {
        return Err(format!("bit-pop not found at {:?}", bit_pop_path));
    }

    app_handle.emit("run-started", "Mapping reads...").unwrap();

    let top_n_val = if use_top_n { top_n.to_string() } else { "1".to_string() };

    let mut args = vec![
        "map".to_string(),
        "-i".to_string(), index.clone(),
        "-r".to_string(), reads.clone(),
        "-o".to_string(), output.clone(),
        "-a".to_string(), align_mode.clone(),
        "-t".to_string(), threads.to_string(),
        "--top-n".to_string(), top_n_val,
    ];

    if use_chunk_pct {
        args.push("--chunk-pct".to_string());
        args.push(chunk_pct.to_string());
    }

    let mut child = TokioCommand::new(&bit_pop_path)
        .current_dir(&repo_root)
        .args(&args)
        .env("BITPOP_PROGRESS", "atomic")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to spawn: {}", e))?;

    let stdout = child.stdout.take().ok_or("Failed to take stdout")?;
    let stderr = child.stderr.take().ok_or("Failed to take stderr")?;

    let stdout_handle = task::spawn({
        let ah = app_handle.clone();
        async move { stream_output(ah, stdout).await }
    });

    let stderr_handle = task::spawn({
        let ah = app_handle.clone();
        async move { stream_stderr(ah, stderr).await }
    });

    let status = child.wait().await.map_err(|e| format!("Failed to wait: {}", e))?;

    stdout_handle.abort();
    stderr_handle.abort();

    if status.success() {
        app_handle.emit("run-finished", "Mapping completed!").unwrap();
        Ok("Done!".to_string())
    } else {
        let msg = format!("Mapping failed with code: {:?}", status.code());
        app_handle.emit("run-error", &msg).unwrap();
        Err(msg)
    }
}

#[tauri::command]
async fn run_concon(
    app_handle: tauri::AppHandle,
    indexes: Vec<String>,
    reads: String,
    output: String,
    threads: u32,
    use_top_n: bool,
    top_n: u32,
    consensus_top_n: u32,
    use_chunk_pct: bool,
    chunk_pct: f64,
) -> Result<String, String> {
    let bit_pop_path = get_bitpop_path();
    let repo_root = get_repo_root();

    if !bit_pop_path.exists() {
        return Err(format!("bit-pop not found at {:?}", bit_pop_path));
    }

    app_handle.emit("run-started", "Running consensus...").unwrap();

    let mut args = vec!["concon".to_string()];

    for idx in &indexes {
        args.push("-i".to_string());
        args.push(idx.clone());
    }

    args.push("-r".to_string());
    args.push(reads.clone());
    args.push("-o".to_string());
    args.push(output.clone());
    args.push("-t".to_string());
    args.push(threads.to_string());

    if use_top_n {
        args.push("--top-n".to_string());
        args.push(top_n.to_string());
        args.push("--consensus-top-n".to_string());
        args.push(consensus_top_n.to_string());
    }

    args.push("--strategy".to_string());
    args.push("weighted_score".to_string());

    if use_chunk_pct {
        args.push("--chunk-pct".to_string());
        args.push(chunk_pct.to_string());
    }

   let mut child = TokioCommand::new(&bit_pop_path)
        .current_dir(&repo_root)
        .args(&args)
        .env("BITPOP_PROGRESS", "atomic")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to spawn: {}", e))?;

    let stdout = child.stdout.take().ok_or("Failed to take stdout")?;
    let stderr = child.stderr.take().ok_or("Failed to take stderr")?;

    let stdout_handle = task::spawn({
        let ah = app_handle.clone();
        async move { stream_output(ah, stdout).await }
    });

    let stderr_handle = task::spawn({
        let ah = app_handle.clone();
        async move { stream_stderr(ah, stderr).await }
    });

    let status = child.wait().await.map_err(|e| format!("Failed to wait: {}", e))?;

    stdout_handle.abort();
    stderr_handle.abort();

    if status.success() {
        app_handle.emit("run-finished", "Consensus completed!").unwrap();
        Ok("Done!".to_string())
    } else {
        let msg = format!("Consensus failed with code: {:?}", status.code());
        app_handle.emit("run-error", &msg).unwrap();
        Err(msg)
    }
}

fn parse_progress(line: &str) -> Option<f64> {
    // Try "count/total" format: e.g., "999462/999462"
    if let Some(slash_pos) = line.rfind('/') {
        if slash_pos > 0 {
            // Extract total (after last /)
            let total_str = line[slash_pos + 1..].trim().split(|c: char| !c.is_ascii_digit() && c != '/').next()?.trim();
            // Extract count (before last /)
            let before_slash = &line[..slash_pos];
            let count_str = before_slash.split(|c: char| !c.is_ascii_digit() && c != '/').last()?.trim();

            if let (Ok(count), Ok(total)) = (count_str.parse::<f64>(), total_str.parse::<f64>()) {
                if total > 0.0 && count >= 0.0 {
                    return Some(count / total);
                }
            }
        }
    }
    // Fallback: try "% " format
    if let Some(start) = line.find('%') {
        if start > 0 {
            let num_end = start;
            let num_start = line[..num_end].rfind(' ').unwrap_or(0);
            if let Ok(pct) = line[num_start..num_end].trim().parse::<f64>() {
                if pct >= 0.0 && pct <= 100.0 {
                    return Some(pct / 100.0);
                }
            }
        }
    }
    None
}

#[derive(serde::Serialize, Clone)]
pub struct SamStats {
    pub total: usize,
    pub mapped: usize,
    pub unmapped: usize,
    pub genomes: Vec<GenomeStat>,
}

#[derive(serde::Serialize, Clone)]
pub struct GenomeStat {
    pub genome_name: String,
    pub mapped: usize,
    pub pct: f64,
    pub avg_score: f64,
}

#[derive(serde::Serialize, Clone)]
pub struct SamRow {
    pub read_name: String,
    pub genome_name: String,
    pub score: f64,
    pub position: String,
    pub cigar: String,
    pub status: String,
}

#[derive(serde::Serialize, Clone)]
pub struct SamRowsResult {
    pub rows: Vec<SamRow>,
}

fn parse_sam_line(line: &str) -> Option<SamRow> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return None;
    }

    let read_name = fields[0].to_string();
    let flag: u32 = fields[1].parse().ok()?;
    let is_unmapped = (flag & 0x4) != 0 || fields[2] == "*";
    let genome_name = if is_unmapped {
        "unclassified".to_string()
    } else {
        fields[2].to_string()
    };
    let cigar = if is_unmapped {
        String::new()
    } else {
        fields[5].to_string()
    };
    let position = if is_unmapped {
        "*".to_string()
    } else {
        fields[3].to_string()
    };

    // Parse optional tags - concatenated (e.g., AS:f:0.8777RK:f:0.25)
    let mut score: f64 = 0.0;
    if fields.len() > 11 {
        let tags = fields[11..].join("");
        // Find AS:f: or AS:i: and extract value until next tag (2 uppercase + colon)
        for prefix in &["AS:f:", "AS:i:"] {
            if let Some(start) = tags.find(prefix) {
                let rest = &tags[start + prefix.len()..];
                // Find next tag: 2 consecutive uppercase letters followed by colon
                let mut found_end = false;
                for i in 0..rest.len().saturating_sub(3) {
                    let ch = rest.chars().nth(i).unwrap();
                    if ch.is_ascii_uppercase() {
                        let ch2 = rest.chars().nth(i + 1).unwrap_or(' ');
                        let ch3 = rest.chars().nth(i + 2).unwrap_or(' ');
                        if ch2.is_ascii_uppercase() && ch3 == ':' {
                            let val_str = &rest[..i];
                            if let Ok(s) = val_str.parse::<f64>() {
                                score = s;
                            }
                            found_end = true;
                            break;
                        }
                    }
                }
                if !found_end {
                    if let Ok(s) = rest.parse::<f64>() {
                        score = s;
                    }
                }
                break;
            }
        }
    }

    Some(SamRow {
        read_name,
        genome_name,
        score,
        position,
        cigar,
        status: if is_unmapped { "unmapped" } else { "mapped" }.to_string(),
    })
}

#[tauri::command]
async fn parse_sam_stats(app_handle: tauri::AppHandle, path: String) -> Result<SamStats, String> {
    let file = std::fs::File::open(&path).map_err(|e| format!("Cannot open file: {}", e))?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut unique_reads: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut mapped_reads: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut genome_counts: std::collections::HashMap<String, (usize, f64)> = std::collections::HashMap::new();

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }

        if let Some(row) = parse_sam_line(&line) {
            unique_reads.insert(row.read_name.clone());
            if row.status == "mapped" {
                mapped_reads.insert(row.read_name.clone());
                let entry = genome_counts.entry(row.genome_name).or_insert((0, 0.0));
                entry.0 += 1;
                entry.1 += row.score;
            }
        }
    }

    let total = unique_reads.len();
    let mapped = mapped_reads.len();
    let unmapped = total - mapped;

    let mut genomes: Vec<GenomeStat> = genome_counts
        .into_iter()
        .map(|(name, (count, total_score))| {
            let avg_score = if count > 0 { total_score / count as f64 } else { 0.0 };
            let pct = if total > 0 { (count as f64 / total as f64) * 100.0 } else { 0.0 };
            GenomeStat {
                genome_name: name,
                mapped: count,
                pct,
                avg_score,
            }
        })
        .collect();

    genomes.sort_by(|a, b| b.mapped.cmp(&a.mapped));

    app_handle.emit("sam-loaded", "SAM file loaded").unwrap();

    Ok(SamStats {
        total,
        mapped,
        unmapped,
        genomes,
    })
}

#[tauri::command]
async fn parse_sam_rows(
    path: String,
    page: usize,
    page_size: usize,
    filter_text: String,
    filter_min_score: f64,
    filter_mapped: String,
    sort_col: String,
    sort_asc: bool,
) -> Result<SamRowsResult, String> {
    let file = std::fs::File::open(&path).map_err(|e| format!("Cannot open file: {}", e))?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut all_rows: Vec<SamRow> = Vec::new();
    let offset = page * page_size;
    let limit = page_size;

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }

        if let Some(row) = parse_sam_line(&line) {
            if !filter_text.is_empty() && !row.read_name.contains(&filter_text) {
                continue;
            }
            if row.score < filter_min_score {
                continue;
            }
            if filter_mapped != "all" && row.status != filter_mapped {
                continue;
            }

            if all_rows.len() <= offset + limit {
                all_rows.push(row);
            }
        }
    }

    all_rows.sort_by(|a, b| {
        let cmp = match sort_col.as_str() {
            "read_name" => a.read_name.cmp(&b.read_name),
            "genome_name" => a.genome_name.cmp(&b.genome_name),
            "score" => a.score.partial_cmp(&b.score).unwrap_or(std::cmp::Ordering::Equal),
            "status" => a.status.cmp(&b.status),
            _ => std::cmp::Ordering::Equal,
        };
        if sort_asc { cmp } else { cmp.reverse() }
    });

    let start = offset.min(all_rows.len());
    let rows: Vec<SamRow> = all_rows.into_iter().skip(start).take(page_size).collect();

    Ok(SamRowsResult { rows })
}

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            build_index,
            run_map,
            run_concon,
            parse_sam_stats,
            parse_sam_rows,
        ])
        .run(tauri::generate_context!())
        .expect("failed to run app");
}