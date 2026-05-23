# CAMI Benchmark Test Suite
# Runs all tests in parallel with 5 hour total timeout

$BIT_POP = "C:\Users\Daddy\Documents\GitHub\bit-pop\target\release\bit-pop.exe"
$INDEX = "D:\CAMI_low\cami_index.bitpop"
$READS = "D:\CAMI_low\cami_reads_100k_random.fastq"
$GT = "D:\CAMI_low\gs_read_mapping.binning\gs_read_mapping.tsv"
$GENOMES = "D:\CAMI_low\source_genomes_low\source_genomes"
$OUTDIR = "D:\CAMI_low"
$ACC_SCRIPT = "C:\Users\Daddy\Documents\GitHub\bit-pop\scripts\cami_accuracy.py"

$RESULTS = @{}
$START = Get-Date

# Helper to run accuracy check
function Get-Accuracy {
    param([string]$samFile)
    $proc = Start-Process python -ArgumentList $ACC_SCRIPT, $samFile, $GT, $GENOMES -NoNewWindow -PassThru -Wait
    # Read stdout from the process - we'll parse it differently
    return $null
}

# Parse accuracy from script output
function Parse-Accuracy {
    param([string]$output)
    $match = [regex]::Match($output, 'Correct predictions:\s+\d+\s+\((\d+\.?\d*)%\s*\)')
    if ($match.Success) {
        return [double]$match.Groups[1].Value
    }
    return $null
}

function Parse-Mapped {
    param([string]$output)
    $match = [regex]::Match($output, 'Mapping complete:\s+(\d+)\/(\d+) reads mapped')
    if ($match.Success) {
        return [int]$match.Groups[1].Value
    }
    return $null
}

# ============================
# TEST DEFINITIONS
# ============================
$tests = @(
    # === BASELINE ===
    @{ Name="base_xor_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_base.sam`" -a xor -t 16"; Desc="xor top-n 1" },
    @{ Name="base_xor_tn2"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_base_tn2.sam`" -a xor -t 16 --top-n 2"; Desc="xor top-n 2" },
    @{ Name="base_xor_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_base_tn3.sam`" -a xor -t 16 --top-n 3"; Desc="xor top-n 3" },
    @{ Name="base_hybrid_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_hybrid.sam`" -a hybrid -t 16"; Desc="hybrid top-n 1" },
    @{ Name="base_hybrid_tn2"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_hybrid_tn2.sam`" -a hybrid -t 16 --top-n 2"; Desc="hybrid top-n 2" },
    @{ Name="base_hybrid_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_hybrid_tn3.sam`" -a hybrid -t 16 --top-n 3"; Desc="hybrid top-n 3" },

    # === ALIGN MODES ===
    @{ Name="xor_sw_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_xor_sw.sam`" -a sw -t 16"; Desc="sw top-n 1" },
    @{ Name="xor_sw_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_xor_sw_tn3.sam`" -a sw -t 16 --top-n 3"; Desc="sw top-n 3" },

    # === SPACED SEEDS ===
    @{ Name="ss_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_ss.sam`" -a xor -t 16 --spaced-seed --top-n 1"; Desc="spaced seed top-n 1" },
    @{ Name="ss_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_ss_tn3.sam`" -a xor -t 16 --spaced-seed --top-n 3"; Desc="spaced seed top-n 3" },
    @{ Name="ss_custom_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_ss_custom.sam`" -a xor -t 16 --spaced-seed --spaced-seed-pattern `"(1111111111){0,1}(1111111111){0,1}(1111111111){0,1}`" --top-n 3"; Desc="custom spaced seed" },

    # === SNAPSHOT ===
    @{ Name="snp_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_snp.sam`" -a xor -t 16 --snp-detect --snp-min-support 3"; Desc="snp top-n 1" },
    @{ Name="snp_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_snp_tn3.sam`" -a xor -t 16 --snp-detect --snp-min-support 3 --top-n 3"; Desc="snp top-n 3" },

    # === GOLDEN ANCHORS ===
    @{ Name="golden_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_golden.sam`" -a xor -t 16 --golden-anchors --top-n 1"; Desc="golden anchors" },
    @{ Name="golden_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_golden_tn3.sam`" -a xor -t 16 --golden-anchors --top-n 3"; Desc="golden anchors tn3" },

    # === SEARCH RADIUS ===
    @{ Name="sr5_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_sr5.sam`" -a xor -t 16 --top-n 3 --search-radius 5"; Desc="search-radius 5" },
    @{ Name="sr10_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_sr10.sam`" -a xor -t 16 --top-n 3 --search-radius 10"; Desc="search-radius 10" },
    @{ Name="sr20_tn3"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_sr20.sam`" -a xor -t 16 --top-n 3 --search-radius 20"; Desc="search-radius 20" },

    # === CHAIN MODE ===
    @{ Name="chain_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_chain.sam`" -a chain -t 16"; Desc="chain mode" },

    # === EM POST-PROCESSING ===
    @{ Name="em_base_t01"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base.sam`" -o `"$OUTDIR\test_em_t01.sam`" --temperature 0.1"; Desc="EM t=0.1" },
    @{ Name="em_base_t05"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base.sam`" -o `"$OUTDIR\test_em_t05.sam`" --temperature 0.5"; Desc="EM t=0.5" },
    @{ Name="em_base_t10"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base.sam`" -o `"$OUTDIR\test_em_t10.sam`" --temperature 1.0"; Desc="EM t=1.0" },
    @{ Name="em_tn3_t01"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base_tn3.sam`" -o `"$OUTDIR\test_em_tn3_t01.sam`" --temperature 0.1"; Desc="EM tn3 t=0.1" },
    @{ Name="em_tn3_t01_ct09"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base_tn3.sam`" -o `"$OUTDIR\test_em_tn3_t01_ct09.sam`" --temperature 0.1 --confidence-threshold 0.9"; Desc="EM tn3 t=0.1 ct=0.9" },
    @{ Name="em_tn3_t01_ct095"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base_tn3.sam`" -o `"$OUTDIR\test_em_tn3_t01_ct095.sam`" --temperature 0.1 --confidence-threshold 0.95"; Desc="EM tn3 t=0.1 ct=0.95" },
    @{ Name="em_tn3_t05_ct095"; Cmd="$BIT_POP em -i `"$OUTDIR\test_base_tn3.sam`" -o `"$OUTDIR\test_em_tn3_t05_ct095.sam`" --temperature 0.5 --confidence-threshold 0.95"; Desc="EM tn3 t=0.5 ct=0.95" },

    # === AUTO K ===
    @{ Name="autok_tn1"; Cmd="$BIT_POP map -i `"$INDEX`" -r `"$READS`" -o `"$OUTDIR\test_autok.sam`" -a xor -t 16 --auto-k"; Desc="auto-k" },

    # === CHUNK-CONSENSUS ===
    @{ Name="cc_1_10_50"; Cmd="$BIT_POP chunk-consensus -i `"$INDEX`" -c 0.01,0.10,0.50 -r `"$READS`" -o `"$OUTDIR\test_cc.sam`" -t 16"; Desc="chunk-consensus 1/10/50" },
    @{ Name="cc_1_5_10"; Cmd="$BIT_POP chunk-consensus -i `"$INDEX`" -c 0.01,0.05,0.10 -r `"$READS`" -o `"$OUTDIR\test_cc2.sam`" -t 16"; Desc="chunk-consensus 1/5/10" },

    # === HF (needs HF index) ===
    @{ Name="hf_tn3"; Cmd="$BIT_POP map -i `"$OUTDIR\cami_index_hf.bitpop`" -r `"$READS`" -o `"$OUTDIR\test_hf.sam`" -a xor -t 16 --top-n 3 --hf --hf-min 3"; Desc="hf tn3" },
    @{ Name="hf_snp_tn3"; Cmd="$BIT_POP map -i `"$OUTDIR\cami_index_hf.bitpop`" -r `"$READS`" -o `"$OUTDIR\test_hf_snp.sam`" -a xor -t 16 --top-n 3 --hf --hf-min 3 --snp-detect --snp-min-support 3"; Desc="hf+snp tn3" },
)

Write-Host "=== CAMI Benchmark Test Suite ===" -ForegroundColor Cyan
Write-Host "Tests: $($tests.Count)" -ForegroundColor Cyan
Write-Host "Started: $($START.ToString('HH:mm:ss'))" -ForegroundColor Cyan
Write-Host ""

# Run tests with limited concurrency (8 parallel)
$MAX_PARALLEL = 8
$running = @()
$completed = @{}
$testIndex = 0

while ($testIndex -lt $tests.Count -or $running.Count -gt 0) {
    # Start new tests if under limit
    while ($running.Count -lt $MAX_PARALLEL -and $testIndex -lt $tests.Count) {
        $test = $tests[$testIndex]
        Write-Host "[$($running.Count+1)/$MAX_PARALLEL] Starting: $($test.Name) - $($test.Desc)" -ForegroundColor Yellow

        $proc = Start-Process powershell -ArgumentList "-NoProfile", "-Command", $test.Cmd -NoNewWindow -PassThru
        $running += @{
            Proc = $proc
            Test = $test
            StartTime = Get-Date
        }
        $testIndex++
    }

    # Check for completed tests
    $stillRunning = @()
    foreach ($r in $running) {
        if ($rProc.HasExited) {
            $test = $r.Test
            $samFile = [regex]::Match($test.Cmd, '-o\s+"([^"]+)"').Groups[1].Value
            $elapsed = (Get-Date) - $r.StartTime

            Write-Host "[DONE] $($test.Name) - $($test.Desc) ($([math]::Round($elapsed.TotalSeconds))s)" -ForegroundColor Green

            # Run accuracy check
            if ($samFile -and (Test-Path $samFile)) {
                $accOut = & python $ACC_SCRIPT $samFile $GT $GENOMES 2>&1
                $mapped = Parse-Mapped $accOut
                $accuracy = Parse-Accuracy $accOut
                $RESULTS[$test.Name] = @{
                    Desc = $test.Desc
                    Mapped = $mapped
                    Accuracy = $accuracy
                    Time = [math]::Round($elapsed.TotalSeconds)
                }
            }
            $completed[$test.Name] = $true
        } else {
            $stillRunning += $r
        }
    }
    $running = $stillRunning

    if ($running.Count -gt 0) {
        Start-Sleep -Milliseconds 5000
    }

    # Safety check - 5 hour timeout
    $totalElapsed = (Get-Date) - $START
    if ($totalElapsed.TotalHours -ge 5) {
        Write-Host "" -ForegroundColor Red
        Write-Host "=== 5 HOUR TIMEOUT REACHED ===" -ForegroundColor Red
        Write-Host "Killing remaining processes..." -ForegroundColor Red
        foreach ($r in $running) { $r.Proc | Stop-Process -Force }
        break
    }
}

# ============================
# RESULTS SUMMARY
# ============================
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "RESULTS SUMMARY" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$sorted = $RESULTS.GetEnumerator() | Sort-Object { $_.Value.Accuracy } -Descending

foreach ($r in $sorted) {
    $accStr = if ($r.Value.Accuracy) { "{0:N2}%" -f $r.Value.Accuracy } else { "N/A" }
    $mappedStr = if ($r.Value.Mapped) { $r.Value.Mapped.ToString() } else { "N/A" }
    $timeStr = if ($r.Value.Time) { "${0}s" -f $r.Value.Time } else { "N/A" }
    Write-Host "$($r.Key,-25) $($r.Value.Desc,-35) Mapped: $($mappedStr,-12) Accuracy: $($accStr,-10) Time: $($timeStr)"
}

$END = Get-Date
$TOTAL = $END - $START
Write-Host ""
Write-Host "Total time: $([math]::Round($TOTAL.TotalMinutes)) minutes" -ForegroundColor Cyan
Write-Host "Completed: $($completed.Count)/$($tests.Count)" -ForegroundColor Cyan

# Save results
$RESULTS | ConvertTo-Json | Out-File "$OUTDIR\test_results.json" -Encoding utf8
Write-Host "Results saved to $OUTDIR\test_results.json" -ForegroundColor Green
