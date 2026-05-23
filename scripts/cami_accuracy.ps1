param(
    [string]$SamFile = "D:\CAMI_low\mapped_base.sam",
    [string]$GroundTruth = "D:\CAMI_low\gs_read_mapping.binning\gs_read_mapping.tsv"
)

Write-Host "=== CAMI Accuracy Test ===" -ForegroundColor Cyan
Write-Host ""

# Parse ground truth - read -> genome mapping
Write-Host "Parsing ground truth..." -ForegroundColor Yellow
$gt = @{}
Get-Content $GroundTruth | Where-Object { $_ -notmatch '^@' -and $_.Trim() } | ForEach-Object {
    $parts = $_.Split("`t")
    if ($parts.Count -ge 2) {
        $readId = $parts[0].Trim()
        $genome = $parts[1].Trim()
        $gt[$readId] = $genome
    }
}
Write-Host "  Loaded $($gt.Count) ground truth mappings" -ForegroundColor Yellow

# Parse SAM - extract read name (strip /1 /2) and genome name
Write-Host "Parsing SAM output..." -ForegroundColor Yellow
$mapped = @{}
$lines = 0
Get-Content $SamFile | Where-Object { $_ -notmatch '^@' } | ForEach-Object {
    $lines++
    $fields = $_.Split()
    if ($fields.Count -lt 3) { return }
    
    $readName = $fields[0]
    $flag = [int]$fields[1]
    
    # Skip unmapped
    if ($flag -band 4) { return }
    
    # Strip /1 /2 suffix
    $readName = $readName -replace '/[12]$', ''
    
    $genome = $fields[2]
    
    # Keep first mapping per read
    if (-not $mapped.ContainsKey($readName)) {
        $mapped[$readName] = $genome
    }
}
Write-Host "  Processed $lines SAM lines" -ForegroundColor Yellow
Write-Host "  Loaded $($mapped.Count) unique read mappings" -ForegroundColor Yellow

# Compare
Write-Host ""
Write-Host "=== Accuracy Results ===" -ForegroundColor Cyan

$correct = 0
$wrong = 0
$unmapped = 0
$gtReads = 0
$wrongPredictions = @{}

foreach ($readId in $gt.Keys) {
    $gtReads++
    $expected = $gt[$readId]
    
    if ($mapped.ContainsKey($readId)) {
        $predicted = $mapped[$readId]
        if ($predicted -eq $expected) {
            $correct++
        } else {
            $wrong++
            if (-not $wrongPredictions.ContainsKey($expected)) {
                $wrongPredictions[$expected] = @{}
            }
            if (-not $wrongPredictions[$expected].ContainsKey($predicted)) {
                $wrongPredictions[$expected][$predicted] = 0
            }
            $wrongPredictions[$expected][$predicted]++
        }
    } else {
        $unmapped++
    }
}

Write-Host "Ground truth reads:    $gtReads"
Write-Host "Mapped reads:          $($gtReads - $unmapped) ($([math]::Round(($gtReads - $unmapped)/$gtReads*100, 1))%)"
Write-Host "Unmapped reads:        $unmapped ($([math]::Round($unmapped/$gtReads*100, 1))%)"
Write-Host ""
Write-Host "Correct predictions:   $correct ($([math]::Round($correct/$gtReads*100, 2))%)"
Write-Host "Wrong predictions:     $wrong ($([math]::Round($wrong/$gtReads*100, 2))%)"
Write-Host ""

# Per-genome breakdown
Write-Host "=== Per-Genome Breakdown ===" -ForegroundColor Cyan
$genomeStats = @{}
foreach ($readId in $gt.Keys) {
    $expected = $gt[$readId]
    if (-not $genomeStats.ContainsKey($expected)) {
        $genomeStats[$expected] = @{total=0; correct=0}
    }
    $genomeStats[$expected].total++
    
    if ($mapped.ContainsKey($readId)) {
        if ($mapped[$readId] -eq $expected) {
            $genomeStats[$expected].correct++
        }
    }
}

Write-Host "$('Genome'.PadRight(40)) $('Total'.PadRight(8)) $('Correct'.PadRight(10)) $('Accuracy')" -ForegroundColor Yellow
Write-Host ("-" * 75)
foreach ($genome in $genomeStats.Keys | Sort-Object) {
    $s = $genomeStats[$genome]
    $acc = [math]::Round($s.correct/$s.total*100, 2)
    Write-Host "$($genome.PadRight(40)) $($s.total.ToString().PadRight(8)) $($s.correct.ToString().PadRight(10)) ${acc}%"
}

# Wrong predictions breakdown
if ($wrongPredictions.Count -gt 0) {
    Write-Host ""
    Write-Host "=== Top Wrong Predictions ===" -ForegroundColor Cyan
    foreach ($expected in $wrongPredictions.Keys | Select-Object -First 10) {
        Write-Host "$expected ->" -ForegroundColor Yellow -NoNewline
        foreach ($predicted in $wrongPredictions[$expected].Keys | Sort-Object { $wrongPredictions[$expected][$_] } -Descending | Select-Object -First 3) {
            Write-Host " $($predicted)($($wrongPredictions[$expected][$predicted]))" -NoNewline
        }
        Write-Host ""
    }
}
