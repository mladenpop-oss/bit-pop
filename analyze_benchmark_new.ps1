$files = @(
    @{Name="E. coli"; File="bench_ecoli_new.sam"; Prefix="eco"},
    @{Name="S. aureus"; File="bench_aureus_new.sam"; Prefix="aureus"},
    @{Name="S. cerevisiae"; File="bench_cerevisiae_new.sam"; Prefix="cerevisiae"}
)

$genome_map = @{
    "eco" = "NC_000913.3"
    "aureus" = "CP029198"
    "cerevisiae" = "Sac_cerevisiae"
}

Write-Host "=== BENCHMARK REZULTATI (k=10) ===" -ForegroundColor Cyan
Write-Host ""

$total_correct = 0
$total_mapped = 0
$total_reads = 0

foreach ($f in $files) {
    $lines = Get-Content $f.File
    $mapped = 0; $correct = 0; $wrong = 0
    $target_genome = $genome_map[$f.Prefix]
    
    foreach ($line in $lines) {
        if ($line -match '^@') { continue }
        if ($line -match '^[ACGT]+$') { continue }
        
        $fields = $line -split "`t"
        $name = $fields[0]
        $ref = $fields[2]
        
        if ($ref -eq '*') { continue }
        
        $mapped++
        $total_reads++
        
        if ($ref -match $target_genome) {
            $correct++
            $total_correct++
        } else {
            $wrong++
            $total_mapped++
            if ($wrong -le 3) {
                Write-Host "  WRONG: $($f.Name) read -> $ref"
            }
        }
    }
    
    $accuracy = if ($mapped -gt 0) { [math]::Round($correct/$mapped*100, 1) } else { 0 }
    Write-Host "$($f.Name): $mapped mapped / accuracy: $accuracy% ($correct/$mapped correct, $wrong wrong)"
}

Write-Host ""
Write-Host "=== UKUPNO ===" -ForegroundColor Yellow
Write-Host "Total mapped: $total_mapped"
Write-Host "Total correct: $total_correct"
Write-Host "Overall accuracy: $([math]::Round($total_correct/$total_mapped*100, 1))%"
