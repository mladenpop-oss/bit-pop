param(
    [Parameter(Mandatory=$true)]
    [string]$SamFile
)

Write-Host "Analyzing: $SamFile" -ForegroundColor Cyan
$size = [math]::Round((Get-Item $SamFile).Length/1MB, 1)
Write-Host "File size: ${size} MB"

$headers = 0
$mapped = 0
$unmapped = 0
$rc = 0
$supplementary = 0
$readNames = @{}

$i = 0
Get-Content $SamFile | ForEach-Object {
    $i++
    if ($_ -match '^@') {
        $headers++
        return
    }
    
    $fields = $_.Split()
    $name = $fields[0]
    $flag = [int]$fields[1]
    
    if (-not $readNames.ContainsKey($name)) {
        $readNames[$name] = $true
    }
    
    if ($flag -band 4) {
        $unmapped++
    } else {
        $mapped++
        if ($flag -band 16) { $rc++ }
        if ($flag -band 2048) { $supplementary++ }
    }
    
    if ($i % 200000 -eq 0) {
        Write-Host "  Processed: $i lines..." -NoNewline
        Write-Host " (mapped: $mapped, unmapped: $unmapped)" -NoNewline
    }
}

$totalReads = $readNames.Count
Write-Host ""
Write-Host "=== Results ===" -ForegroundColor Green
Write-Output "Header lines:     $headers"
Write-Output "Total SAM lines:  $i"
Write-Output "Unique reads:     $totalReads"
Write-Output "Mapped reads:     $mapped ($([math]::Round($mapped/$totalReads*100, 1))%)"
Write-Output "Unmapped reads:   $unmapped ($([math]::Round($unmapped/$totalReads*100, 1))%)"
Write-Output "Reverse comp:     $rc"
Write-Output "Supplementary:    $supplementary"
