# Fix FASTA headers and group genomes into 3 equal groups
# Usage: .\fix_headers_and_group.ps1

$sourceDir = "G:\FA"
$targetBase = "G:\fn_fixed"

# Get all genome files
$genomes = Get-ChildItem $sourceDir -Filter "*.fa" -File
Write-Host "Found $($genomes.Count) genome files"

# Shuffle genomes for even distribution
$shuffled = $genomes | Sort-Object { [guid]::NewGuid() }

# Create 3 groups (roughly equal)
$numGroups = 3
$groupSize = [math]::Ceiling($shuffled.Count / $numGroups)
$groups = @()
for ($i = 0; $i -lt $numGroups; $i++) {
    $start = $i * $groupSize
    $end = [math]::Min($start + $groupSize, $shuffled.Count)
    $groups += , $shuffled[$start..($end - 1)]
}

# Create target directories and process
for ($g = 1; $g -le $numGroups; $g++) {
    $groupDir = Join-Path $targetBase "group$g"
    New-Item -ItemType Directory -Path $groupDir -Force | Out-Null
    
    $genomeList = $groups[$g - 1]
    Write-Host "Group ${g}: $($genomeList.Count) genomes"
    
    foreach ($genome in $genomeList) {
        # Extract genome name from filename (remove .fa extension)
        $genomeName = $genome.BaseName
        
        # Read the file content
        $content = Get-Content $genome.FullName -Raw
        
        # Replace the first line (header) with the genome name
        $lines = $content -split "`r?`n"
        $lines[0] = ">$genomeName"
        $newContent = $lines -join "`r`n"
        
        # Write to target directory
        $targetPath = Join-Path $groupDir "$genomeName.fa"
        Set-Content -Path $targetPath -Value $newContent -NoNewline -Encoding UTF8
    }
}

Write-Host ""
Write-Host "Done! Files saved to: $targetBase"
Write-Host "Total groups: $numGroups"
$totalFiles = (Get-ChildItem $targetBase -Recurse -File).Count
Write-Host "Total files: $totalFiles"
