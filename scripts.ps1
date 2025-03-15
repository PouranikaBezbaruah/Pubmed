# run-pubmed-fetcher.ps1

param(
    [Parameter(Mandatory=$true)]
    [string]$Query,
    
    [Parameter(Mandatory=$false)]
    [string]$OutputFile,
    
    [Parameter(Mandatory=$false)]
    [switch]$Debug,
    
    [Parameter(Mandatory=$false)]
    [int]$MaxResults = 100
)

# Construct the command
$command = "get-papers-list `"$Query`""

if ($OutputFile) {
    $command += " -f `"$OutputFile`""
}

if ($Debug) {
    $command += " -d"
}

if ($MaxResults -ne 100) {
    $command += " -m $MaxResults"
}

# Run the command using Poetry
Write-Host "Running command: $command"
poetry run $command