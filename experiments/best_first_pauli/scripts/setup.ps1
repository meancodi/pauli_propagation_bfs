param([string]$Julia = "", [string]$Python = "")
$ErrorActionPreference = 'Stop'
$experimentRoot = Split-Path -Parent $PSScriptRoot
$repoRoot = Split-Path -Parent (Split-Path -Parent $experimentRoot)
if (-not $env:JULIA_DEPOT_PATH) { $env:JULIA_DEPOT_PATH = Join-Path $repoRoot '.julia-depot' }
$env:JULIA_NUM_THREADS = '1'
$env:JULIA_NUM_PRECOMPILE_TASKS = '2'
if (-not $Julia) {
    $candidate = Get-Command julia -ErrorAction SilentlyContinue
    if ($candidate -and $candidate.Source -notlike '*WindowsApps*') { $Julia = $candidate.Source }
    else {
        $installed = Get-ChildItem -LiteralPath (Join-Path $env:USERPROFILE '.julia\juliaup') -Filter julia.exe -Recurse -ErrorAction SilentlyContinue | Sort-Object FullName -Descending | Select-Object -First 1
        if ($installed) { $Julia = $installed.FullName }
    }
}
if (-not $Julia) { throw 'Julia 1.10+ was not found. Install Julia via juliaup, or pass -Julia C:\path\to\julia.exe.' }
if (-not $Python) {
    $candidate = Get-Command python -ErrorAction SilentlyContinue
    if ($candidate -and $candidate.Source -notlike '*WindowsApps*') { $Python = $candidate.Source }
    else {
        $installed = Get-ChildItem -Path (Join-Path $env:LOCALAPPDATA 'Programs\Python\Python*\python.exe') -ErrorAction SilentlyContinue | Sort-Object FullName -Descending | Select-Object -First 1
        if ($installed) { $Python = $installed.FullName }
    }
}
if (-not $Python) { throw 'Python is needed for Matplotlib figures. Pass -Python C:\path\to\python.exe.' }
$env:PYTHON = $Python
$env:MPLCONFIGDIR = Join-Path $experimentRoot 'results\work\matplotlib'
& $Julia --startup-file=no --threads=1 (Join-Path $PSScriptRoot 'setup.jl')
if ($LASTEXITCODE -ne 0) { throw 'Julia setup failed; see diagnostics above.' }
& $Python -m pip install -r (Join-Path $experimentRoot 'requirements-plots.txt')
if ($LASTEXITCODE -ne 0) { throw 'Plot dependency installation failed.' }
Write-Host "Setup complete. JULIA_DEPOT_PATH=$env:JULIA_DEPOT_PATH"
Write-Host "Julia: $Julia"
Write-Host "Python: $Python"
Write-Host "Run: & '$Julia' --startup-file=no --threads=1 --project='$experimentRoot' '$PSScriptRoot\run_tests.jl'"
