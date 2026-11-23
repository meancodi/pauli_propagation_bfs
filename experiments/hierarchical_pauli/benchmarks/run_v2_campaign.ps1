param([ValidateSet('smoke','moderate','stress')][string]$Tier='moderate')
$ErrorActionPreference='Stop'
$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
julia --project="$root" --threads=1 "$PSScriptRoot/run_v2_campaign.jl" "--$Tier"
