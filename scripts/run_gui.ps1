param(
  [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"

$Root = Resolve-Path "$PSScriptRoot\.."

if (-not $OutputDir) {
  $OutputDir = "$Root\outputs"
}

python "$Root\gui\laboltz_viewer.py" "$OutputDir"
