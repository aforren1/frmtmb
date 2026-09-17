# R CMD check --as-cran on frmtmb, ONCE, on the final pass.
#
# pandoc and TinyTeX have to be on PATH or the check reports a bogus
# pdflatex ERROR; the CRAN-incoming remote check is off because frmtmb
# is not on CRAN and its WARNING is pre-existing. No --no-manual: the
# manual sections are where an unescaped % in an Rd surfaces.
# NOT "Stop": R CMD build writes progress to stderr, and PowerShell 5.1
# wraps a native command's stderr in an ErrorRecord, so "Stop" aborted
# this script after a SUCCESSFUL build and never ran the check.
$ErrorActionPreference = "Continue"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$pandoc = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:RSTUDIO_PANDOC = $pandoc
$env:PATH = "$pandoc;" +
  "C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/priorform-lib;" +
  "C:/Users/adf44/source/r/pinlib;" +
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = ""
$env:FRMTMB_STAN_CACHE = Join-Path $root "dev/stan-cache"
# --output= needs the directory to EXIST. Deleting it before a rerun
# and forgetting this gives "Error in setwd(outdir)" after a build that
# took four minutes.
$out = Join-Path $root "dev/priorform-check"
if (-not (Test-Path $out)) { New-Item -ItemType Directory $out | Out-Null }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build .
$tar = Get-ChildItem -Path $root -Filter "frmtmb_*.tar.gz" |
  Sort-Object LastWriteTime | Select-Object -Last 1
Write-Output "built $($tar.Name)"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$out $tar.FullName
