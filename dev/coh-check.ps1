# R CMD check --as-cran on frmtmb.coupling, ONCE, on the final pass.
#
# pandoc and TinyTeX have to be on PATH or the check reports a bogus
# pdflatex ERROR; the CRAN-incoming remote check is off because frmtmb
# is not on CRAN and its WARNING is pre-existing.
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$pandoc = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:RSTUDIO_PANDOC = $pandoc
$env:PATH = "$pandoc;" +
  "C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/coh-lib;" +
  "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;" +
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = ""
$out = Join-Path $root "dev/coh-check"
if (-not (Test-Path $out)) { New-Item -ItemType Directory $out | Out-Null }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build extensions/frmtmb.coupling
$tar = Get-ChildItem -Path $root -Filter "frmtmb.coupling_*.tar.gz" |
  Sort-Object LastWriteTime | Select-Object -Last 1
Write-Output "built $($tar.Name)"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$out $tar.FullName
