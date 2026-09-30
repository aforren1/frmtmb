# Run one R script of this lane with the toolchain on PATH and TMP set.
# Usage: powershell -File dev/aterms2-run.ps1 <script.R> <log.txt>
param([string]$Script, [string]$Log)
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:TMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:TEMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/stan-cache"
$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$ErrorActionPreference = "Continue"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" $Script 2>&1 |
  ForEach-Object { "$_" } |
  Out-File -FilePath $Log -Encoding utf8
