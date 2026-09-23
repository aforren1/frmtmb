# R CMD check --as-cran on core, with vignettes and the manual.
# Adapted from dev/lane-rules.md: pandoc and TinyTeX must be on PATH or
# the run reports a bogus pdflatex ERROR, and --no-manual is NOT passed
# because the manual sections are where an Rd defect surfaces.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-skewinit"
$OUT  = "$ROOT/dev/skewinit-log/check"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"

if (-not (Test-Path $R)) { throw "R missing: $R" }
New-Item -ItemType Directory -Force $OUT | Out-Null

$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
# the check's own R must find the lane build's dependencies, and its
# own copy of frmtmb, without ever writing to a shared library
$env:R_LIBS = "C:/Users/adf44/source/r/skewinit-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"

Push-Location $ROOT
$ErrorActionPreference = "Continue"
& $R CMD build . 2>&1 | Tee-Object -FilePath "$OUT/build.log"
$tar = Get-ChildItem -Path $ROOT -Filter "frmtmb_*.tar.gz" |
       Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $tar) { Pop-Location; throw "no tarball built" }
& $R CMD check --as-cran -o $OUT $tar.FullName 2>&1 |
  Tee-Object -FilePath "$OUT/check.log"
Pop-Location
Write-Output ("CHECK DONE " + $tar.Name)
