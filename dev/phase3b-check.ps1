# R CMD build and R CMD check --as-cran, with vignettes and the manual,
# for one extension. The check tree is kept under the lane's check
# directory; its tests/testthat.Rout is copied into dev/phase3b-log
# before anything else touches it.
# Usage: powershell -File dev/phase3b-check.ps1 <pkg>
param([string]$pkg)
$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-phase3b"
$CHK  = "C:/Users/adf44/source/r/phase3b-check/$pkg"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$env:R_LIBS = "C:/Users/adf44/source/r/phase3b-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:R_LIBS_USER = "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
New-Item -ItemType Directory -Force $CHK | Out-Null
Set-Location $CHK
Get-ChildItem -Filter "$pkg*.tar.gz" -ErrorAction SilentlyContinue | Remove-Item
$ErrorActionPreference = "Continue"
& $R CMD build "$ROOT/extensions/$pkg" *> "$CHK/build.log"
$tar = Get-ChildItem -Filter "$pkg*.tar.gz" | Select-Object -First 1
if ($null -eq $tar) { "NO TARBALL" | Out-File "$CHK/check-status.txt"; exit 1 }
& $R CMD check --as-cran $tar.Name *> "$CHK/check.log"
$rout = "$CHK/$pkg.Rcheck/tests/testthat.Rout"
if (Test-Path $rout) { Copy-Item $rout "$ROOT/dev/phase3b-log/check-$pkg-testthat.Rout" }
if (Test-Path "$rout.fail") { Copy-Item "$rout.fail" "$ROOT/dev/phase3b-log/check-$pkg-testthat.Rout.fail" }
Copy-Item "$CHK/check.log" "$ROOT/dev/phase3b-log/check-$pkg.log"
Copy-Item "$CHK/build.log" "$ROOT/dev/phase3b-log/build-$pkg.log"
