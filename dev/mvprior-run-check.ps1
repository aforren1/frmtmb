# Lane wt-mvprior: dev/release/run-check.ps1 on the four packages this
# lane touched (frmtmb, frmtmb.sample, frmtmb.latent, frmtmb.eam), with
# the lane library first.
# Built WITH vignettes and checked WITH the manual, as the release is.
# Each package's tests/testthat.Rout is copied into dev/mvprior-log
# BEFORE the check tree can be deleted, because a passing check echoes
# no test counts into its log. The destination is cleared first, so a
# stale copy cannot stand in for this run's.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-mvprior"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$LOG  = "$ROOT/dev/mvprior-log/check.log"
$OUT  = "C:/Users/adf44/source/r/mvprior-check"

if (-not (Test-Path $R)) { throw "R missing: $R" }

$env:R_LIBS = "C:/Users/adf44/source/r/mvprior-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue

$pkgs = New-Object System.Collections.ArrayList
[void]$pkgs.Add(@{ d = $ROOT; n = "frmtmb" })
[void]$pkgs.Add(@{ d = "$ROOT/extensions/frmtmb.sample"; n = "frmtmb.sample" })
[void]$pkgs.Add(@{ d = "$ROOT/extensions/frmtmb.latent"; n = "frmtmb.latent" })
[void]$pkgs.Add(@{ d = "$ROOT/extensions/frmtmb.eam"; n = "frmtmb.eam" })

$ErrorActionPreference = "Continue"
$done = 0
foreach ($p in $pkgs) {
  $nm = $p.n
  Add-Content -Path $LOG -Value ("===== CHECK " + $nm + " =====")
  Push-Location $OUT
  Get-ChildItem -Path $OUT -Filter "*.tar.gz" | Remove-Item
  & $R CMD build $p.d 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
  $tgz = Get-ChildItem -Path $OUT -Filter "$nm`_*.tar.gz" | Select-Object -Last 1
  if ($tgz) {
    & $R CMD check --as-cran $tgz.FullName 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
    $done = $done + 1
    $dst = "$ROOT/dev/mvprior-log/$nm-testthat.Rout"
    Remove-Item $dst -ErrorAction SilentlyContinue
    $src = "$OUT/$nm.Rcheck/tests/testthat.Rout"
    if (Test-Path $src) {
      Copy-Item $src $dst
    } else {
      $fail = "$OUT/$nm.Rcheck/tests/testthat.Rout.fail"
      if (Test-Path $fail) { Copy-Item $fail $dst } else {
        Add-Content -Path $LOG -Value ("NO testthat.Rout for " + $nm)
      }
    }
  } else {
    Add-Content -Path $LOG -Value "NO TARBALL BUILT"
  }
  Pop-Location
}
Add-Content -Path $LOG -Value ("CHECKED " + $done + " of " + $pkgs.Count)
Write-Output ("CHECKED " + $done + " of " + $pkgs.Count)
