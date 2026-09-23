# Which frmtmb.sample test files actually leave a Rplots.pdf behind?
#
# WHY THIS IS MEASURED. Punch round 1 named six files from a grep for
# drawing calls. A grep cannot tell a call that runs from one that does
# not: a `pp_check()` inside `quote()` never evaluates, and one inside
# `expect_error()` may throw before it draws. The list is going to be
# quoted when the device-guard fix is taken, so it has to come from a
# run.
#
# The gated environment is ON, because that is the superset: a drawing
# call inside a gated block draws only there, and `R CMD build` usually
# follows a full run. One file per process, as this project requires,
# and the file is deleted BEFORE each one so a leftover cannot be read
# as a hit.
#
#   powershell -File dev/correct-punch4-rplots.ps1

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-correct"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/correct-run-tests.R"
$DIR  = "$ROOT/extensions/frmtmb.sample/tests/testthat"
$PDF  = "$DIR/Rplots.pdf"
$LOG  = "$ROOT/dev/correct-log/punch4-rplots.txt"

if (-not (Test-Path $R)) { throw "Rscript missing: $R" }
if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $DIR)) { throw "test dir missing: $DIR" }

$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue
Add-Content -Path $LOG -Value ("CORRECT_LIB=" + $env:CORRECT_LIB)

$files = Get-ChildItem -Path $DIR -Filter "test-*.R" | Sort-Object Name
$ErrorActionPreference = "Continue"
$hits = 0
foreach ($f in $files) {
  $nm = $f.Name
  Remove-Item $PDF -Force -ErrorAction SilentlyContinue
  if (Test-Path $PDF) { throw "could not clear $PDF before $nm" }
  $out = & $R $S "frmtmb.sample" $f.FullName 2>&1
  $res = ($out | Select-String -Pattern "^RESULT").Line
  if (-not $res) { $res = "NO RESULT LINE" }
  if (Test-Path $PDF) {
    $hits = $hits + 1
    Add-Content -Path $LOG -Value ("DRAWS   " + $nm + "  " + $res)
  } else {
    Add-Content -Path $LOG -Value ("clean   " + $nm + "  " + $res)
  }
}
Remove-Item $PDF -Force -ErrorAction SilentlyContinue
Add-Content -Path $LOG -Value ("FILES " + $files.Count + " DRAWS " + $hits)
Write-Output ("FILES " + $files.Count + " DRAWS " + $hits)
