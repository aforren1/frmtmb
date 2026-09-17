# The frmtmb.sample suite, one R process per test file, in one arm.
# Usage: powershell -File tmbstan121-05-suite.ps1 <arm>
param([string]$arm)
$ErrorActionPreference = "Stop"
$WT  = "C:/Users/adf44/source/r/frmtmb-wt-tmbstan121"
$R   = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S   = "$WT/dev/tmbstan121-05-run-tests.R"
$LOG = "$WT/dev/tmbstan121-suite-$arm.log"
$DIR = "$WT/extensions/frmtmb.sample/tests/testthat"
if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $DIR)) { throw "test dir missing: $DIR" }

$env:NOT_CRAN = "true"
# Each arm reads its own copy of the cache, so a program one arm
# compiles cannot be served to the other and hide a compile failure.
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/tmbstan121-lib/stan-cache-$arm"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue
$ErrorActionPreference = "Continue"
$files = Get-ChildItem -Path $DIR -Filter "test-*.R" | Sort-Object Name
Add-Content -Path $LOG -Value ("== frmtmb.sample " + $files.Count + " files arm " + $arm + " ==")
$ran = 0
foreach ($f in $files) {
  $t0 = Get-Date
  $out = & $R $S $arm "frmtmb.sample" $f.FullName 2>&1
  $sec = [math]::Round(((Get-Date) - $t0).TotalSeconds)
  $hit = $out | Select-String -Pattern "^RESULT"
  if ($hit) {
    $ran = $ran + 1
    Add-Content -Path $LOG -Value ($hit.Line + " secs=" + $sec)
    $det = $out | Select-String -Pattern "^DETAIL"
    foreach ($x in $det) { Add-Content -Path $LOG -Value $x.Line }
  } else {
    Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $f.FullName)
    Add-Content -Path $LOG -Value ($out | Select-Object -Last 8)
  }
}
Add-Content -Path $LOG -Value ("== frmtmb.sample END ==")
Add-Content -Path $LOG -Value ("SUITE ran " + $ran + " of " + $files.Count)
