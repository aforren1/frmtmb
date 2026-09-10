# The scale tier: one realistic-scale fit per extension, gated behind
# FRMTMB_SCALE_TESTS. These rows are MEASUREMENTS, so nothing else
# should be running while this does.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/release/run-tests.R"
$LOG  = "$ROOT/dev/release/scale.log"

if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:NOT_CRAN = "true"
$env:FRMTMB_SCALE_TESTS = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue

$jobs = New-Object System.Collections.ArrayList
foreach ($n in @("frmtmb.eam", "frmtmb.ode", "frmtmb.latent",
                 "frmtmb.spline", "frmtmb.coupling", "frmtmb.sample",
                 "frmtmb.learn")) {
  [void]$jobs.Add(@{ n = $n; f = "$ROOT/extensions/$n/tests/testthat/test-scale.R" })
}
foreach ($j in $jobs) {
  if (-not (Test-Path $j.f)) { throw ("test file missing: " + $j.f) }
}

Add-Content -Path $LOG -Value ("== scale " + $jobs.Count + " files ==")
# stderr from a native command is an ErrorRecord under "Stop" in
# PowerShell 5.1, and an ordinary rstan message would then kill the run
$ErrorActionPreference = "Continue"
$ran = 0
foreach ($j in $jobs) {
  Add-Content -Path $LOG -Value ("-- " + $j.n + " start " + (Get-Date -Format "HH:mm:ss"))
  $out = & $R $S $j.n $j.f 2>&1
  $hit = $out | Select-String -Pattern "^RESULT"
  $rows = $out | Select-String -Pattern "^SCALE row="
  if ($rows) { Add-Content -Path $LOG -Value $rows.Line }
  if ($hit) {
    $ran = $ran + 1
    Add-Content -Path $LOG -Value $hit.Line
  } else {
    Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $j.f)
    Add-Content -Path $LOG -Value ($out | Select-Object -Last 5)
  }
}
Add-Content -Path $LOG -Value ("SCALE ran " + $ran + " of " + $jobs.Count)
Write-Output ("SCALE ran " + $ran + " of " + $jobs.Count)
