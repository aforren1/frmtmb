# Lane wt-predfix copy of dev/release/run-suite.ps1: this worktree, the lane runner
# (lane library first), and each FAILED line kept under its RESULT line.
# The whole ungated suite, one R process per test file, all eight
# packages. This log is what dev/suite-baseline.tsv is regenerated
# from, so it must live somewhere that survives a session.

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-predfix"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/predfix-run-tests.R"
$LOG  = "$ROOT/dev/predfix-log/suite.log"

if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:NOT_CRAN = "true"
$env:PREDFIX_ARM = "lane"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue

$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add(@{ n = "frmtmb"; d = "$ROOT/tests/testthat" })
foreach ($n in @("frmtmb.eam", "frmtmb.sample", "frmtmb.coupling",
                 "frmtmb.spline", "frmtmb.learn", "frmtmb.latent",
                 "frmtmb.ode")) {
  [void]$jobs.Add(@{ n = $n; d = "$ROOT/extensions/$n/tests/testthat" })
}
foreach ($j in $jobs) {
  if (-not (Test-Path $j.d)) { throw ("test dir missing: " + $j.d) }
}

$ErrorActionPreference = "Continue"
$ran = 0
$want = 0
foreach ($j in $jobs) {
  $files = Get-ChildItem -Path $j.d -Filter "test-*.R" | Sort-Object Name
  $want = $want + $files.Count
  Add-Content -Path $LOG -Value ("== " + $j.n + " " + $files.Count + " files ==")
  foreach ($f in $files) {
    $out = & $R $S $j.n $f.FullName 2>&1
    $hit = $out | Select-String -Pattern "^RESULT"
    if ($hit) {
      $ran = $ran + 1
      Add-Content -Path $LOG -Value $hit.Line
    $bad = $out | Select-String -Pattern "^(FAILED|WARNING)"
    foreach ($b in $bad) { Add-Content -Path $LOG -Value ("  " + $b.Line) }
    } else {
      Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $f.FullName)
      Add-Content -Path $LOG -Value ($out | Select-Object -Last 5)
    }
  }
  Add-Content -Path $LOG -Value ("== " + $j.n + " END ==")
}
Add-Content -Path $LOG -Value ("SUITE ran " + $ran + " of " + $want)
Write-Output ("SUITE ran " + $ran + " of " + $want)
