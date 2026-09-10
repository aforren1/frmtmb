# The gated tiers: brms fit tests, the BCM chapters, the RL identity,
# the fuzz grammar tier, and the two sampler files. One R process per
# file. Reports how many files actually produced a RESULT line, because
# a driver that prints a completion marker over zero work has happened
# three times in this project.

# "Stop" only while checking preconditions. It must NOT be in force
# around the R calls: PowerShell 5.1 wraps every stderr line from a
# native command in an ErrorRecord, so an ordinary rstan message such
# as "the number of chains is less than 1" becomes terminating and
# kills the run 13 files in. See dev/lane-rules.md, shell traps.
$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/release/run-tests.R"
$LOG  = "$ROOT/dev/release/gated.log"

if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$env:FRMTMB_FUZZ = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
# backslashes: a forward-slash PATH left `cmd` unresolvable once, and
# the driver ran nothing while still writing its completion marker
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Remove-Item $LOG -ErrorAction SilentlyContinue
$core = "$ROOT/tests/testthat"

$jobs = New-Object System.Collections.ArrayList
foreach ($f in (Get-ChildItem -Path $core -Filter "test-bcm-*.R" | Sort-Object Name)) {
  [void]$jobs.Add(@{ n = "frmtmb"; f = $f.FullName })
}
foreach ($nm in @("test-brms-agreement.R", "test-brms-likelihood.R",
                  "test-brms-methods.R", "test-brms-priors.R",
                  "test-brms-port.R", "test-rl-example.R",
                  "test-fuzz.R")) {
  [void]$jobs.Add(@{ n = "frmtmb"; f = "$core/$nm" })
}
[void]$jobs.Add(@{ n = "frmtmb.learn"; f = "$ROOT/extensions/frmtmb.learn/tests/testthat/test-stan-identity.R" })
[void]$jobs.Add(@{ n = "frmtmb.sample"; f = "$ROOT/extensions/frmtmb.sample/tests/testthat/test-loo.R" })
[void]$jobs.Add(@{ n = "frmtmb.sample"; f = "$ROOT/extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R" })

foreach ($j in $jobs) {
  if (-not (Test-Path $j.f)) { throw ("test file missing: " + $j.f) }
}

Add-Content -Path $LOG -Value ("== gated " + $jobs.Count + " files ==")
$ErrorActionPreference = "Continue"
$ran = 0
foreach ($j in $jobs) {
  $out = & $R $S $j.n $j.f 2>&1
  $hit = $out | Select-String -Pattern "^RESULT"
  if ($hit) {
    $ran = $ran + 1
    Add-Content -Path $LOG -Value $hit.Line
  } else {
    Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $j.f)
    Add-Content -Path $LOG -Value ($out | Select-Object -Last 5)
  }
}
Add-Content -Path $LOG -Value ("GATED ran " + $ran + " of " + $jobs.Count)
Write-Output ("GATED ran " + $ran + " of " + $jobs.Count)
