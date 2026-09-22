# The lane's suite driver: one R process per test file, against the
# lane's private library. Modeled on dev/release/run-suite.ps1.
#   powershell -File dev/shapes-suite.ps1 [-Only frmtmb] [-Gated]

param([string]$Only = "", [switch]$Gated, [string]$Tag = "suite")

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-shapes"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/shapes-run.R"
$LOG  = "$ROOT/dev/shapes-log/$Tag.log"

if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "$ROOT/dev/stan-cache"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
if ($Gated) { $env:FRMTMB_BRMS_FIT_TESTS = "true" }
else { Remove-Item Env:\FRMTMB_BRMS_FIT_TESTS -ErrorAction SilentlyContinue }

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
  if ($Only -ne "" -and $j.n -ne $Only) { continue }
  $files = Get-ChildItem -Path $j.d -Filter "test-*.R" | Sort-Object Name
  $want = $want + $files.Count
  Add-Content -Path $LOG -Value ("== " + $j.n + " " + $files.Count + " files ==")
  foreach ($f in $files) {
    $out = & $R $S $j.n $f.FullName 2>&1
    $hit = $out | Select-String -Pattern "^RESULT"
    if ($hit) {
      $ran = $ran + 1
      Add-Content -Path $LOG -Value ($j.n + " " + $hit.Line)
      $bad = $out | Select-String -Pattern "^  "
      if ($bad) { Add-Content -Path $LOG -Value $bad.Line }
    } else {
      Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $f.FullName)
      Add-Content -Path $LOG -Value ($out | Select-Object -Last 8)
    }
  }
  Add-Content -Path $LOG -Value ("== " + $j.n + " END ==")
}
Add-Content -Path $LOG -Value ("RAN " + $ran + " of " + $want + " files")
Write-Output ("RAN " + $ran + " of " + $want + " files; log " + $LOG)
