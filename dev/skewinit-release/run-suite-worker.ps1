# One worker of the parallel ungated run. Takes every Nth file of the
# shared job list, still one R process per file, and writes its own
# log. Parameters: <slot 0-based> <stride>.
#
# The serial driver dev/release/run-suite.ps1 is unchanged and is what
# a release uses. This exists only because three other lanes were
# running on the same box and the serial pass was measured at about
# two minutes per file, which is nine hours for 270 files.

param([int]$Slot, [int]$Stride)

$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-skewinit"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$S    = "$ROOT/dev/skewinit-release/run-tests.R"
$LOG  = "$ROOT/dev/skewinit-log/suite-$Slot.log"

if (-not (Test-Path $S)) { throw "runner missing: $S" }
if (-not (Test-Path $R)) { throw "Rscript missing: $R" }

$env:NOT_CRAN = "true"
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

# One flat, deterministic list, so slot k and stride N partition it
$all = New-Object System.Collections.ArrayList
foreach ($j in $jobs) {
  $files = Get-ChildItem -Path $j.d -Filter "test-*.R" | Sort-Object Name
  foreach ($f in $files) {
    [void]$all.Add(@{ n = $j.n; f = $f.FullName })
  }
}

$ErrorActionPreference = "Continue"
$ran = 0
$want = 0
for ($i = $Slot; $i -lt $all.Count; $i = $i + $Stride) {
  $want = $want + 1
  $it = $all[$i]
  $out = & $R $S $it.n $it.f 2>&1
  $hit = $out | Select-String -Pattern "^RESULT"
  if ($hit) {
    $ran = $ran + 1
    Add-Content -Path $LOG -Value ($it.n + " " + $hit.Line)
  } else {
    Add-Content -Path $LOG -Value ("NO RESULT LINE for " + $it.f)
    Add-Content -Path $LOG -Value ($out | Select-Object -Last 5)
  }
}
Add-Content -Path $LOG -Value ("WORKER $Slot ran " + $ran + " of " + $want)
Write-Output ("WORKER $Slot ran " + $ran + " of " + $want)
