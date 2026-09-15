# Lane `learnhier`: run replicates under a HARD cap on how many R
# processes this lane holds at once.
#
#   dev/learnhier-launch2.ps1 -kind run -design rlddm -n 60 `
#       -outdir <dir> -cap 6
#   dev/learnhier-launch2.ps1 -kind imp -design bandit -n 40 -draws 100 `
#       -outdir <dir> -cap 6
#
# WHY THIS REPLACES learnhier-launch.ps1. That one counted the jobs IT
# started, which is not the same thing as the processes this lane is
# holding: a replicate started outside it, or a straggler it had already
# written off, is invisible to it, and the lane ran four concurrent
# rlddm fits while its launcher believed it was running three. The cap
# here is read from the machine every cycle, by matching `learnhier` in
# the command line of every live Rscript, so it is a cap on the REAL
# number.
#
# The unit is a PROCESS and not a job. `Rscript.exe` re-execs its own
# x64 build, so one replicate is two entries in the process table and a
# cap of 6 is three concurrent replicates.
#
# NEVER EDIT dev/learnhier-run.R WHILE THIS IS RUNNING. Rscript reads a
# script incrementally rather than parsing it whole, so a running
# replicate re-reads the file from a byte offset that an edit has
# invalidated. That cost this lane a 44-minute rlddm fit, which died at
# the end with `unexpected symbol in "rec$total_s <- as.numeric(difftim
# comparisons"`: the start of one line spliced to the tail of another.
param(
  [Parameter(Mandatory = $true)][string]$kind,
  [Parameter(Mandatory = $true)][string]$design,
  [Parameter(Mandatory = $true)][int]$n,
  [Parameter(Mandatory = $true)][string]$outdir,
  [int]$cap = 6,
  [int]$draws = 100,
  [int]$seedbase = 20261000,
  [int]$ns = 100,
  [int]$nt = 200,
  [int]$stagger = 60,
  [string]$runner = "dev/learnhier-run.R"
)
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$logdir = Join-Path $outdir "logs"
New-Item -ItemType Directory -Force $outdir | Out-Null
New-Item -ItemType Directory -Force $logdir | Out-Null
Set-Location $tree

function Get-LaneProcs {
  @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
      Where-Object { $_.CommandLine -like '*learnhier*' }).Count
}

# A seed whose replicate is ALREADY IN FLIGHT, started by an earlier
# launcher or by hand. The R side skips a seed whose result file exists,
# which is not the same test: a replicate that is running has produced
# no file yet, and starting a second copy of it wastes a slot on work
# that is already being done. This lane did exactly that once.
function Test-SeedRunning([int]$sd) {
  $hit = @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
             Where-Object { $_.CommandLine -like ("*learnhier*" + $sd + "*") })
  $hit.Count -gt 0
}

$queue = New-Object System.Collections.ArrayList
for ($i = 1; $i -le $n; $i++) { [void]$queue.Add($seedbase + 10 * $i) }

while ($queue.Count -gt 0) {
  $live = Get-LaneProcs
  if (($live + 2) -gt $cap) {
    Start-Sleep -Seconds 20
    continue
  }
  $sd = $queue[0]
  if (Test-SeedRunning $sd) {
    Start-Sleep -Seconds 20
    continue
  }
  $queue.RemoveAt(0)
  if ($kind -eq "run") {
    $tag = "learnhier-" + $design + "-" + $sd
    $argl = @($runner, $design, "$sd", $outdir, "$ns", "$nt")
  } else {
    $tag = "learnhier-imp-" + $design + "-" + $draws + "-" + $sd
    $argl = @("dev/learnhier-imp.R", $design, "$sd", "$draws", $outdir)
  }
  $p = Start-Process -FilePath $rs -PassThru -NoNewWindow `
    -WorkingDirectory $tree -ArgumentList $argl `
    -RedirectStandardOutput (Join-Path $logdir ($tag + ".log")) `
    -RedirectStandardError (Join-Path $logdir ($tag + ".err"))
  Write-Output ("START " + $tag + " pid=" + $p.Id + " lane_procs_before=" +
                $live + " left=" + $queue.Count + " at " +
                (Get-Date -Format "HH:mm:ss"))
  # apart, so that two tape builds do not peak together
  Start-Sleep -Seconds $stagger
}
# do not exit until the lane is quiet, so a caller can chain arms
while ((Get-LaneProcs) -gt 0) { Start-Sleep -Seconds 20 }
Write-Output ("ALL DONE " + $kind + " " + $design + " n=" + $n)
