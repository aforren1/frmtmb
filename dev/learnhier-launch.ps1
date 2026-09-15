# Lane `learnhier`: run N replicates of one design, one per PROCESS,
# `par` of them at a time.
#
#   dev/learnhier-launch.ps1 -design rlddm -n 60 -outdir <dir> -par 2
#
# ONE REPLICATE PER PROCESS. A 100 by 200 rlddm fit peaked at 8.0 GB of
# working set in the Phase 0 scale row, so `par` is a memory budget and
# not a core count: this machine has 12 physical cores and 31.7 GB, and
# the cores are not the binding constraint.
#
# Seeds are spaced by 10 because the simulators offset the seed by 1 and
# by 2 internally, so consecutive seeds would make one replicate's task
# design share a stream with the next replicate's deviations.
#
# `$LASTEXITCODE`, not `$?`: R writes package-loading messages to stderr
# and `$?` goes false on that alone.
param(
  [Parameter(Mandatory = $true)][string]$design,
  [Parameter(Mandatory = $true)][int]$n,
  [Parameter(Mandatory = $true)][string]$outdir,
  [int]$par = 4,
  [int]$seedbase = 20261000,
  [int]$ns = 100,
  [int]$nt = 200,
  [int]$stagger = 0
)

$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$logdir = Join-Path $outdir "logs"
New-Item -ItemType Directory -Force $outdir | Out-Null
New-Item -ItemType Directory -Force $logdir | Out-Null

# An ArrayList, because PowerShell flattens @(@(a,b)) into @(a,b) and a
# job list built with += has silently lost a whole suite on this project
# before.
$jobs = New-Object System.Collections.ArrayList
$queue = New-Object System.Collections.ArrayList
for ($i = 1; $i -le $n; $i++) { [void]$queue.Add($seedbase + 10 * $i) }

$done = 0
$started = 0
Set-Location $tree
while ($done -lt $n) {
  while (($jobs.Count -lt $par) -and ($queue.Count -gt 0)) {
    $sd = $queue[0]
    $queue.RemoveAt(0)
    $lo = Join-Path $logdir ("learnhier-" + $design + "-" + $sd + ".log")
    $le = Join-Path $logdir ("learnhier-" + $design + "-" + $sd + ".err")
    $p = Start-Process -FilePath $rs -PassThru -NoNewWindow `
      -WorkingDirectory $tree `
      -ArgumentList @("dev/learnhier-run.R", $design, "$sd", $outdir,
                      "$ns", "$nt") `
      -RedirectStandardOutput $lo -RedirectStandardError $le
    [void]$jobs.Add($p)
    $started++
    Write-Output ("START " + $design + " seed=" + $sd + " pid=" + $p.Id +
                  " running=" + $jobs.Count + " started=" + $started +
                  "/" + $n + " at " + (Get-Date -Format "HH:mm:ss"))
    # The peak of an rlddm fit is a TRANSIENT during the tape build, not
    # the steady state: 4485 MB against about 2500 MB once the tape is
    # up. Starting concurrent jobs apart keeps those transients apart,
    # which is what keeps this lane from taking the machine's commit
    # away from another one.
    if ($stagger -gt 0) { Start-Sleep -Seconds $stagger }
  }
  Start-Sleep -Seconds 5
  $live = New-Object System.Collections.ArrayList
  foreach ($j in $jobs) {
    if ($j.HasExited) {
      $done++
      Write-Output ("EXIT  pid=" + $j.Id + " code=" + $j.ExitCode +
                    " done=" + $done + "/" + $n + " at " +
                    (Get-Date -Format "HH:mm:ss"))
    } else {
      [void]$live.Add($j)
    }
  }
  $jobs = $live
}
Write-Output ("ALL DONE " + $design + " n=" + $n)
