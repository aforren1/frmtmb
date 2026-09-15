# Lane `learnhier`: run ONE replicate and watch what it costs the
# machine.
#
#   dev/learnhier-memwatch.ps1 -design rlddm -seed 20261010 -outdir <d> `
#                              -ns 100 -nt 200
#
# WHY THIS EXISTS. The Phase 0 scale row recorded a peak process of
# 7995 MB for the 100 by 200 rlddm fit, which is the largest figure in
# that table and was measured on a NON-CONVERGENT arm at frmtmb 0.55.0.
# How many of these can run at once is a memory budget, and other lanes
# are on this machine at the same time, so the budget is measured before
# it is spent rather than after.
#
# It samples the working set of the Rscript process and its x64 child,
# and the machine's committed bytes beside it, because the number that
# stops another lane is the COMMIT and not the resident set.
param(
  [Parameter(Mandatory = $true)][string]$design,
  [Parameter(Mandatory = $true)][int]$seed,
  [Parameter(Mandatory = $true)][string]$outdir,
  [int]$ns = 100,
  [int]$nt = 200,
  [int]$every = 5
)
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
New-Item -ItemType Directory -Force $outdir | Out-Null
$log = Join-Path $outdir ("learnhier-memwatch-" + $design + "-" + $seed + ".log")
$p = Start-Process -FilePath $rs -PassThru -NoNewWindow -WorkingDirectory $tree `
  -ArgumentList @("dev/learnhier-run.R", $design, "$seed", $outdir, "$ns", "$nt") `
  -RedirectStandardOutput (Join-Path $outdir ("learnhier-mw-" + $seed + ".out")) `
  -RedirectStandardError (Join-Path $outdir ("learnhier-mw-" + $seed + ".err"))
$peak = 0
$peakCommit = 0
$t0 = Get-Date
while (-not $p.HasExited) {
  Start-Sleep -Seconds $every
  $mine = Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
    Where-Object { $_.CommandLine -like ("*learnhier-run.R*" + $seed + "*") }
  $ws = 0
  foreach ($m in $mine) { $ws += $m.WorkingSetSize }
  $ws = [math]::Round($ws / 1MB)
  $cm = [math]::Round(((Get-Counter '\Memory\Committed Bytes').CounterSamples[0].CookedValue) / 1GB, 2)
  if ($ws -gt $peak) { $peak = $ws }
  if ($cm -gt $peakCommit) { $peakCommit = $cm }
  $el = [math]::Round(((Get-Date) - $t0).TotalSeconds)
  Add-Content -Path $log -Value ("$el s  ws=$ws MB  peak=$peak MB  commit=$cm GB")
}
Write-Output ("DONE " + $design + " seed=" + $seed + " exit=" + $p.ExitCode +
              " peak_ws_MB=" + $peak + " peak_commit_GB=" + $peakCommit +
              " elapsed_s=" + [math]::Round(((Get-Date) - $t0).TotalSeconds))
Get-Content (Join-Path $outdir ("learnhier-mw-" + $seed + ".out"))
