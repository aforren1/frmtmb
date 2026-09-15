# Lane `learnhier`: run N importance-correction replicates, one per
# process, `par` at a time.
#
#   dev/learnhier-launch-imp.ps1 -design bandit -n 40 -draws 100 `
#                                -outdir <dir> -par 3
#
# The seeds are the SAME seeds the recovery arm used, so a replicate's
# correction and its Laplace recovery row are the same dataset.
param(
  [Parameter(Mandatory = $true)][string]$design,
  [Parameter(Mandatory = $true)][int]$n,
  [Parameter(Mandatory = $true)][int]$draws,
  [Parameter(Mandatory = $true)][string]$outdir,
  [int]$par = 3,
  [int]$seedbase = 20261000
)
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$logdir = Join-Path $outdir "logs"
New-Item -ItemType Directory -Force $outdir | Out-Null
New-Item -ItemType Directory -Force $logdir | Out-Null

$jobs = New-Object System.Collections.ArrayList
$queue = New-Object System.Collections.ArrayList
for ($i = 1; $i -le $n; $i++) { [void]$queue.Add($seedbase + 10 * $i) }

$done = 0
Set-Location $tree
while ($done -lt $n) {
  while (($jobs.Count -lt $par) -and ($queue.Count -gt 0)) {
    $sd = $queue[0]
    $queue.RemoveAt(0)
    $tag = "learnhier-imp-" + $design + "-" + $draws + "-" + $sd
    $p = Start-Process -FilePath $rs -PassThru -NoNewWindow `
      -WorkingDirectory $tree `
      -ArgumentList @("dev/learnhier-imp.R", $design, "$sd", "$draws",
                      $outdir) `
      -RedirectStandardOutput (Join-Path $logdir ($tag + ".log")) `
      -RedirectStandardError (Join-Path $logdir ($tag + ".err"))
    [void]$jobs.Add($p)
    Write-Output ("START imp seed=" + $sd + " pid=" + $p.Id)
  }
  Start-Sleep -Seconds 10
  $live = New-Object System.Collections.ArrayList
  foreach ($j in $jobs) {
    if ($j.HasExited) {
      $done++
      Write-Output ("EXIT  pid=" + $j.Id + " code=" + $j.ExitCode +
                    " done=" + $done + "/" + $n + " at " +
                    (Get-Date -Format "HH:mm:ss"))
    } else { [void]$live.Add($j) }
  }
  $jobs = $live
}
Write-Output ("ALL DONE imp " + $design + " draws=" + $draws + " n=" + $n)
