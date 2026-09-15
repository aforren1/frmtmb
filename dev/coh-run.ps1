# Runs dev/coh-recovery.R as several processes, one seed block each.
#
#   powershell -File dev/coh-run.ps1 -Arm main -Reps 200 -Workers 6
#
# One process per block, because a fit belongs in a process that has
# not already built a tape of its own, and because six of these fit in
# memory where one of them takes half an hour.
param(
  [string]$Arm = "main",
  [int]$Reps = 200,
  [int]$Workers = 6,
  [string]$Out = ""
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
if (-not $Out) { $Out = "dev/coh-recovery-$Arm.tsv" }

$rscript = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$per = [math]::Ceiling($Reps / $Workers)
$procs = New-Object System.Collections.ArrayList
for ($w = 0; $w -lt $Workers; $w++) {
  $from = $w * $per + 1
  $to = [math]::Min($Reps, ($w + 1) * $per)
  if ($from -gt $to) { continue }
  # Start-Process needs two files, and the per-replicate lines go to
  # STDERR because dev/coh-recovery.R reports them with message(): the
  # ".err" file is the one with the run in it.
  $log = "dev/coh-run-$Arm-$w-log.txt"
  $env:COH_ARM = $Arm
  $env:COH_FROM = "$from"
  $env:COH_TO = "$to"
  $env:COH_OUT = "$Out.$w"
  $env:NOT_CRAN = "true"
  $p = Start-Process -FilePath $rscript `
    -ArgumentList "dev/coh-recovery.R" `
    -RedirectStandardOutput $log -RedirectStandardError "$log.err" `
    -NoNewWindow -PassThru
  [void]$procs.Add($p)
  Write-Output "worker $w reps $from..$to pid $($p.Id)"
}
foreach ($p in $procs) { $p.WaitForExit() }
# One table out of the per-worker files, so that a summary reads one
# path. The parts are kept: a rerun of one block must not need all six.
$parts = Get-ChildItem -Path "$Out.*" | Sort-Object Name
$all = New-Object System.Collections.ArrayList
foreach ($f in $parts) {
  foreach ($ln in (Get-Content $f.FullName)) { [void]$all.Add($ln) }
}
[System.IO.File]::WriteAllText(
  (Join-Path $root $Out), ($all -join "`n") + "`n",
  (New-Object System.Text.UTF8Encoding($false)))
Write-Output "wrote $Out with $($all.Count) rows"
