# Reviewer driver: run a TSV job list, N at a time, one R process per job.
# Each line: <tag> <TAB> <core-lib> <TAB> <test-file> [<TAB> <pkg>]
# Logs land in dev/gradcheck-rev-log/<tag>.txt
param(
  [Parameter(Mandatory = $true)][string]$JobList,
  [int]$Par = 3
)
$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-gradcheck"
$logd = Join-Path $root "dev\gradcheck-rev-log"
if (-not (Test-Path $logd)) { New-Item -ItemType Directory $logd | Out-Null }
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$runner = Join-Path $root "dev\gradcheck-rev-runfile.R"

$jobs = New-Object System.Collections.ArrayList
foreach ($line in (Get-Content $JobList)) {
  if ($line.Trim().Length -eq 0) { continue }
  if ($line.StartsWith("#")) { continue }
  [void]$jobs.Add(($line -split "`t"))
}
Write-Host ("jobs: " + $jobs.Count)

$running = New-Object System.Collections.ArrayList
$done = 0
foreach ($j in $jobs) {
  while ($running.Count -ge $Par) {
    Start-Sleep -Seconds 2
    $still = New-Object System.Collections.ArrayList
    foreach ($p in $running) {
      if (-not $p.HasExited) { [void]$still.Add($p) } else { $done++ }
    }
    $running = $still
  }
  $tag = $j[0]
  $log = Join-Path $logd ($tag + ".txt")
  $args = @($runner, $j[1], $j[2])
  if ($j.Count -ge 4) { $args += $j[3] }
  Write-Host ("start " + $tag)
  $p = Start-Process -FilePath $rs -ArgumentList $args `
    -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") `
    -NoNewWindow -PassThru
  [void]$running.Add($p)
}
foreach ($p in $running) { $p.WaitForExit() }
Write-Host "ALL DONE"
foreach ($j in $jobs) {
  $log = Join-Path $logd ($j[0] + ".txt")
  $hit = Select-String -Path $log -Pattern "^RESULT" -ErrorAction SilentlyContinue
  if ($hit) { Write-Host ($j[0] + "  " + $hit.Line) }
  else { Write-Host ($j[0] + "  NO RESULT LINE") }
}
