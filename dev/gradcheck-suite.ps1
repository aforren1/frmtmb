# Run the core test suite, ONE FILE PER R PROCESS, eight at a time.
# Counts come from dev/gradcheck-runtests.R, which reads `error` and
# `skipped` as well as `failed`.
#
#   powershell -File dev/gradcheck-suite.ps1 lane
#
# Writes one log per test file under dev/gradcheck-log/suite-<arm>/ so a
# file that aborts leaves its own evidence, and the driver's count comes
# from the RESULT lines rather than from the launches.

param([string]$Arm = "lane")

$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-gradcheck"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$out = Join-Path $root "dev\gradcheck-log\suite-$Arm"
New-Item -ItemType Directory -Force $out | Out-Null
Set-Location $root

$files = Get-ChildItem (Join-Path $root "tests\testthat") -Filter "test-*.R" |
  Sort-Object Name
Write-Output ("files: " + $files.Count)

$jobs = New-Object System.Collections.ArrayList
foreach ($f in $files) {
  $nm = $f.Name
  $log = Join-Path $out ($nm + ".txt")
  while (@(Get-Job -State Running).Count -ge 8) { Start-Sleep -Seconds 2 }
  $j = Start-Job -ScriptBlock {
    param($rs, $root, $rel, $log, $arm)
    Set-Location $root
    $env:NOT_CRAN = "true"
    & $rs "dev/gradcheck-runtests.R" $arm $rel *> $log
  } -ArgumentList $rs, $root, ("tests/testthat/" + $nm), $log, $Arm
  [void]$jobs.Add($j)
}
Wait-Job -Job $jobs | Out-Null
Remove-Job -Job $jobs -Force

$ran = 0
foreach ($f in $files) {
  $log = Join-Path $out ($f.Name + ".txt")
  if (-not (Test-Path $log)) {
    Write-Output ("NO LOG " + $f.Name)
    continue
  }
  $hit = Select-String -Path $log -Pattern "pass=|ABORTED" -SimpleMatch
  if ($hit) {
    $ran = $ran + 1
    foreach ($h in $hit) { Write-Output $h.Line }
  } else {
    Write-Output ("NO RESULT " + $f.Name)
    Get-Content $log -Tail 3 | ForEach-Object { Write-Output ("    " + $_) }
  }
}
Write-Output ("SUITE ran " + $ran + " of " + $files.Count)
