# One test file per R process, two at a time, over the whole
# frmtmb.sample suite. Counts come from the runner, which sums failures
# AND errors (dev/lane-rules.md).
$ErrorActionPreference = "Stop"
$root = "C:/Users/adf44/source/r/frmtmb-wt-brmsmatch"
$pkg = "$root/extensions/frmtmb.sample"
$log = "$root/dev/brmsmatch-suite-log"
New-Item -ItemType Directory -Force -Path $log | Out-Null
$env:FRMTMB_STAN_CACHE = "$root/dev/stan-cache"
$env:NOT_CRAN = "true"
$env:FRMTMB_BRMS_FIT_TESTS = "true"
$files = Get-ChildItem "$pkg/tests/testthat" -Filter "test-*.R" |
  Sort-Object Name
$rscript = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
$jobs = New-Object System.Collections.ArrayList
foreach ($f in $files) {
  $nm = $f.Name
  $full = $f.FullName -replace "\\", "/"
  while (@(Get-Job -State Running).Count -ge 2) { Start-Sleep -Seconds 2 }
  $j = Start-Job -ScriptBlock {
    param($rs, $runner, $full, $out, $root)
    Set-Location $root
    & $rs $runner "C:/Users/adf44/source/r/brmsmatch-lib" $full `
      "C:/Users/adf44/source/r/rellib-r3" *>&1 |
      Out-File -FilePath $out -Encoding utf8
  } -ArgumentList $rscript, "$root/dev/brmsmatch-runtest.R", $full,
      "$log/$nm.log", $root
  [void]$jobs.Add($j)
}
Get-Job | Wait-Job | Out-Null
Get-Job | Remove-Job
Write-Output "ALL FILES DONE"
