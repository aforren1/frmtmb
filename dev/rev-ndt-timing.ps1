# rev-ndt: interleaved arms, one process each, three rounds.
#
# PowerShell variables are case insensitive, so a loop counter named $r
# overwrites an executable path named $R. That silently turned every
# invocation into `& 3` on the first attempt.
$Rexe = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$Script = "C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/rev-ndt-timing.R"
Remove-Item Env:\R_LIBS -ErrorAction SilentlyContinue
Remove-Item Env:\FRMTMB_SCALE_TESTS -ErrorAction SilentlyContinue
$env:NOT_CRAN = "true"
$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add("new-plain")
[void]$jobs.Add("new-unbounded")
foreach ($round in 1..3) {
  [void]$jobs.Add("new-group")
  [void]$jobs.Add("old-plain")
}
foreach ($a in $jobs) {
  $env:REV_ARM = $a
  Write-Output "==== ARM $a  $(Get-Date -Format o) ===="
  & $Rexe --vanilla $Script 2>&1 | Out-String -Width 500
}
Write-Output "==== TIMING DONE $(Get-Date -Format o) ===="
