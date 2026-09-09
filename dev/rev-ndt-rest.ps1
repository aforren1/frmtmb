# rev-ndt: the remaining measurements, serial, after the timing job.
$Rexe = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$Dir = "C:/Users/adf44/source/r/frmtmb-wt-ndt/dev"
Remove-Item Env:\R_LIBS -ErrorAction SilentlyContinue
Remove-Item Env:\FRMTMB_SCALE_TESTS -ErrorAction SilentlyContinue
$env:NOT_CRAN = "true"

Write-Output "==== idempotent $(Get-Date -Format o) ===="
& $Rexe --vanilla "$Dir/rev-ndt-idempotent.R" 2>&1 | Out-String -Width 400

foreach ($a in @("new-group", "old-plain", "new-group", "old-plain")) {
  $env:REV_ARM = $a
  Write-Output "==== iters $a $(Get-Date -Format o) ===="
  & $Rexe --vanilla "$Dir/rev-ndt-iters.R" 2>&1 | Out-String -Width 400
}
Remove-Item Env:\REV_ARM -ErrorAction SilentlyContinue

Write-Output "==== floors $(Get-Date -Format o) ===="
& $Rexe --vanilla "$Dir/rev-ndt-floors.R" 2>&1 | Out-String -Width 400
Write-Output "==== REST DONE $(Get-Date -Format o) ===="
