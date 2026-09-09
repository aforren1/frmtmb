# ndt lane: the three eam scale rows, one per fresh process, serial.
$env:R_LIBS = "C:/Users/adf44/source/r/ndt-lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:FRMTMB_SCALE_TESTS = "true"
$env:FRMTMB_SCALE_OUT = "C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/ndt-scripts/ndt-scale.tsv"
Remove-Item Env:\FRMTMB_SCALE_SMALL -ErrorAction SilentlyContinue
Set-Location C:/Users/adf44/source/r/frmtmb-wt-ndt/extensions/frmtmb.eam
foreach ($row in @("eam", "eam-unbounded", "eam-sv")) {
  $env:FRMTMB_SCALE_ROW = $row
  Write-Output "==== ROW $row  $(Get-Date -Format o) ===="
  & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" --vanilla -e "suppressMessages({library(frmtmb);library(frmtmb.eam);library(testthat)}); testthat::test_file('tests/testthat/test-scale.R')" 2>&1 | Out-String -Width 400
}
Write-Output "==== ALL ROWS DONE $(Get-Date -Format o) ===="
