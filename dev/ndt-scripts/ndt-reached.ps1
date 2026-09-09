# ndt lane, punch round 1: the test files this round's changes reach.
#
# NOT the whole suite. Eight of the nine test files edited in round 1
# are reverted to their 92e9330 state and the review has current
# numbers for them; what this runs is the files the punch round
# touched, plus the ones whose subject is the ndt link or the family
# objects, because the scalar bound moved back into the link.
param([string]$Lib = "C:/Users/adf44/source/r/ndt-lib")
$env:R_LIBS = "$Lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
Remove-Item Env:\FRMTMB_SCALE_TESTS -ErrorAction SilentlyContinue
Set-Location C:/Users/adf44/source/r/frmtmb-wt-ndt/extensions/frmtmb.eam
$files = @("test-ndt-bound.R", "test-family.R", "test-defects.R",
           "test-brms-parity.R", "test-variability.R", "test-lba.R",
           "test-rdm-gng.R", "test-simulate-density.R",
           "test-extension-api.R", "test-gddm-family.R",
           "test-bracket-access.R", "test-message-uniqueness.R",
           "test-surface.R", "test-sampling.R")
foreach ($nm in $files) {
  $expr = "suppressMessages({library(frmtmb);library(frmtmb.eam);library(testthat)}); testthat::set_max_fails(Inf); testthat::test_file('tests/testthat/$nm', package = 'frmtmb.eam')"
  $out = & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" --vanilla -e $expr 2>&1 | Out-String -Width 400
  $code = $LASTEXITCODE
  $sum = ($out -split "`n" | Select-String -Pattern "^\[ FAIL" | Select-Object -Last 1)
  if ($null -eq $sum) { $sum = "NO SUMMARY LINE" }
  Write-Output ("{0,-30} exit={1} {2}" -f $nm, $code, ($sum -replace "`r", ""))
  if ($sum -notmatch "FAIL 0 " -or $code -ne 0) {
    Write-Output "---- detail for $nm ----"
    Write-Output ($out -split "`n" | Select-String -Pattern "Failure|Error|FAILURE|ERROR" -Context 0,7 | Out-String -Width 400)
  }
}
Write-Output "==== REACHED FILES DONE ===="
