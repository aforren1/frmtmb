# ndt lane: the whole frmtmb.eam suite, ONE TEST FILE PER PROCESS.
#
# `package =` is what puts the test env inside the package NAMESPACE.
# Without it a file that touches an internal (ddm_lpdf_both, ddm_floor,
# ddm_accepts) errors with "could not find function", which looks
# exactly like a real failure and is not.
#
# The reporter's default failure cap is lifted, because a capped count
# is not a count. Each line reports FAIL, WARN, SKIP and PASS as
# testthat gives them plus the process exit code, so a file that
# aborted halfway cannot print a clean line.
param([string]$Only = "")
$env:R_LIBS = "C:/Users/adf44/source/r/ndt-lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
Remove-Item Env:\FRMTMB_SCALE_TESTS -ErrorAction SilentlyContinue
Set-Location C:/Users/adf44/source/r/frmtmb-wt-ndt/extensions/frmtmb.eam
$files = Get-ChildItem tests/testthat -Filter "test-*.R" | Sort-Object Name
foreach ($f in $files) {
  $nm = $f.Name
  if ($Only -ne "" -and $nm -notlike $Only) { continue }
  $expr = "suppressMessages({library(frmtmb);library(frmtmb.eam);library(testthat)}); testthat::set_max_fails(Inf); testthat::test_file('tests/testthat/$nm', package = 'frmtmb.eam')"
  $out = & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" --vanilla -e $expr 2>&1 | Out-String -Width 400
  $code = $LASTEXITCODE
  $sum = ($out -split "`n" | Select-String -Pattern "^\[ FAIL" | Select-Object -Last 1)
  if ($null -eq $sum) { $sum = "NO SUMMARY LINE" }
  Write-Output ("{0,-32} exit={1} {2}" -f $nm, $code, ($sum -replace "`r", ""))
  if ($sum -notmatch "FAIL 0 " -or $code -ne 0) {
    Write-Output "---- detail for $nm ----"
    Write-Output ($out -split "`n" | Select-String -Pattern "Failure|Error|FAILURE|ERROR" -Context 0,7 | Out-String -Width 400)
  }
}
Write-Output "==== SUITE DONE ===="
