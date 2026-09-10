# Lane `latent`: roxygenise frmtmb.latent and install it into THIS
# lane's private library. Nothing else is ever installed: frmtmb and the
# other extensions come from the round's shared reference build, which
# R_LIBS below makes visible to R CMD INSTALL as a READ path only.
#
# `$LASTEXITCODE`, NOT `$?`. Windows PowerShell sets `$?` to false the
# moment a native command writes anything to stderr, and roxygenise()
# writes "Loading required package: frmtmb" there on every run. Gating
# the install on `$?` therefore skipped the install every time while
# reporting nothing, and two probes were then run against a stale
# library.
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$lib = "C:\Users\adf44\source\r\latent-lib"
$env:R_LIBS = "C:/Users/adf44/source/r/latent-lib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
Set-Location $tree
Rscript dev/latent-roxy.R
if ($LASTEXITCODE -eq 0) {
  & "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch "$tree\extensions\frmtmb.latent"
  Write-Output "INSTALL exit=$LASTEXITCODE"
} else {
  Write-Output "ROXYGEN FAILED exit=$LASTEXITCODE; nothing installed"
}
