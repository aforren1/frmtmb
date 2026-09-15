# Lane `learnhier`: roxygenise frmtmb.learn and install it into THIS
# lane's private library. Nothing else is ever installed: frmtmb, the
# eam sibling and the other extensions come from the round's shared
# reference build, which R_LIBS below makes visible to R CMD INSTALL as
# a READ path only.
#
# `$LASTEXITCODE`, NOT `$?`. Windows PowerShell sets `$?` to false the
# moment a native command writes anything to stderr, and roxygenise()
# writes "Loading required package: frmtmb" there on every run.
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$lib = "C:\Users\adf44\source\r\learnhier-lib"
$env:R_LIBS = "C:/Users/adf44/source/r/learnhier-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
Set-Location $tree
Rscript dev/learnhier-roxy.R
if ($LASTEXITCODE -eq 0) {
  & "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL --library=$lib --no-multiarch "$tree\extensions\frmtmb.learn"
  Write-Output "INSTALL exit=$LASTEXITCODE"
} else {
  Write-Output "ROXYGEN FAILED exit=$LASTEXITCODE; nothing installed"
}
