# Lane `learnhier`: R CMD check --as-cran on frmtmb.learn, ONCE.
#
# Built WITH vignettes and checked WITHOUT --no-manual, for the reasons
# dev/release/run-check.ps1 gives: --no-build-vignettes manufactures two
# WARNINGs and a NOTE, and --no-manual skips the HTML-manual check where
# this machine's environmental V8 math-rendering NOTE lives, which makes
# the result look cleaner than other lanes' and not comparable to them.
# The eamhier lane skipped the manual checks and shipped five unescaped
# `%` in an Rd table.
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$tree = "C:\Users\adf44\source\r\frmtmb-wt-learnhier"
$out  = "C:\Users\adf44\source\r\learnhier-check"
$env:R_LIBS = "C:/Users/adf44/source/r/learnhier-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
if (Test-Path $out) { Remove-Item -Recurse -Force $out }
New-Item -ItemType Directory -Force $out | Out-Null
Set-Location $out
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build "$tree\extensions\frmtmb.learn"
Write-Output "BUILD exit=$LASTEXITCODE"
$tar = (Get-ChildItem "$out\frmtmb.learn_*.tar.gz" | Select-Object -First 1).FullName
Write-Output "TARBALL $tar"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran $tar
Write-Output "CHECK exit=$LASTEXITCODE"
