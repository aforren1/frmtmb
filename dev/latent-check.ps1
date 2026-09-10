# Lane `latent`: `R CMD check --as-cran` on frmtmb.latent, ONCE, on the
# final pass.
#
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR (dev/lane-rules.md). No --no-build-vignettes and no
# --no-manual: an earlier round manufactured two WARNINGs and a NOTE on
# all eight packages with the first, and skipped the environmental V8
# note with the second.
#
#   powershell -File dev/latent-check.ps1
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\R\R-4.6.1\bin;C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/latent-lib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$tree = "C:\Users\adf44\source\r\frmtmb-wt-latent"
$out = "C:\Users\adf44\source\r\latent-check"
New-Item -ItemType Directory -Force $out | Out-Null
Set-Location $out
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build "$tree\extensions\frmtmb.latent"
Write-Output "BUILD exit=$LASTEXITCODE"
if ($LASTEXITCODE -ne 0) { exit 1 }
$tgz = (Get-ChildItem "$out\frmtmb.latent_*.tar.gz" | Sort-Object LastWriteTime | Select-Object -Last 1).Name
Write-Output "TARBALL $tgz"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran --no-multiarch $tgz
Write-Output "CHECK exit=$LASTEXITCODE"
Get-Content "$out\frmtmb.latent.Rcheck\00check.log" | Select-String -Pattern "NOTE|WARNING|ERROR|OK$" | Select-Object -Last 40
