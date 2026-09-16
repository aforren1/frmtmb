# The lane's ONE R CMD check --as-cran, without --no-manual. pandoc and
# TinyTeX go on PATH or the check reports a bogus pdflatex ERROR
# (dev/lane-rules.md).
$ErrorActionPreference = "Stop"
$root = "C:/Users/adf44/source/r/frmtmb-wt-brmsmatch"
$out = "C:/Users/adf44/source/r/brmsmatch-check"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/brmsmatch-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "$root/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
Set-Location $out
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build "$root/extensions/frmtmb.sample"
$tar = Get-ChildItem "$out/frmtmb.sample_*.tar.gz" |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$out $tar.FullName
Write-Output "CHECK DONE"
