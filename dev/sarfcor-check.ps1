# The lane's ONE R CMD check --as-cran per changed package, without
# --no-manual, built and checked inside dev/sarfcor-check/<pkg>.
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR (dev/lane-rules.md).
#   powershell -File dev/sarfcor-check.ps1 frmtmb|frmtmb.spatial|frmtmb.sample
param([string]$which)
$ErrorActionPreference = "Stop"
$root = "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
$out = "$root/dev/sarfcor-check/$which"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/wt-sarfcor-lib;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$root/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
Set-Location $out
if ($which -eq "frmtmb") { $src = $root } else { $src = "$root/extensions/$which" }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build $src
$tar = Get-ChildItem "$out/$($which)_*.tar.gz" |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$out $tar.FullName
Write-Output "CHECK DONE"
