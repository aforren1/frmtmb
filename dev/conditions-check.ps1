# The lane's ONE R CMD check --as-cran per package, without --no-manual.
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR (dev/lane-rules.md).
#   powershell -File dev/conditions-check.ps1 core|sample
param([string]$which)
$ErrorActionPreference = "Stop"
$root = "C:/Users/adf44/source/r/frmtmb-wt-conditions"
$out = "C:/Users/adf44/source/r/conditions-check-$which"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/conditions-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:FRMTMB_STAN_CACHE = "$root/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
Set-Location $out
if ($which -eq "core") { $src = $root; $pat = "frmtmb_*.tar.gz" }
else { $src = "$root/extensions/frmtmb.sample"; $pat = "frmtmb.sample_*.tar.gz" }
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build $src
$tar = Get-ChildItem "$out/$pat" |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran `
  --output=$out $tar.FullName
Write-Output "CHECK DONE"
