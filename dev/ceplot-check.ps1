# Lane ceplot: R CMD build and R CMD check --as-cran of one package of
# the lane, built and checked inside dev/ceplot-check/<pkg>, once.
#   powershell -File dev/ceplot-check.ps1 frmtmb
#   powershell -File dev/ceplot-check.ps1 frmtmb.sample
param([string]$pkg = "frmtmb")
$ErrorActionPreference = "Continue"
$WT = "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
$R = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$src = if ($pkg -eq "frmtmb") { $WT } else { "$WT/extensions/$pkg" }
$out = "$WT/dev/ceplot-check/$pkg"
New-Item -ItemType Directory -Force $out | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/wt-ceplot-lib;C:/Users/adf44/source/r/rellib-r4;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
Remove-Item Env:NOT_CRAN -ErrorAction SilentlyContinue
Remove-Item Env:FRMTMB_BRMS_FIT_TESTS -ErrorAction SilentlyContinue
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
Set-Location $out
Get-ChildItem -Filter "$pkg*.tar.gz" | Remove-Item -Force
& $R CMD build $src *> "$out/build.log"
$tar = Get-ChildItem -Filter "$pkg*.tar.gz" | Select-Object -First 1
if ($null -eq $tar) { "BUILD FAILED" | Out-File "$out/status.txt"; exit 1 }
& $R CMD check --as-cran $tar.Name *> "$out/check.log"
Select-String -Path "$out/check.log" -Pattern "^Status:" | ForEach-Object { $_.Line } | Out-File "$out/status.txt"
