# Lane hollow: R CMD build and check --as-cran of one package, inside
# dev/hollow-check/. Usage: powershell -File dev/hollow-check.ps1 core|sample
param([string]$what = "core")
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-hollow"
$lib = "C:/Users/adf44/source/r/wt-hollow-lib"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
$src = if ($what -eq "core") { $wt } else { "$wt/extensions/frmtmb.sample" }
$out = "$wt/dev/hollow-check/$what"
New-Item -ItemType Directory -Force $out | Out-Null
Set-Location $out
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:R_LIBS = "$lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:R_LIBS_USER = ""
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$wt/dev/stan-cache"
Get-ChildItem -Filter "*.tar.gz" | Remove-Item -Force -Confirm:$false
& $r CMD build $src
Write-Output "BUILD EXIT $LASTEXITCODE"
$tar = Get-ChildItem -Filter "*.tar.gz" | Select-Object -First 1
& $r CMD check --as-cran $tar.Name
Write-Output "CHECK EXIT $LASTEXITCODE"
