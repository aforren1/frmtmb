# R CMD build + check --as-cran for one package, inside dev/formula2-check.
param([string]$pkg)
$wt = "C:\Users\adf44\source\r\frmtmb-wt-formula2"
$src = if ($pkg -eq "frmtmb") { $wt } else { "$wt\extensions\$pkg" }
$out = "$wt\dev\formula2-check\$pkg"
New-Item -ItemType Directory -Force $out | Out-Null
Set-Location $out
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:R_LIBS = "C:/Users/adf44/source/r/wt-formula2-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:R_LIBS_USER = $env:R_LIBS
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:NOT_CRAN = "false"
$R = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $R CMD build $src *> build.log
$tar = Get-ChildItem "$pkg*.tar.gz" | Select-Object -First 1
& $R CMD check --as-cran $tar.Name *> check.log
