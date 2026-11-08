# Lane bigx: build and check core with --as-cran inside dev/bigx-check.
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-bigx"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-bigx-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = "false"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
New-Item -ItemType Directory -Force "$wt/dev/bigx-check" | Out-Null
Set-Location "$wt/dev/bigx-check"
& $r CMD build $wt *> build.txt
$tb = Get-ChildItem -Filter "frmtmb_*.tar.gz" | Select-Object -First 1
& $r CMD check --as-cran $tb.Name *> check.txt
