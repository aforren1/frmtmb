# Environment for lane splinecurve: toolchain on PATH, private library first.
param([string]$Lib = "C:/Users/adf44/source/r/wt-splinecurve-lib")
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:R_LIBS = "$Lib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:R_LIBS_USER = "C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = "true"
