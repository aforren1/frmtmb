# Local validation of dev/release/build-docs.R on the reference library.
# Scratch only: builds into dev/docsci-site/, never into the committed
# docs/ tree.
$ErrorActionPreference = "Continue"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-docsci"
$R    = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"

$env:R_LIBS = "C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\rtools45\usr\bin;" + $env:PATH

& $R "$ROOT/dev/release/build-docs.R" `
  "--dest=$ROOT/dev/docsci-site" `
  "--log=$ROOT/dev/docsci-build.log"
Write-Output ("EXIT " + $LASTEXITCODE)
