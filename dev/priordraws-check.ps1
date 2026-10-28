# R CMD check --as-cran of one package this lane changed, built WITH
# vignettes, in dev/priordraws-check/<name>, against the lane's library
# first. Modeled on dev/release/run-check.ps1.
#   powershell -File dev/priordraws-check.ps1 <package dir> <name>
param([string]$Pkg, [string]$Name)
$ErrorActionPreference = "Stop"
$R   = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$OUT = "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/priordraws-check/$Name"
$LOG = "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/priordraws-check/$Name.log"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-priordraws-lib;C:/Users/adf44/source/r/rellib-r14;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/stan-cache"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue
$ErrorActionPreference = "Continue"
Push-Location $OUT
& $R CMD build $Pkg 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
$tar = Get-ChildItem -Filter "*.tar.gz" | Select-Object -First 1
& $R CMD check --as-cran $tar.Name 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
Pop-Location
Add-Content -Path $LOG -Value "CHECK DONE"
