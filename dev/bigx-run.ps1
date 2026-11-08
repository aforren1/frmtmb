# Run one of lane bigx's R scripts with the lane's environment.
#   powershell -File dev/bigx-run.ps1 <arm> <script> <log> [args...]
param([string]$arm, [string]$script, [string]$log)
$env:PATH = $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:BIGX_ARM = $arm
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/stan-cache"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" $script @args *> $log
