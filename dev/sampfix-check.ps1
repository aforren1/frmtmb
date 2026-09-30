# Lane sampfix: R CMD build and check --as-cran of ONE package, inside
# dev/sampfix-check/<name>, against the private library first.
#   powershell -File dev/sampfix-check.ps1 <package path> <name>
param([string]$Pkg, [string]$Name)
$ErrorActionPreference = "Stop"
$WT  = "C:/Users/adf44/source/r/frmtmb-wt-sampfix"
$R   = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$OUT = "$WT/dev/sampfix-check/$Name"
$LOG = "$WT/dev/sampfix-check/$Name.log"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-sampfix-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:FRMTMB_STAN_CACHE = "$WT/dev/stan-cache"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
Remove-Item $LOG -ErrorAction SilentlyContinue
$ErrorActionPreference = "Continue"
Push-Location $OUT
& $R CMD build $Pkg 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
$tgz = Get-ChildItem -Path $OUT -Filter "*.tar.gz" | Select-Object -Last 1
if ($tgz) {
  & $R CMD check --as-cran $tgz.FullName 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
} else {
  Add-Content -Path $LOG -Value "NO TARBALL BUILT"
}
Pop-Location
Add-Content -Path $LOG -Value "CHECK DONE"
