# R CMD check --as-cran on one package (core or frmtmb.sample), built
# and checked inside dev/loodraws-check/<pkg>/.
#   powershell -File dev/loodraws-check.ps1 <frmtmb|frmtmb.sample>
# Log: dev/loodraws-check/<pkg>/check.log
param([string]$Pkg = "frmtmb.sample")
$ErrorActionPreference = "Stop"
$WT  = "C:/Users/adf44/source/r/frmtmb-wt-loodraws"
$R   = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$OUT = "$WT/dev/loodraws-check/$Pkg"
$SRC = if ($Pkg -eq "frmtmb") { $WT } else { "$WT/extensions/$Pkg" }
$LOG = "$OUT/check.log"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-loodraws-lib;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
$ErrorActionPreference = "Continue"
Push-Location $OUT
& $R CMD build "$SRC" 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
$tgz = Get-ChildItem -Path $OUT -Filter "*.tar.gz" | Select-Object -Last 1
if ($tgz) {
  & $R CMD check --as-cran $tgz.FullName 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
} else {
  Add-Content -Path $LOG -Value "NO TARBALL BUILT"
}
Pop-Location
