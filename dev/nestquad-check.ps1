# Lane nestquad: R CMD check --as-cran on core frmtmb, built with
# vignettes, inside dev/nestquad-check/ (dev/release/run-check.ps1's
# settings, the lane's library first).
$ErrorActionPreference = "Stop"
$ROOT = "C:/Users/adf44/source/r/frmtmb-wt-nestquad"
$R    = "C:/Program Files/R/R-4.6.1/bin/R.exe"
$OUT  = "$ROOT/dev/nestquad-check"
$LOG  = "$OUT/check.log"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-nestquad-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = "true"
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
if (Test-Path $OUT) { Remove-Item -Recurse -Force $OUT }
New-Item -ItemType Directory -Force $OUT | Out-Null
$ErrorActionPreference = "Continue"
Push-Location $OUT
& $R CMD build $ROOT 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
$tgz = Get-ChildItem -Path $OUT -Filter "*.tar.gz" | Select-Object -Last 1
if ($tgz) {
  & $R CMD check --as-cran $tgz.FullName 2>&1 | Out-File -FilePath $LOG -Append -Encoding utf8
} else {
  Add-Content -Path $LOG -Value "NO TARBALL BUILT"
}
Pop-Location
Write-Output "CHECK DONE"
