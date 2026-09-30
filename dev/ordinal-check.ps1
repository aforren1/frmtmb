# Build and check one package with --as-cran inside dev/ordinal-check.
# Built WITH vignettes and the manual (dev/release/run-check.ps1 says
# why). Usage: powershell -File dev/ordinal-check.ps1 <pkgdir> <tag>
param([string]$PkgDir, [string]$Tag)
$WT = "C:\Users\adf44\source\r\frmtmb-wt-ordinal"
$CK = "$WT\dev\ordinal-check\$Tag"
if (Test-Path $CK) { Remove-Item -Recurse -Force $CK }
New-Item -ItemType Directory -Force $CK | Out-Null
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:TMP = "C:\Users\adf44\AppData\Local\Temp\1"
$env:TEMP = $env:TMP
$env:R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:R_LIBS = "C:/Users/adf44/source/r/wt-ordinal-lib;C:/Users/adf44/source/r/rellib-r4;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$env:NOT_CRAN = ""
$ErrorActionPreference = "Continue"
Set-Location $CK
$R = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $R CMD build $PkgDir 2>&1 | ForEach-Object { "$_" } |
  Out-File -FilePath "$CK\build.log" -Encoding utf8
$tar = Get-ChildItem "$CK\*.tar.gz" | Sort-Object LastWriteTime |
  Select-Object -Last 1
& $R CMD check --as-cran $tar.FullName 2>&1 | ForEach-Object { "$_" } |
  Out-File -FilePath "$CK\check.log" -Encoding utf8
