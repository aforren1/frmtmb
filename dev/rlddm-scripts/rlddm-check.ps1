# R CMD check --as-cran on one package, ONCE, on the final pass.
#
#   powershell -File dev\rlddm-scripts\rlddm-check.ps1 <pkg>
#
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR, and _R_CHECK_CRAN_INCOMING_REMOTE_ is off because the
# not-on-CRAN WARNING is pre-existing. No --no-build-vignettes: that
# flag manufactured two WARNINGs and a NOTE on all eight packages in the
# 0.55.2 round, including four the round never touched.

param([string]$Pkg)

$root = "C:\Users\adf44\source\r\frmtmb-wt-rlddm"
$lib = "C:\Users\adf44\source\r\rlddm-lib"
$work = "C:\Users\adf44\source\r\rlddm-check"
$env:PATH = "C:\rtools45\usr\bin;" +
  "C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" +
  "C:\Program Files\R\R-4.6.1\bin;" +
  "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;" +
  "C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH
$env:RSTUDIO_PANDOC =
  "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:R_LIBS = "$lib;C:\Users\adf44\source\r\pinlib;" +
  "C:\Users\adf44\source\r\rellib-0552;" +
  "C:\Users\adf44\AppData\Local\R\win-library\4.6"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
$env:NOT_CRAN = "false"

New-Item -ItemType Directory -Force $work | Out-Null
Set-Location $work
Get-ChildItem "$work\$Pkg`_*.tar.gz" -ErrorAction SilentlyContinue |
  Remove-Item -Force
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build "$root\extensions\$Pkg"
$tar = (Get-ChildItem "$work\$Pkg`_*.tar.gz" | Select-Object -First 1).Name
Write-Output "built $tar"
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran $tar
Write-Output "--- Status ---"
Get-Content "$work\$Pkg.Rcheck\00check.log" |
  Select-String -Pattern "^Status|WARNING|NOTE|ERROR"
