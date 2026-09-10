# lane gddm: build the tarball and run R CMD check --as-cran once, on
# the final pass. pandoc and TinyTeX have to be on PATH or the check
# reports a bogus pdflatex ERROR.
$ErrorActionPreference = "Continue"
$root = "C:\Users\adf44\source\r\frmtmb-wt-gddm"
$out = "C:\Users\adf44\source\r\gddm-check"
$lib = "C:\Users\adf44\source\r\gddm-lib"
$ref = "C:\Users\adf44\source\r\rellib-0552"
$pin = "C:\Users\adf44\source\r\pinlib"
$usr = "C:\Users\adf44\AppData\Local\R\win-library\4.6"
if (-not (Test-Path $out)) { New-Item -ItemType Directory $out | Out-Null }
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;C:\Program Files\R\R-4.6.1\bin;" + $env:PATH
$env:R_LIBS = "$lib;$ref;$pin;$usr"
$env:NOT_CRAN = "true"
$env:_R_CHECK_CRAN_INCOMING_REMOTE_ = "FALSE"
Set-Location $out
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD build "$root\extensions\frmtmb.eam"
$tar = Get-ChildItem -Path $out -Filter "frmtmb.eam_*.tar.gz" |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
"tarball: " + $tar.FullName
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD check --as-cran --library=$lib $tar.FullName
"CHECK-DONE"
