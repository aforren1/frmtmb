# Lane silent2: roxygenise and install core and the named extensions into
# the lane library. Usage: powershell -File dev/silent2-build.ps1 [ext ...]
param([string[]]$exts = @())
$ErrorActionPreference = "Continue"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-silent2"
$lib = "C:/Users/adf44/source/r/wt-silent2-lib"
$env:PATH = $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$env:R_LIBS = "$lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $rs -e "roxygen2::roxygenise('$wt')"
& $r CMD INSTALL "--library=$lib" --no-multiarch $wt
foreach ($e in $exts) {
  & $rs -e "roxygen2::roxygenise('$wt/extensions/$e')"
  & $r CMD INSTALL "--library=$lib" --no-multiarch "$wt/extensions/$e"
}
