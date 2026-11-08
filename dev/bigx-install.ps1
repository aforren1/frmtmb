# Roxygenise and install lane bigx's core (and, with arguments, the named
# extensions) into the lane's private library.
#   powershell -File dev/bigx-install.ps1 [ext ...]
$ErrorActionPreference = "Continue"
$env:PATH = $env:PATH + ";C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
$wt = "C:/Users/adf44/source/r/frmtmb-wt-bigx"
$lib = "C:/Users/adf44/source/r/wt-bigx-lib"
$env:R_LIBS = "$lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$r = "C:\Program Files\R\R-4.6.1\bin\R.exe"
& $rs -e "roxygen2::roxygenise('$wt')"
& $r CMD INSTALL --library=$lib --no-multiarch $wt
foreach ($e in $args) {
  & $rs -e "roxygen2::roxygenise('$wt/extensions/$e')"
  & $r CMD INSTALL --library=$lib --no-multiarch "$wt/extensions/$e"
}
