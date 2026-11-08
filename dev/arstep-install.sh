#!/usr/bin/env bash
# Lane arstep: roxygenise core and install it into the lane library.
export TMP="C:/Users/adf44/AppData/Local/Temp" TEMP="C:/Users/adf44/AppData/Local/Temp"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin"
cd /c/Users/adf44/source/r/frmtmb-wt-arstep
"$R/Rscript.exe" -e '.libPaths(c("C:/Users/adf44/source/r/wt-arstep-lib", "C:/Users/adf44/source/r/rellib-r18", "C:/Users/adf44/AppData/Local/R/win-library/4.6")); roxygen2::roxygenise(".")'
R_LIBS="C:/Users/adf44/source/r/wt-arstep-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6" "$R/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/wt-arstep-lib --no-multiarch ${1:-.}
