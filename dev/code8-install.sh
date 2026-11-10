#!/usr/bin/env bash
# Lane code8: roxygenise core and install it into the lane library.
set -e
cd "$(dirname "$0")/.."
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
LIB=C:/Users/adf44/source/r/wt-code8-lib
"$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r19','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('.')" 2>&1 | tail -5
"$R" CMD INSTALL --library=$LIB --no-multiarch . 2>&1 | tail -3
