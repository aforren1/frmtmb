#!/bin/sh
# Lane loglik: roxygenise and install core, then the named extensions,
# into the lane's private library.
#   sh dev/loglik-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-loglik
LIB=C:/Users/adf44/source/r/wt-loglik-lib
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r13;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd "$WT"
"$R/Rscript.exe" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$R/Rscript.exe" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$WT/extensions/$e')"
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
