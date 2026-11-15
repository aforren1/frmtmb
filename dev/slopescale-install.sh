#!/usr/bin/env bash
# Lane slopescale: roxygenise and install core (and the named extensions)
# into the lane's private library.
# Usage: bash dev/slopescale-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-slopescale
LIB=C:/Users/adf44/source/r/wt-slopescale-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
