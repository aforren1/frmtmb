#!/usr/bin/env bash
# Lane quadesc: roxygenise and install core (and the named extensions)
# into the lane's private library.
#   bash dev/quadesc-install.sh [extension...]
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-quadesc
LIB=C:/Users/adf44/source/r/wt-quadesc-lib
cd "$WT"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$WT/extensions/$e')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
