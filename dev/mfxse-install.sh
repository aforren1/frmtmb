#!/usr/bin/env bash
# Roxygenise core, then install core (and the extensions named as
# arguments) from the lane's worktree into the lane's library.
#   bash dev/mfxse-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-mfxse
LIB=C:/Users/adf44/source/r/wt-mfxse-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
