#!/usr/bin/env bash
# Roxygenise and install core (and, with an argument, one extension)
# from the cdfs worktree into the lane's private library.
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-cdfs
LIB=C:/Users/adf44/source/r/wt-cdfs-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6"
pkg=${1:-core}
if [ "$pkg" = core ]; then
  "$RS" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r10', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
else
  "$RS" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r10', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT/extensions/$pkg')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$pkg"
fi
