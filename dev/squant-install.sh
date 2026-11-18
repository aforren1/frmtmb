#!/bin/sh
# Lane squant: roxygenise core and install it into the lane library.
# Run: sh dev/squant-install.sh  (from the worktree)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
WT=C:/Users/adf44/source/r/frmtmb-wt-squant
LIB=C:/Users/adf44/source/r/wt-squant-lib
R="/c/Program Files/R/R-4.6.1/bin"
"$R/Rscript.exe" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r22','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT')" || exit 1
"$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch $WT
