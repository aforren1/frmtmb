#!/bin/sh
# Reviewer: roxygenise and install the trial merge (release tree + lane
# nanse) into the reviewer's own library. Usage: sh dev/nanse-rev-install.sh
T=C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/nanse-rev-merge-check/tree
LIB=C:/Users/adf44/source/r/nanse-rev-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r6','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$T')" || exit 1
"$R" CMD INSTALL --library=$LIB --no-multiarch $T || exit 1
