#!/bin/bash
# Roxygenise core and install it into the lane library.
#   bash dev/nanse-install.sh [pkgdir]   (default: core)
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
W=/c/Users/adf44/source/r/frmtmb-wt-nanse
P=${1:-.}
cd $W
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e ".libPaths(c('C:/Users/adf44/source/r/wt-nanse-lib','C:/Users/adf44/source/r/rellib-r5','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$P')" 2>&1 | tail -5
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/wt-nanse-lib --no-multiarch $P 2>&1 | tail -3
