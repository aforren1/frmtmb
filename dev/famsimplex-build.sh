#!/bin/sh
# Roxygenise core and install it (and any named extensions) into the lane
# library. Usage: sh dev/famsimplex-build.sh [frmtmb.sample ...]
WT=/c/Users/adf44/source/r/frmtmb-wt-famsimplex
LIB=C:/Users/adf44/source/r/wt-famsimplex-lib
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RCMD="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$R" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r8','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-famsimplex')" || exit 1
"$RCMD" CMD INSTALL --library=$LIB --no-multiarch "$WT" || exit 1
for p in "$@"; do
  "$R" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r8','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-famsimplex/extensions/$p')" || exit 1
  "$RCMD" CMD INSTALL --library=$LIB --no-multiarch "$WT/extensions/$p" || exit 1
done
