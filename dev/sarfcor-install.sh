#!/bin/sh
# Roxygenise and install core (and the named extensions) into the lane
# library. Usage: sh dev/sarfcor-install.sh [core] [ext ...]
WT=C:/Users/adf44/source/r/frmtmb-wt-sarfcor
LIB=C:/Users/adf44/source/r/wt-sarfcor-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
for p in "$@"; do
  if [ "$p" = core ]; then dir=$WT; else dir=$WT/extensions/$p; fi
  "$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r12','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$dir')" || exit 1
  "$R" CMD INSTALL --library=$LIB --no-multiarch "$dir" || exit 1
done
