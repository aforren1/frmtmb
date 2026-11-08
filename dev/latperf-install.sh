#!/usr/bin/env bash
# Lane latperf: roxygenise core and install it into the lane library.
# Usage: bash dev/latperf-install.sh [extension ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-latperf
LIB=C:/Users/adf44/source/r/wt-latperf-lib
RB="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
"$RB/Rscript.exe" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r18', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT')"
"$RB/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for ext in "$@"; do
  "$RB/Rscript.exe" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r18', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT/extensions/$ext')"
  R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
    "$RB/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$ext"
done
