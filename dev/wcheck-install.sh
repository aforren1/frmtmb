#!/usr/bin/env bash
# Lane wcheck: roxygenise and install core (and any extension named)
# into the lane library. Rtools goes at the END of PATH (lane-env).
#   bash dev/wcheck-install.sh [extension ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-wcheck
LIB=C:/Users/adf44/source/r/wt-wcheck-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch "$WT"
for p in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$p')"
  "$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch "$WT/extensions/$p"
done
