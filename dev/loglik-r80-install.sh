#!/usr/bin/env bash
# Round r080, lane loglik: roxygenise and install core (unless
# SKIP_CORE=1) and the named extensions into the lane library, against
# rellib-r18. Usage: bash dev/loglik-r80-install.sh [ext ...]
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-loglik
LIB=C:/Users/adf44/source/r/wt-loglik-lib
RB="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
cd "$WT"
LP=".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]])"
if [ "${SKIP_CORE:-0}" != 1 ]; then
  "$RB/Rscript.exe" -e "$LP; roxygen2::roxygenise('$WT')"
  "$RB/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
fi
for e in "$@"; do
  "$RB/Rscript.exe" -e "$LP; roxygen2::roxygenise('$WT/extensions/$e')"
  "$RB/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
