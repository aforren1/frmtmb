#!/usr/bin/env bash
# Lane arspat: roxygenise and install core (and the named extensions)
# into the lane's private library. Usage:
#   bash dev/arspat-install.sh [core] [frmtmb.spatial] [frmtmb.sample] ...
set -euo pipefail
WT=C:/Users/adf44/source/r/frmtmb-wt-arspat
LIB=C:/Users/adf44/source/r/wt-arspat-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in "$@"; do
  if [ "$p" = core ]; then dir="$WT"; else dir="$WT/extensions/$p"; fi
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$dir')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$dir"
done
