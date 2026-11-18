#!/usr/bin/env bash
# Lane nlband: roxygenise core and install it (and, with arguments, the
# named extensions) into the lane's private library.
# Usage: bash dev/nlband-install.sh <log> [ext ...]
log="$1"; shift
WT=/c/Users/adf44/source/r/frmtmb-wt-nlband
LIB=C:/Users/adf44/source/r/wt-nlband-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
{
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-nlband')" &&
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT" &&
  for e in "$@"; do
    "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e" || exit 1
  done
  echo "INSTALL EXIT $?"
} > "$log" 2>&1
tail -3 "$log"
