#!/usr/bin/env bash
# Lane slopetrig: roxygenise core and install it (and any extensions
# named) into the lane library.
# Usage: bash dev/slopetrig-install.sh <log> [ext ...]
log="$1"; shift
WT=C:/Users/adf44/source/r/frmtmb-wt-slopetrig
LIB=C:/Users/adf44/source/r/wt-slopetrig-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
{
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')"
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
  for e in "$@"; do
    "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$e')"
    "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
  done
} > "$log" 2>&1
echo "EXIT $?" >> "$log"
