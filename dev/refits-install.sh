#!/usr/bin/env bash
# Roxygenise core and install it (and the named extensions) into the lane
# library. Usage: bash dev/refits-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-refits
LIB=C:/Users/adf44/source/r/wt-refits-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$WT"
"$RS" -e "roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
