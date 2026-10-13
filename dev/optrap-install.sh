#!/usr/bin/env bash
# Lane optrap: roxygenise and install core (and any extensions named) into
# the lane's private library.
#   bash dev/optrap-install.sh [ext ...]
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-optrap
LIB=C:/Users/adf44/source/r/wt-optrap-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r7;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd "$WT"
"$RS" -e "roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
