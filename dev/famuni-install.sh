#!/bin/sh
# Roxygenise and install core (and the named extensions) into the lane's
# private library. Usage: sh dev/famuni-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-famuni
LIB=C:/Users/adf44/source/r/wt-famuni-lib
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RCMD="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r8;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R" -e "roxygen2::roxygenise('$WT')"
"$RCMD" CMD INSTALL --library=$LIB --no-multiarch $WT
for e in "$@"; do
  "$R" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$RCMD" CMD INSTALL --library=$LIB --no-multiarch $WT/extensions/$e
done
