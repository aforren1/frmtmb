#!/bin/bash
# Roxygenise and install core (unless "noroot") and the named extensions
# into the lane library. Usage: bash dev/reslot-install.sh [noroot] [ext ...]
WT=C:/Users/adf44/source/r/frmtmb-wt-reslot
LIB=C:/Users/adf44/source/r/wt-reslot-lib
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r14;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
one() {
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$1')" || exit 1
  "$R" CMD INSTALL --library=$LIB --no-multiarch "$1" || exit 1
}
if [ "$1" = noroot ]; then shift; else one $WT; fi
for e in "$@"; do one $WT/extensions/$e; done
