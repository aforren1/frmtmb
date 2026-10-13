#!/usr/bin/env bash
# Lane seflags: roxygenise core and install it into the lane's private
# library (never anywhere else). Extensions are installed by name:
#   bash dev/seflags-install.sh [ext ...]
set -eu
W=C:/Users/adf44/source/r/frmtmb-wt-seflags
LIB=C:/Users/adf44/source/r/wt-seflags-lib
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
cd "/c/Users/adf44/source/r/frmtmb-wt-seflags"
"$RS" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r7', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$W')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$W"
for e in "$@"; do
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$W/extensions/$e"
done
