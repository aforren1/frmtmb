#!/bin/sh
# Roxygenise and install core (and any extension named) into the lane's
# private library.  sh dev/nb2dm-install.sh [ext ...]
set -e
cd "$(dirname "$0")/.."
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
LIB=C:/Users/adf44/source/r/wt-nb2dm-lib
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$R" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r10','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('.')"
"$RC" CMD INSTALL --library=$LIB --no-multiarch .
for e in "$@"; do
  "$R" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r10','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('extensions/$e')"
  "$RC" CMD INSTALL --library=$LIB --no-multiarch extensions/$e
done
