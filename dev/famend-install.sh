#!/usr/bin/env bash
# Lane famend: roxygenise and install core (and, with arguments, the
# named extensions) into the lane's private library.
#   bash dev/famend-install.sh [extension ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-famend
LIB=C:/Users/adf44/source/r/wt-famend-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for ext in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$ext')"
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$ext"
done
