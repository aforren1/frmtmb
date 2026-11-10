#!/usr/bin/env bash
# Lane support: roxygenise core and install it, then each extension
# named on the command line, into the lane's private library.
#
#   bash dev/support-install.sh [extension ...]
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-support
LIB=C:/Users/adf44/source/r/wt-support-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$RS" -e "roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
echo "INSTALL DONE"
