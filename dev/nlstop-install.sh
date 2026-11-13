#!/usr/bin/env bash
# Lane nlstop: roxygenise core and install it into the lane's library.
# Usage: bash dev/nlstop-install.sh [extension ...]
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
LIB=C:/Users/adf44/source/r/wt-nlstop-lib
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$RS" -e 'roxygen2::roxygenise(".")' || exit 1
"$R" CMD INSTALL --library="$LIB" --no-multiarch . || exit 1
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('extensions/$e')" || exit 1
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "extensions/$e" || exit 1
done
echo INSTALL DONE
