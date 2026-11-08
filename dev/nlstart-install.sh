#!/usr/bin/env bash
# Roxygenise core and install it into the lane's library.
# Usage: bash dev/nlstart-install.sh [extension ...]
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
LIB=C:/Users/adf44/source/r/wt-nlstart-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$RS" -e 'roxygen2::roxygenise("C:/Users/adf44/source/r/frmtmb-wt-nlstart")' || exit 1
"$R" CMD INSTALL --library="$LIB" --no-multiarch . || exit 1
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-nlstart/extensions/$e')" || exit 1
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "extensions/$e" || exit 1
done
echo INSTALL OK
