#!/bin/sh
# Roxygenise and install core frmtmb (and frmtmb.sample when asked) into
# the lane's private library. Usage: sh dev/aterms2-install.sh [sample]
WT=/c/Users/adf44/source/r/frmtmb-wt-aterms2
LIB=C:/Users/adf44/source/r/wt-aterms2-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$WT" || exit 1
"$R/Rscript.exe" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-aterms2')" || exit 1
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch C:/Users/adf44/source/r/frmtmb-wt-aterms2 || exit 1
if [ "$1" = "sample" ]; then
  "$R/Rscript.exe" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-aterms2/extensions/frmtmb.sample')" || exit 1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch C:/Users/adf44/source/r/frmtmb-wt-aterms2/extensions/frmtmb.sample || exit 1
fi
