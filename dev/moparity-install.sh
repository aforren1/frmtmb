#!/bin/sh
# Roxygenise and install core frmtmb (and frmtmb.sample when asked) into
# the lane's private library. Usage: sh dev/moparity-install.sh [sample|only-sample]
WT=C:/Users/adf44/source/r/frmtmb-wt-moparity
LIB=C:/Users/adf44/source/r/wt-moparity-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r7;C:/Users/adf44/AppData/Local/R/win-library/4.6"
if [ "$1" != "only-sample" ]; then
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')" || exit 1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT" || exit 1
fi
if [ "$1" = "sample" ] || [ "$1" = "only-sample" ]; then
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/frmtmb.sample')" || exit 1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/frmtmb.sample" || exit 1
fi
