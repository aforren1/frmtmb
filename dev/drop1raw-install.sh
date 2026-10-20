#!/bin/sh
# Roxygenise core and install it (and any extension named as an
# argument) into the lane library. Log: dev/drop1raw-out/install.txt
cd /c/Users/adf44/source/r/frmtmb-wt-drop1raw || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
R="/c/Program Files/R/R-4.6.1/bin"
LIB=C:/Users/adf44/source/r/wt-drop1raw-lib
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e "roxygen2::roxygenise('.')" || exit 1
"$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch . || exit 1
for p in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('extensions/$p')" || exit 1
  "$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch extensions/$p || exit 1
done
