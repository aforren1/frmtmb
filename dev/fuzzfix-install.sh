#!/usr/bin/env bash
# roxygenise and install core (and any extension named) into the lane lib
set -e
WT=/c/Users/adf44/source/r/frmtmb-wt-fuzzfix
LIB=C:/Users/adf44/source/r/wt-fuzzfix-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R/Rscript.exe" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-fuzzfix')"
"$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch "C:/Users/adf44/source/r/frmtmb-wt-fuzzfix"
for p in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-fuzzfix/extensions/$p')"
  "$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch "C:/Users/adf44/source/r/frmtmb-wt-fuzzfix/extensions/$p"
done
