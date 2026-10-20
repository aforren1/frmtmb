#!/bin/bash
# Lane illse2: roxygenise core and install it (and, with args, the named
# extensions) into the lane library. Log: dev/illse2-log/install-<n>.txt
set -e
WT=/c/Users/adf44/source/r/frmtmb-wt-illse2
LIB=C:/Users/adf44/source/r/wt-illse2-lib
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$RS" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-illse2')"
"$R" CMD INSTALL --library=$LIB --no-multiarch "$WT"
for p in "$@"; do
  "$RS" -e "roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-illse2/extensions/$p')"
  "$R" CMD INSTALL --library=$LIB --no-multiarch "$WT/extensions/$p"
done
