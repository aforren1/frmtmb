#!/usr/bin/env bash
# Reviewer: install the A-only tree (core, then frmtmb.spatial, then
# frmtmb.sample) into the review library.
set -eu
WT=/c/Users/adf44/source/r/frmtmb-wt-sarfcor
SRC=$WT/dev/review-sarfcor-A
LIB=C:/Users/adf44/source/r/rev-sarfcor-lib
mkdir -p /c/Users/adf44/source/r/rev-sarfcor-lib
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
cd "$SRC"
"$R" CMD INSTALL --library=$LIB --no-multiarch --no-test-load .
cd "$SRC/extensions"
"$R" CMD INSTALL --library=$LIB --no-multiarch frmtmb.spatial
"$R" CMD INSTALL --library=$LIB --no-multiarch frmtmb.sample
echo INSTALL-DONE
