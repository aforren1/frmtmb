#!/usr/bin/env bash
# roxygenise core and install it into the lane library; log
# dev/quadseam-log/install.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam
LIB=C:/Users/adf44/source/r/wt-quadseam-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
mkdir -p $WT/dev/quadseam-log
R=/c/Program\ Files/R/R-4.6.1/bin
{
"$R/Rscript.exe" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r20','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-quadseam')"
"$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch C:/Users/adf44/source/r/frmtmb-wt-quadseam
} > $WT/dev/quadseam-log/install.txt 2>&1
echo "exit $?"; tail -3 $WT/dev/quadseam-log/install.txt
