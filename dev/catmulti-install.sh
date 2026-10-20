#!/bin/sh
# Roxygenise and install lane catmulti's core (and, with an argument,
# frmtmb.sample) into the lane's private library.
#   sh dev/catmulti-install.sh [sample]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-catmulti
LIB=C:/Users/adf44/source/r/wt-catmulti-lib
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RCMD="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$WT"
if [ "$1" != "sample" ]; then
  "$R" -e "roxygen2::roxygenise('$WT')" > dev/catmulti-out/rox-core.txt 2>&1
  "$RCMD" CMD INSTALL --library="$LIB" --no-multiarch "$WT" \
    > dev/catmulti-out/install-core.txt 2>&1
  tail -3 dev/catmulti-out/install-core.txt
fi
if [ "$1" = "sample" ] || [ "$2" = "sample" ]; then
  "$R" -e "roxygen2::roxygenise('$WT/extensions/frmtmb.sample')" \
    > dev/catmulti-out/rox-sample.txt 2>&1
  "$RCMD" CMD INSTALL --library="$LIB" --no-multiarch \
    "$WT/extensions/frmtmb.sample" > dev/catmulti-out/install-sample.txt 2>&1
  tail -3 dev/catmulti-out/install-sample.txt
fi
