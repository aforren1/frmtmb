#!/usr/bin/env bash
# Lane gpstop: roxygenise core and install it into the lane's library.
# Args: extra package dirs (relative to the worktree) to install after core.
set -e
. "$(dirname "$0")/gpstop-env.sh"
WT=C:/Users/adf44/source/r/frmtmb-wt-gpstop
LIB=C:/Users/adf44/source/r/wt-gpstop-lib
"$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r20','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WT')"
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=$LIB --no-multiarch $WT
for p in "$@"; do
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=$LIB --no-multiarch $WT/$p
done
