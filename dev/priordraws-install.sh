#!/usr/bin/env bash
# Roxygenise and install core and frmtmb.sample from this worktree into
# the lane's private library. Usage: bash dev/priordraws-install.sh [core|sample|all]
set -euo pipefail
WT=/c/Users/adf44/source/r/frmtmb-wt-priordraws
LIB=C:/Users/adf44/source/r/wt-priordraws-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r14;C:/Users/adf44/AppData/Local/R/win-library/4.6"
what=${1:-all}
rox() {
  "$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r14','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$1')"
}
if [ "$what" = core ] || [ "$what" = all ]; then
  rox "C:/Users/adf44/source/r/frmtmb-wt-priordraws"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "C:/Users/adf44/source/r/frmtmb-wt-priordraws"
fi
if [ "$what" = sample ] || [ "$what" = all ]; then
  rox "C:/Users/adf44/source/r/frmtmb-wt-priordraws/extensions/frmtmb.sample"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "C:/Users/adf44/source/r/frmtmb-wt-priordraws/extensions/frmtmb.sample"
fi
