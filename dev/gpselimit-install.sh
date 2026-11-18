#!/bin/bash
# Lane gpselimit: roxygenise core and install it (and, with "ext", every
# extension) from this worktree into the lane's private library.
# Usage: bash dev/gpselimit-install.sh [ext]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-gpselimit
LIB=C:/Users/adf44/source/r/wt-gpselimit-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$WT"
"$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')"
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
if [ "$1" = "ext" ]; then
  for p in coupling eam latent learn ode sample spatial survival; do
    "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/frmtmb.$p')"
    "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/frmtmb.$p"
  done
fi
