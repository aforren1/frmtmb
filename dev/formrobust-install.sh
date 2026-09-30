#!/usr/bin/env bash
# Roxygenise and install core (and any extension named as an argument)
# into the lane's private library. Usage:
#   bash dev/formrobust-install.sh [frmtmb.sample ...]
set -eu
WT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
LIB=C:/Users/adf44/source/r/wt-formrobust-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r4;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r4','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-formrobust')" 2>&1 | tail -5
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT" 2>&1 | tail -3
for p in "$@"; do
  "$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r4','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('C:/Users/adf44/source/r/frmtmb-wt-formrobust/extensions/$p')" 2>&1 | tail -5
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$p" 2>&1 | tail -3
done
