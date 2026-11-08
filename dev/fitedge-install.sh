#!/usr/bin/env bash
# Roxygenise and install core (and any extension named) from the worktree
# into the lane's private library only.
#   bash dev/fitedge-install.sh [extension ...]
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-fitedge
LIB=C:/Users/adf44/source/r/wt-fitedge-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$RS" -e "roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e"
done
