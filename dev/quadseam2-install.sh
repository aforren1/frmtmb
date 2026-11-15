#!/usr/bin/env bash
# Roxygenise and install core (and any extensions named) from this
# worktree into wt-quadseam2-lib. Usage:
#   bash dev/quadseam2-install.sh [extension ...]
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
LIB=C:/Users/adf44/source/r/wt-quadseam2-lib
WTW=C:/Users/adf44/source/r/frmtmb-wt-quadseam2
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
mkdir -p "$WT/dev/quadseam2-log"
cd "$WT"
"$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r21','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WTW')" \
  > dev/quadseam2-log/install-rox.txt 2>&1 || { echo "roxygenise failed"; tail -20 dev/quadseam2-log/install-rox.txt; exit 1; }
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT" > dev/quadseam2-log/install-core.txt 2>&1 \
  || { echo "core install failed"; tail -20 dev/quadseam2-log/install-core.txt; exit 1; }
echo "core installed"
for ext in "$@"; do
  "$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r21','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$WTW/extensions/$ext')" \
    > "dev/quadseam2-log/install-rox-$ext.txt" 2>&1 || { echo "rox $ext failed"; exit 1; }
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$ext" \
    > "dev/quadseam2-log/install-$ext.txt" 2>&1 \
    || { echo "$ext install failed"; tail -20 "dev/quadseam2-log/install-$ext.txt"; exit 1; }
  echo "$ext installed"
done
