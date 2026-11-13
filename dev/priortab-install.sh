#!/bin/bash
# Roxygenise core (and any extension named), install into the lane lib.
# bash dev/priortab-install.sh [ext ...]   e.g. frmtmb.sample frmtmb.spatial
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-priortab
LIB=C:/Users/adf44/source/r/wt-priortab-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p "$WT/dev/priortab-out"
if [ "${SKIP_CORE:-0}" != "1" ]; then
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')" > "$WT/dev/priortab-out/rox-core.log" 2>&1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT" > "$WT/dev/priortab-out/inst-core.log" 2>&1
  tail -2 "$WT/dev/priortab-out/inst-core.log"
fi
for e in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$e')" > "$WT/dev/priortab-out/rox-$e.log" 2>&1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e" > "$WT/dev/priortab-out/inst-$e.log" 2>&1
  tail -2 "$WT/dev/priortab-out/inst-$e.log"
done
