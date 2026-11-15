#!/bin/bash
# Roxygenise core (and any extension named), install into the lane lib.
# bash dev/priors2-install.sh [ext ...]   e.g. frmtmb.sample frmtmb.spatial
# SKIP_CORE=1 installs only the extensions named.
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-priors2
LIB=C:/Users/adf44/source/r/wt-priors2-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p "$WT/dev/priors2-out"
if [ "${SKIP_CORE:-0}" != "1" ]; then
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')" > "$WT/dev/priors2-out/rox-core.log" 2>&1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT" > "$WT/dev/priors2-out/inst-core.log" 2>&1
  tail -2 "$WT/dev/priors2-out/inst-core.log"
fi
for e in "$@"; do
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/$e')" > "$WT/dev/priors2-out/rox-$e.log" 2>&1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/$e" > "$WT/dev/priors2-out/inst-$e.log" 2>&1
  tail -2 "$WT/dev/priors2-out/inst-$e.log"
done
