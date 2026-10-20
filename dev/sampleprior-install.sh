#!/usr/bin/env bash
# Lane sampleprior: roxygenise and install core and frmtmb.sample (or
# the packages named) into the lane library.
#   bash dev/sampleprior-install.sh [core] [sample] [eam] ...
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-sampleprior
LIB=C:/Users/adf44/source/r/wt-sampleprior-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
export TEMP="$TMP"
pkgs="${*:-core sample}"
for p in $pkgs; do
  if [ "$p" = core ]; then dir="$WT"; else dir="$WT/extensions/frmtmb.$p"; fi
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$dir')" 2>&1 | tail -3
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$dir" 2>&1 | tail -2
done
