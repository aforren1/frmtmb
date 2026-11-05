#!/bin/bash
# Roxygenise and install core and the named extensions into the lane library.
#   bash dev/predread-install.sh [core] [ext ...]
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
WT=C:/Users/adf44/source/r/frmtmb-wt-predread
LIB=C:/Users/adf44/source/r/wt-predread-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in "$@"; do
  if [ "$p" = core ]; then d="$WT"; else d="$WT/extensions/$p"; fi
  "$RS" -e "roxygen2::roxygenise('$d')" 2>&1 | grep -v "^Writing\|^Loading\|^ℹ\|^i Loading" | tail -5
  "$R" CMD INSTALL --library="$LIB" --no-multiarch --no-docs "$d" 2>&1 | tail -2
done
