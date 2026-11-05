#!/usr/bin/env bash
# Roxygenise and install core, and optionally named extensions, into the lane
# library. Usage: bash dev/ncgp-install.sh [core] [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-ncgp
LIB=C:/Users/adf44/source/r/wt-ncgp-lib
RBIN="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
mkdir -p "$LIB"
for p in "$@"; do
  if [ "$p" = core ]; then dir="$WT"; else dir="$WT/extensions/$p"; fi
  "$RBIN/Rscript.exe" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r16','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$dir')" > "$WT/dev/ncgp-out/rox-$p.txt" 2>&1 || { tail -20 "$WT/dev/ncgp-out/rox-$p.txt"; exit 1; }
  R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6" "$RBIN/R.exe" CMD INSTALL --library="$LIB" --no-multiarch --no-docs --no-test-load ${NCGP_FAST:+--no-byte-compile} "$dir" > "$WT/dev/ncgp-out/inst-$p.txt" 2>&1 || { tail -20 "$WT/dev/ncgp-out/inst-$p.txt"; exit 1; }
  tail -1 "$WT/dev/ncgp-out/inst-$p.txt"
done
