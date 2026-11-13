#!/usr/bin/env bash
# Lane names: roxygenise and install core and the named extensions from
# this worktree into the lane library. Usage:
#   bash dev/names-install.sh [core] [ext ...]
cd /c/Users/adf44/source/r/frmtmb-wt-names || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
LIB="C:/Users/adf44/source/r/wt-names-lib"
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in "$@"; do
  if [ "$p" = core ]; then dir=.; else dir="extensions/$p"; fi
  "$R/Rscript.exe" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$dir')" > "dev/names-out/rox-$p.log" 2>&1 || { echo "ROX FAIL $p"; tail -20 "dev/names-out/rox-$p.log"; exit 1; }
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$dir" > "dev/names-out/inst-$p.log" 2>&1 || { echo "INSTALL FAIL $p"; tail -20 "dev/names-out/inst-$p.log"; exit 1; }
  echo "installed $p"
done
