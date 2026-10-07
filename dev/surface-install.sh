#!/usr/bin/env bash
# Lane surface: roxygenise and install core, then the named extensions,
# into the lane's private library only.
#
#   bash dev/surface-install.sh [extension ...]
set -eu
wt="$(cd "$(dirname "$0")/.." && pwd -W)"
lib="C:/Users/adf44/source/r/wt-surface-lib"
rs="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
r="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$lib;C:/Users/adf44/source/r/rellib-r6;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$rs" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$wt')"
"$r" CMD INSTALL --library="$lib" --no-multiarch "$wt"
for e in "$@"; do
  "$rs" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$wt/extensions/$e')"
  "$r" CMD INSTALL --library="$lib" --no-multiarch "$wt/extensions/$e"
done
echo "INSTALL DONE"
