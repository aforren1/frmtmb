#!/usr/bin/env bash
# Roxygenise and install core (and the packages named as arguments,
# under extensions/) into the lane library.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd -W)"
LIB=C:/Users/adf44/source/r/wt-spatial-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in "$@"; do
  if [ "$p" = core ]; then d="$ROOT"; else d="$ROOT/extensions/$p"; fi
  "$RS" -e "roxygen2::roxygenise('$d')"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$d"
done
