#!/usr/bin/env bash
# Roxygenise and install core (and any extensions named) from this
# worktree into the lane library.
#   bash dev/hygiene-install.sh [core] [frmtmb.ode ...]
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd -W)"
LIB=C:/Users/adf44/source/r/wt-hygiene-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for p in "$@"; do
  if [ "$p" = core ]; then dir="$ROOT"; else dir="$ROOT/extensions/$p"; fi
  echo "== roxygenise $dir"
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$dir')"
  echo "== install $dir"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$dir"
done
