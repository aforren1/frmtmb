#!/bin/sh
# Lane wt-defects: roxygenise and install core (and, with "sample",
# frmtmb.sample) into the lane's private library only.
#   sh dev/defects-install.sh [core] [sample]
cd "$(dirname "$0")/.." || exit 1
W="$(pwd -W)"
LIB=C:/Users/adf44/source/r/wt-defects-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="$R_LIBS"
[ $# -eq 0 ] && set -- core
for p in "$@"; do
  case "$p" in
    core) d="$W" ;;
    sample) d="$W/extensions/frmtmb.sample" ;;
    *) echo "unknown $p"; exit 1 ;;
  esac
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$d')" || exit 1
  "$RC" CMD INSTALL --library="$LIB" --no-multiarch "$d" || exit 1
done
