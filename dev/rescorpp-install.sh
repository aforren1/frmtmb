#!/bin/sh
# Lane rescorpp: roxygenise and install core and frmtmb.sample into the
# lane library. Usage: sh dev/rescorpp-install.sh [core] [sample]
WT=C:/Users/adf44/source/r/frmtmb-wt-rescorpp
LIB=C:/Users/adf44/source/r/wt-rescorpp-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r15;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p "$LIB"
for p in "$@"; do
  case $p in
    core) dir=$WT ;;
    sample) dir=$WT/extensions/frmtmb.sample ;;
    *) dir=$WT/extensions/$p ;;
  esac
  echo "== roxygenise $dir"
  "$R/Rscript.exe" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$dir')" || exit 1
  echo "== install $dir"
  "$R/R.exe" CMD INSTALL --library=$LIB --no-multiarch "$dir" || exit 1
done
