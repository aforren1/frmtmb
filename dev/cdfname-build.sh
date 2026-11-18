#!/usr/bin/env bash
# Roxygenise and install core and the named extensions into the lane's
# library.  bash dev/cdfname-build.sh [core] [ext ...]
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-cdfname
LIB=C:/Users/adf44/source/r/wt-cdfname-lib
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
for p in "$@"; do
  if [ "$p" = core ]; then d=$WT; else d=$WT/extensions/$p; fi
  "$RS" -e ".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r22','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$d')" || exit 1
  "$R" CMD INSTALL --library=$LIB --no-multiarch "$d" || exit 1
  echo "BUILT $p"
done
