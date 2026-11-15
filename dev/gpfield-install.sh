#!/usr/bin/env bash
# Lane gpfield: roxygenise core and install it into wt-gpfield-lib.
# Usage: gpfield-install.sh <log> [extension ...]
WT=C:/Users/adf44/source/r/frmtmb-wt-gpfield
LIB=C:/Users/adf44/source/r/wt-gpfield-lib
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
log=$WT/dev/gpfield-out/$1; shift
mkdir -p $WT/dev/gpfield-out
{
Rscript -e "roxygen2::roxygenise('$WT')"
R CMD INSTALL --library=$LIB --no-multiarch $WT
for e in "$@"; do
  Rscript -e "roxygen2::roxygenise('$WT/extensions/$e')"
  R CMD INSTALL --library=$LIB --no-multiarch $WT/extensions/$e
done
} > $log 2>&1
echo "EXIT $? $log"; tail -3 $log
