#!/bin/sh
# Reviewer (final check after punch round 2): roxygenise and install the
# final trial merge (current release tree + lane nanse) into the
# reviewer's own library: core, then frmtmb.spline (R code merged),
# frmtmb.learn (floor) and frmtmb.sample (comment).
#   sh dev/nanse-rev3-install.sh
T=C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/nanse-rev3-merge-check/tree
LIB=C:/Users/adf44/source/r/nanse-rev-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/nanse-rev-lib;C:/Users/adf44/source/r/rellib-r6;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
LP=".libPaths(c('$LIB','C:/Users/adf44/source/r/rellib-r6','C:/Users/adf44/AppData/Local/R/win-library/4.6'))"
"$RS" -e "$LP; roxygen2::roxygenise('$T')" || exit 1
"$R" CMD INSTALL --library=$LIB --no-multiarch $T || exit 1
for p in frmtmb.spline frmtmb.learn frmtmb.sample; do
  "$RS" -e "$LP; roxygen2::roxygenise('$T/extensions/$p')" || exit 1
  "$R" CMD INSTALL --library=$LIB --no-multiarch $T/extensions/$p || exit 1
done
