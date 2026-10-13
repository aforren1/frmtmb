#!/usr/bin/env bash
# Lane optrap, consolidation fixes: roxygenise and install the scratch
# merge (dev/optrap-relmerge/: the release tree with this lane's
# consolidation diffs applied) into optrap-relmerge-lib, ahead of the
# release's rellib-r8.
#   bash dev/optrap-c3-install-merge.sh
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-optrap/dev/optrap-relmerge
LIB=C:/Users/adf44/source/r/optrap-relmerge-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
mkdir -p "$LIB"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r8;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$RS" -e "roxygen2::roxygenise('$WT')"
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT"
