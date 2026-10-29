#!/bin/bash
# Lane meme: roxygenise and install core, then the named extensions,
# into the lane library. Usage: bash dev/meme-install.sh [ext ...]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-meme
LIB=C:/Users/adf44/source/r/wt-meme-lib
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RCMD="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r15;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"$R" -e "roxygen2::roxygenise('$WT')"
"$RCMD" CMD INSTALL --library=$LIB --no-multiarch "$WT"
for e in "$@"; do
  "$R" -e "roxygen2::roxygenise('$WT/extensions/$e')"
  "$RCMD" CMD INSTALL --library=$LIB --no-multiarch "$WT/extensions/$e"
done
