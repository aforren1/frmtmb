#!/usr/bin/env bash
# Lane optstart: roxygenise core and install it into the lane library.
#   bash dev/optstart-install.sh <tag>
set -u
ROOT=C:/Users/adf44/source/r/frmtmb-wt-optstart
LIB=C:/Users/adf44/source/r/wt-optstart-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
log="$ROOT/dev/optstart-log/install-$1.txt"
"$RS" -e ".libPaths(c('$LIB', 'C:/Users/adf44/source/r/rellib-r8', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('$ROOT')" > "$log" 2>&1
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$ROOT" >> "$log" 2>&1
echo "exit $?" >> "$log"
tail -n 3 "$log"
