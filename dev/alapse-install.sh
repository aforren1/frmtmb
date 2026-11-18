#!/usr/bin/env bash
# Lane alapse: roxygenise and install core (and the named extensions)
# into the lane library.
#   bash dev/alapse-install.sh [ext ...]
set -u
ROOT=C:/Users/adf44/source/r/frmtmb-wt-alapse
LIB=C:/Users/adf44/source/r/wt-alapse-lib
cd "$ROOT"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p dev/alapse-out
"$RS" -e "roxygen2::roxygenise('$ROOT')" > dev/alapse-out/install-rox.txt 2>&1 || { echo ROX FAIL; tail dev/alapse-out/install-rox.txt; exit 1; }
"$R" CMD INSTALL --library="$LIB" --no-multiarch "$ROOT" > dev/alapse-out/install-core.txt 2>&1 || { echo CORE FAIL; tail dev/alapse-out/install-core.txt; exit 1; }
for e in "$@"; do
  "$RS" -e "roxygen2::roxygenise('$ROOT/extensions/$e')" > dev/alapse-out/install-rox-$e.txt 2>&1 || { echo ROX FAIL $e; exit 1; }
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$ROOT/extensions/$e" > dev/alapse-out/install-$e.txt 2>&1 || { echo FAIL $e; tail dev/alapse-out/install-$e.txt; exit 1; }
done
echo INSTALLED
