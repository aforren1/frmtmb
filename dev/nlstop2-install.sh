#!/bin/bash
# Lane nlstop2: roxygenise the worktree's core and install it into the
# lane's library. Usage: bash dev/nlstop2-install.sh <log name>
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
LIB=C:/Users/adf44/source/r/wt-nlstop2-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
log="dev/nlstop2-out/${1:-install.log}"
{
  R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
    "$RS" -e 'roxygen2::roxygenise("C:/Users/adf44/source/r/frmtmb-wt-nlstop2")'
  echo "ROXYGEN EXIT $?"
  R_LIBS="C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
    "$R" CMD INSTALL --library="$LIB" --no-multiarch .
  echo "INSTALL EXIT $?"
} > "$log" 2>&1
tail -3 "$log"
