#!/usr/bin/env bash
# Roxygenise and install core into the lane's private library.
#   bash dev/nlfix-install.sh [--no-roxy]
# Rtools goes at the END of PATH (see dev/equiv-run.sh).
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin"
LIB=C:/Users/adf44/source/r/wt-nlfix-lib2
mkdir -p "$LIB"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
if [ "${1:-}" != "--no-roxy" ]; then
  "$R/Rscript.exe" -e 'roxygen2::roxygenise("C:/Users/adf44/source/r/frmtmb-wt-nlfix")' > dev/nlfix-roxy.log 2>&1 || { tail -20 dev/nlfix-roxy.log; exit 1; }
fi
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch --no-test-load . > dev/nlfix-install.log 2>&1 || { tail -30 dev/nlfix-install.log; exit 1; }
tail -2 dev/nlfix-install.log
