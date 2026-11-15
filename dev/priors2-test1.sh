#!/bin/bash
# One test file, one R process, through dev/release/run-tests.R.
# bash dev/priors2-test1.sh <package> <test file> [arm]
# arm "after" (default) puts the lane library first; "before" runs
# rellib-r21 alone. The log goes to dev/priors2-out/t1-<arm>-<file>.txt
# and its RESULT line is printed.
set -u
p="$1"; f="$2"; arm="${3:-after}"
WT=C:/Users/adf44/source/r/frmtmb-wt-priors2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
export NOT_CRAN=true
if [ "$arm" = "after" ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-priors2-lib
else
  unset FRMTMB_RELLIB
fi
mkdir -p "$WT/dev/priors2-out"
o="$WT/dev/priors2-out/t1-$arm-$(basename "$f" .R).txt"
"$RS" "$WT/dev/release/run-tests.R" "$p" "$f" > "$o" 2>&1
grep -a "^lib:" "$o"
grep -a "^RESULT" "$o" | tail -1
grep -a "^FAIL\|^ERROR\|^WARN" "$o" | head -20
