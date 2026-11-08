#!/usr/bin/env bash
# Fix round: the changed and nearby test files, one process each, in sequence.
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
O=$WT/dev/fitedge-out/r2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-fitedge-lib
for f in "$@"; do
  "$RS" "$WT/dev/fitedge-run-tests.R" frmtmb "$WT/tests/testthat/$f.R" > "$O/t-$f.txt" 2>&1
  grep -h "^lib\|^RESULT\|^DETAIL" "$O/t-$f.txt" | cut -c1-240
done
