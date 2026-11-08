#!/usr/bin/env bash
# Fix round: test-gp-start.R on the one-step build, its mutants, and the
# review's start probe, one process at a time.
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
O=$WT/dev/fitedge-out/r2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-fitedge-lib
"$RS" "$WT/dev/fitedge-run-tests.R" frmtmb "$WT/tests/testthat/test-gp-start.R" > "$O/t-gp-start-1step.txt" 2>&1
for m in gn2 noothers prefix; do
  "$RS" "$WT/dev/fitedge-mutrun.R" $m test-gp-start.R > "$O/mut1-$m-test-gp-start.R.txt" 2>&1
done
"$RS" "$WT/dev/review-fitedge-gnprobe.R" lane > "$O/gnprobe-1step.txt" 2>&1
grep -h "^lib\|^RESULT\|^DETAIL" "$O/t-gp-start-1step.txt" "$O"/mut1-*.txt | cut -c1-240
grep -h "^GN" "$O/gnprobe-1step.txt"
