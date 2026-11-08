#!/usr/bin/env bash
# Review of lane fitedge: the five named test files on the lane library, one
# R process each, sequentially (one process of the review's three).
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/review-fitedge-out"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true R_PARALLEL_PORT=11831
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-fitedge-lib
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/stan-cache
for f in ${FILES:-test-gp-start test-latres test-nl-family test-se-check test-gp-noncentered}; do
  "$RS" "$WT/dev/fitedge-run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
    > "$OUT/t-$f-lane.txt" 2>&1
done
grep -h "^lib\|^RESULT\|^DETAIL" "$OUT"/t-*-lane.txt
echo TESTS-DONE
