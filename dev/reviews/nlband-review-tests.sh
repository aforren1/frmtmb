#!/usr/bin/env bash
# Review of lane nlband, check 7: the final test-fit-restart.R on
# rellib-r22 (before; FRMTMB_RELLIB unset) and on the lane library (after).
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=dev/reviews/nlband-review-out
unset FRMTMB_RELLIB
"$RS" dev/release/run-tests.R frmtmb tests/testthat/test-fit-restart.R \
  > $OUT/t-fit-restart-r22.txt 2>&1
FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib "$RS" \
  dev/release/run-tests.R frmtmb tests/testthat/test-fit-restart.R \
  > $OUT/t-fit-restart-lane.txt 2>&1
echo "TESTS DONE"
