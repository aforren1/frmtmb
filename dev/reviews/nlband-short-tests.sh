#!/usr/bin/env bash
# Short check of lane nlband: test-fit-restart.R on the lane library.
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_PARALLEL_PORT=11931
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib \
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/release/run-tests.R \
  frmtmb tests/testthat/test-fit-restart.R \
  > dev/reviews/nlband-short-out/t-fit-restart.txt 2>&1
echo "EXIT $?" >> dev/reviews/nlband-short-out/t-fit-restart.txt
