#!/usr/bin/env bash
# Lane wcheck2: one test file in one R process.
#   bash dev/wcheck2-tf.sh <lane|base> <package> <test file name> <log>
# "lane" puts the lane library first (FRMTMB_RELLIB); "base" runs
# rellib-r21 alone. The test file is read from this worktree.
WT=C:/Users/adf44/source/r/frmtmb-wt-wcheck2
arm=$1; pkg=$2; f=$3; log=$4
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
if [ "$arm" = lane ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-wcheck2-lib
else
  unset FRMTMB_RELLIB
fi
if [ "$pkg" = frmtmb ]; then dir=$WT/tests/testthat; else
  dir=$WT/extensions/$pkg/tests/testthat; fi
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" $WT/dev/release/run-tests.R \
  "$pkg" "$dir/$f" > "$log" 2>&1
grep -a "^RESULT\|^lib:" "$log"
