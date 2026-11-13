#!/usr/bin/env bash
# run quadseam-short-hetrace.R on one test file: <arm> <file>
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/stan-cache
cd /c/Users/adf44/source/r/frmtmb-wt-quadseam/dev/reviews
arm=$1; f=$2; b=$(basename $f .R)
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" quadseam-short-hetrace.R $arm \
  C:/Users/adf44/source/r/frmtmb-wt-quadseam/tests/testthat/$f \
  quadseam-short-log/he-$arm-$b.txt > quadseam-short-log/he-$arm-$b.out 2>&1
echo "done $f $?"
