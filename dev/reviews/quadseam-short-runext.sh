#!/usr/bin/env bash
# trace obj_he_cached() over an extension's test files, one file per process, sequentially: <pkg>
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/stan-cache
cd /c/Users/adf44/source/r/frmtmb-wt-quadseam/dev/reviews
pkg=$1
for f in ../../extensions/$pkg/tests/testthat/test-*.R; do
  b=$(basename $f .R)
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" quadseam-short-hetrace-ext.R lane \
    C:/Users/adf44/source/r/frmtmb-wt-quadseam/extensions/$pkg/tests/testthat/$b.R \
    quadseam-short-log/hx-$pkg-$b.txt $pkg > quadseam-short-log/hx-$pkg-$b.out 2>&1
  echo "done $pkg $b $?"
done
