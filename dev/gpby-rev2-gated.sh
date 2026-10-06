#!/bin/sh
# Reviewer round 2: every core test file, gated, one R process per file,
# 14 at a time, on the lane build. Logs in dev/gpby-rev2-gated/.
W=/c/Users/adf44/source/r/frmtmb-wt-gpby
cd $W
mkdir -p dev/gpby-rev2-gated
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
ls tests/testthat/test-*.R | /usr/bin/xargs -P 14 -n 1 sh dev/gpby-rev2-one.sh
echo ALLDONE
