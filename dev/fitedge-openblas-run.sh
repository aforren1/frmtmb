#!/usr/bin/env bash
# The changed test files, one process each, sequentially, under an OpenBLAS
# copy of R built by dev/ciharden-openblas.sh. Args: <ver> (lapack mode).
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
RS="$WT/dev/ciharden-out/Rob$1-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-fitedge-lib
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/stan-cache
for f in test-gp-start test-latres test-nl-family test-se-check; do
  "$RS" "$WT/dev/fitedge-run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
    > "$WT/dev/fitedge-out/ob$1-$f.txt" 2>&1
done
grep -h "^lib\|^RESULT\|^DETAIL" "$WT"/dev/fitedge-out/ob$1-*.txt
