#!/usr/bin/env bash
# Lane gpstop: the test files this lane changed or that reach its path, one
# process each, 3 at a time, under an OpenBLAS copy of R built by
# dev/ciharden-openblas.sh. Args: <ver> (lapack mode). Each log opens with
# the R home (the copy) and ends with its RESULT line.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-gpstop
RS="$WT/dev/ciharden-out/Rob$1-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpstop-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/gpstop-out/ob$1
mkdir -p "$OUT"
for f in test-gp-start test-gp-ar-start test-gp-noncentered test-gp-by test-gp-multidim test-nl-family; do
  while [ "$(jobs -rp | wc -l)" -ge "${OB_N:-3}" ]; do sleep 1; done
  ( "$RS" -e "cat('blas:', R.home(), '\n')" > "$OUT/$f.txt" 2>&1
    "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
      >> "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^blas\|^lib\|^RESULT" "$OUT"/*.txt
