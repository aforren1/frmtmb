#!/usr/bin/env bash
# Lane gpfield: the test files this lane changed or that reach its path, one
# process each, OB_N at a time (default 2), under an OpenBLAS copy of R built
# by dev/ciharden-openblas.sh. Args: <ver> (lapack mode). Each log opens with
# the R home (the copy) and ends with its RESULT line.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-gpfield
RS="$WT/dev/ciharden-out/Rob$1-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpfield-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/gpfield-out/ob$1
mkdir -p "$OUT"
for f in test-gp-kernels test-gp-field test-gp-predict-draw test-gp-start test-gp-by test-gp-multidim test-gp-noncentered test-gp-ar-start test-gp-inner-pd test-predict-re-uncertainty; do
  while [ "$(jobs -rp | wc -l)" -ge "${OB_N:-2}" ]; do sleep 1; done
  ( "$RS" -e "cat('blas:', R.home(), '\n')" > "$OUT/$f.txt" 2>&1
    "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
      >> "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^blas\|^lib\|^RESULT" "$OUT"/*.txt
