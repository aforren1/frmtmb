#!/usr/bin/env bash
# Review of lane priors2: four test files under OpenBLAS 0.3.32 + LAPACK,
# once with OpenBLAS's own core detection and once with
# OPENBLAS_CORETYPE=Haswell, 2 processes at a time.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-priors2
RS="$WT/dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-priors2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/reviews/priors2-review-ob-out
mkdir -p "$OUT"
for arm in auto Haswell; do
  for f in test-prior-bounds-ends test-prior-mixture-theta test-mixture-starts test-prior-table-brms; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
    (
      if [ "$arm" = Haswell ]; then export OPENBLAS_CORETYPE=Haswell; fi
      "$RS" -e "cat('blas:', R.home(), extSoftVersion()[['BLAS']], Sys.getenv('OPENBLAS_CORETYPE'), '\n'); print(La_library())" > "$OUT/$arm-$f.txt" 2>&1
      "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" >> "$OUT/$arm-$f.txt" 2>&1
    ) &
  done
done
wait
for x in "$OUT"/*.txt; do echo "== $(basename $x)"; grep -h "^blas\|^lib\|^RESULT" "$x"; done
