#!/usr/bin/env bash
# Lane gpselimit: the test files this lane changed or that reach its path,
# one process each, OB_N at a time (default 2), under an OpenBLAS copy of R
# built by dev/ciharden-openblas.sh <ver> lapack. Each log opens with the R
# home (the copy) and ends with its RESULT line.
# Usage: bash dev/gpselimit-openblas-run.sh <ver>
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-gpselimit
RS="/c/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/ciharden-out/Rob$1-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpselimit-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
OUT=/c/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-out/ob$1
rm -rf "$OUT"; mkdir -p "$OUT"
for f in test-gp-white-limit test-gp-field test-gp-by test-gp-multidim \
         test-marginaleffects-re test-se-tier test-se-check test-se-limit \
         test-gp-start test-gp-kernels; do
  while [ "$(jobs -rp | wc -l)" -ge "${OB_N:-2}" ]; do sleep 1; done
  ( "$RS" -e "cat('blas:', R.home(), extSoftVersion()[['BLAS']], '\n')" \
      > "$OUT/$f.txt" 2>&1
    "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
      >> "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^blas\|^lib\|^RESULT" "$OUT"/*.txt
