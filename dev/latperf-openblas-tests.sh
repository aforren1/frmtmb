#!/usr/bin/env bash
# Lane latperf: the test files this lane's numerical change reaches, under
# the R copy of dev/ciharden-openblas.sh (OpenBLAS BLAS and LAPACK), one file
# per process, one at a time. Usage: bash dev/latperf-openblas-tests.sh <ver>
set -u
VER="${1:-0.3.26}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="$ROOT/dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
OUT="$ROOT/dev/latperf-out/openblas-$VER"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-latperf-lib
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
mkdir -p "$OUT"
for f in test-latres-sparse test-latres test-latres-newrow \
         test-predict-lp-basis test-curve test-latres-brms; do
  while [ "$(jobs -rp | wc -l)" -ge 1 ]; do sleep 1; done
  ( "$RS" "$ROOT/dev/latperf-run-tests.R" frmtmb \
      "$ROOT/tests/testthat/$f.R" > "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^lib:\|^RESULT\|^DETAIL" "$OUT"/*.txt
