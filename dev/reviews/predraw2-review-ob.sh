#!/usr/bin/env bash
# Review of lane predraw2: rev-b and the two draw test files under
# OpenBLAS 0.3.26 / 0.3.32 (with LAPACK), cores Haswell and Sandybridge,
# one R process at a time.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-predraw2
cd "$ROOT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 R_PARALLEL_PORT=11875
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
O=dev/reviews/predraw2-review-ob-out
mkdir -p $O
for v in 0.3.26 0.3.32; do
  RS="$ROOT/dev/ciharden-out/Rob$v-lapack/bin/Rscript.exe"
  for c in Haswell Sandybridge; do
    export OPENBLAS_CORETYPE=$c
    "$RS" dev/reviews/predraw2-review-b.R "ob$v-$c" > $O/b-$v-$c.txt 2>&1
    for f in test-predict-draw-check.R test-predict-draw-edge.R; do
      "$RS" dev/reviews/predraw2-review-run.R frmtmb tests/testthat/$f > $O/t-$v-$c-$f.txt 2>&1
    done
  done
done
echo done
