#!/usr/bin/env bash
# Lane drawdir: dev/drawdir-ulp.R and the two draw test
# files under OpenBLAS 0.3.26 and 0.3.32 with LAPACK, core Haswell, one R
# process at a time. Builds the R copies with dev/ciharden-openblas.sh.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drawdir
cd "$ROOT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
O=dev/drawdir-ob-out
mkdir -p $O
for v in 0.3.26 0.3.32; do
  bash dev/ciharden-openblas.sh $v lapack > $O/build-$v.txt 2>&1
  RS="$ROOT/dev/ciharden-out/Rob$v-lapack/bin/x64/Rscript.exe"
  "$RS" dev/drawdir-ulp.R C:/Users/adf44/source/r/wt-drawdir-lib "ob$v-Haswell" > $O/ulp-$v.txt 2>&1
  for f in test-predict-draw-check test-predict-draw-edge; do
    "$RS" dev/drawdir-run-tests.R frmtmb tests/testthat/$f.R > $O/$v-$f.txt 2>&1
  done
done
echo done
