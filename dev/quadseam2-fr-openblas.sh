#!/usr/bin/env bash
# Fix round copy: OPENBLAS_CORETYPE=Haswell, logs dev/quadseam2-log/fr-ob<ver>-*.txt
# OpenBLAS 0.3.26 and 0.3.32 with LAPACK (dev/ciharden-openblas.sh) on the
# test files this lane added or changed, three processes at a time, then
# dev/ciharden-out/ removed. Logs dev/quadseam2-log/ob<ver>-<file>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
cd "$WT"
export OPENBLAS_CORETYPE=Haswell NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-quadseam2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
L=dev/quadseam2-log
T=tests/testthat
LT=extensions/frmtmb.latent/tests/testthat
for v in 0.3.26 0.3.32; do
  # the build's own BLAS check is an R process: run it alone
  bash dev/ciharden-openblas.sh $v lapack > $L/fr-ob$v-build.txt 2>&1
  RO="dev/ciharden-out/Rob$v-lapack/bin/Rscript.exe"
  run() {
    "$RO" dev/release/run-tests.R "$1" "$2" > "$L/fr-ob$v-$(basename "$2" .R).txt" 2>&1
  }
  run frmtmb $T/test-quadrature-defects.R &
  run frmtmb $T/test-stationary-passes.R &
  run frmtmb $T/test-latres-quadratic.R &
  wait
  run frmtmb $T/test-quad-escape-record.R &
  run frmtmb $T/test-se-tier.R &
  run frmtmb.latent $LT/test-hmm-stationary.R &
  wait
  run frmtmb $T/test-skew-normal-start.R &
  run frmtmb $T/test-quadrature-nested.R &
  run frmtmb $T/test-se-limit.R &
  wait
  grep -h "^RESULT\|^lib" $L/fr-ob$v-test-*.txt | cut -c1-200
  tail -3 $L/fr-ob$v-build.txt
done
rm -rf dev/ciharden-out
echo "OB DONE"
