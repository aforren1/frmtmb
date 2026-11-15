#!/usr/bin/env bash
# Review of lane quadseam2, attack 9: OpenBLAS 0.3.32 with LAPACK
# (dev/ciharden-openblas.sh) under OPENBLAS_CORETYPE=Haswell and
# Sandybridge, on test-se-tier.R, test-quadrature-defects.R and
# test-latres-quadratic.R, lane library; then dev/ciharden-out/ removed.
# Logs dev/reviews/quadseam2-review-log/ob-<coretype>-<file>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
cd "$WT"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-quadseam2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
L=dev/reviews/quadseam2-review-log
T=tests/testthat
bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob-build.txt 2>&1
RO="dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe"
for ct in Haswell Sandybridge; do
  export OPENBLAS_CORETYPE=$ct
  for f in test-se-tier.R test-quadrature-defects.R test-latres-quadratic.R; do
    "$RO" -e "cat('blas:', extSoftVersion()[['BLAS']], 'coretype:', Sys.getenv('OPENBLAS_CORETYPE'), '\n')" \
      > "$L/ob-$ct-$(basename $f .R).txt" 2>&1
    "$RO" dev/release/run-tests.R frmtmb $T/$f >> "$L/ob-$ct-$(basename $f .R).txt" 2>&1 &
  done
  wait
done
grep -h "^RESULT\|^lib\|^blas" $L/ob-*-test-*.txt | cut -c1-200
tail -3 $L/ob-build.txt
rm -rf dev/ciharden-out
echo "OB DONE"
