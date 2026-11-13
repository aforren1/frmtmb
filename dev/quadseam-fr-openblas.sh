#!/usr/bin/env bash
# Fix round: OpenBLAS 0.3.26 and 0.3.32 with LAPACK (dev/ciharden-openblas.sh)
# on the three files that test the changed paths, three processes at a
# time, then dev/ciharden-out/ removed. Logs dev/quadseam-log/fr-ob<ver>-<file>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam
cd "$WT"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-quadseam-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
for v in 0.3.26 0.3.32; do
  bash dev/ciharden-openblas.sh $v lapack > dev/quadseam-log/fr-ob$v-build.txt 2>&1
  RO="dev/ciharden-out/Rob$v-lapack/bin/Rscript.exe"
  for f in test-quadrature-defects test-quadrature-nested test-se-tier; do
    "$RO" dev/release/run-tests.R frmtmb tests/testthat/$f.R \
      > dev/quadseam-log/fr-ob$v-$f.txt 2>&1 &
  done
  wait
  grep -h "RESULT\|DETAIL\|^lib" dev/quadseam-log/fr-ob$v-*.txt | cut -c1-200
  tail -3 dev/quadseam-log/fr-ob$v-build.txt
done
rm -rf dev/ciharden-out
echo "OB DONE"
