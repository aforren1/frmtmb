#!/usr/bin/env bash
# test-gp-start.R again under OpenBLAS after its last row, 0.3.26 and
# 0.3.32. Builds the R copies; deletes them after.
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/ncgp-out/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/ncgp-out/ob32-build.txt 2>&1
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
for v in 26 32; do
  dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb \
    tests/testthat/test-gp-start.R > "dev/ncgp-out/ob$v/test-gp-start.txt" 2>&1
done
rm -rf dev/ciharden-out
echo OB2-DONE
