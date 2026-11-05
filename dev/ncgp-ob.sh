#!/usr/bin/env bash
# The numerical test files under OpenBLAS with LAPACK routed to it, 0.3.26 then
# 0.3.32, one file per process. Builds the R copies; delete dev/ciharden-out after.
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/ncgp-out/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/ncgp-out/ob32-build.txt 2>&1
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
for v in 26 32; do
  RB=dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe
  rm -rf dev/ncgp-out/ob$v; mkdir -p dev/ncgp-out/ob$v
  for spec in "frmtmb tests/testthat/test-gp-noncentered.R" \
              "frmtmb tests/testthat/test-gp-start.R" \
              "frmtmb tests/testthat/test-emmeans.R" \
              "frmtmb tests/testthat/test-se-check.R" \
              "frmtmb tests/testthat/test-brms-likelihood.R" \
              "frmtmb tests/testthat/test-gp-multidim.R" \
              "frmtmb tests/testthat/test-gp-by.R" \
              "frmtmb tests/testthat/test-curve.R" \
              "frmtmb tests/testthat/test-curve-deriv.R" \
              "frmtmb tests/testthat/test-curve-difference.R" \
              "frmtmb tests/testthat/test-fd-chain.R" \
              "frmtmb tests/testthat/test-covstruct-floor.R" \
              "frmtmb tests/testthat/test-predict-re-uncertainty.R" \
              "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-gp-by-draws.R" \
              "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-reparam.R"; do
    set -- $spec
    while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 2; done
    ( "$RB" dev/ncgp-runtest.R "$1" "$2" > "dev/ncgp-out/ob$v/$(basename "$2" .R).txt" 2>&1 ) &
  done
  wait
done
echo OB-DONE
