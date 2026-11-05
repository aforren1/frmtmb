#!/usr/bin/env bash
# OpenBLAS on the final build, the files the last change can reach, then the checks.
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
for v in 26 32; do
  RB=dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe
  for spec in "frmtmb tests/testthat/test-gp-noncentered.R" \
              "frmtmb tests/testthat/test-gp-multidim.R" \
              "frmtmb tests/testthat/test-gp-by.R" \
              "frmtmb tests/testthat/test-curve.R" \
              "frmtmb tests/testthat/test-curve-difference.R" \
              "frmtmb tests/testthat/test-fd-chain.R" \
              "frmtmb tests/testthat/test-emmeans.R" \
              "frmtmb tests/testthat/test-predict-re-uncertainty.R"; do
    set -- $spec
    while [ "$(jobs -r | wc -l)" -ge 3 ]; do sleep 2; done
    ( "$RB" dev/ncgp-runtest.R "$1" "$2" > "dev/ncgp-out/ob$v/$(basename "$2" .R).txt" 2>&1 ) &
  done
  wait
done
echo OB2-DONE
unset FRMTMB_BRMS_FIT_TESTS OPENBLAS_NUM_THREADS
export -n PATH
bash dev/ncgp-check.sh . > dev/ncgp-out/check-core.txt 2>&1
bash dev/ncgp-check.sh extensions/frmtmb.sample > dev/ncgp-out/check-sample.txt 2>&1
echo CHECK-DONE
