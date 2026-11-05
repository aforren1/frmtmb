#!/usr/bin/env bash
# The numerical test files under OpenBLAS with LAPACK routed to it, 0.3.26 then
# 0.3.32, one file per process. Builds the R copies; delete dev/ciharden-out after.
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/gpgap-out/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/gpgap-out/ob32-build.txt 2>&1
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
for v in 26 32; do
  RB=dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe
  rm -rf dev/gpgap-out/ob$v; mkdir -p dev/gpgap-out/ob$v
  for spec in "frmtmb tests/testthat/test-gp-nugget-floor.R" \
              "frmtmb tests/testthat/test-brms-lik-gp-exact.R" \
              "frmtmb tests/testthat/test-gp-multidim.R" \
              "frmtmb tests/testthat/test-gp-by.R" \
              "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-reparam.R"; do
    set -- $spec
    while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 2; done
    ( "$RB" dev/gpgap-runtest.R "$1" "$2" > "dev/gpgap-out/ob$v/$(basename "$2" .R).txt" 2>&1 ) &
  done
  wait
done
