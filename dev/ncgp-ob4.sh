#!/usr/bin/env bash
# test-brms-likelihood.R after row 10a's tolerance: both OpenBLAS at 1 and 4 threads,
# and the reference BLAS. Deletes the R copies after.
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
F=tests/testthat/test-brms-likelihood.R
for v in 26 32; do
  for th in 1 4; do
    OPENBLAS_NUM_THREADS=$th dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb $F > "dev/ncgp-out/ob$v/test-brms-likelihood-t$th.txt" 2>&1 &
  done
done
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/ncgp-runtest.R frmtmb $F > dev/ncgp-out/gated/frmtmb--test-brms-likelihood.txt 2>&1 &
wait
cp dev/ncgp-out/ob26/test-brms-likelihood-t4.txt dev/ncgp-out/ob26/test-brms-likelihood.txt
cp dev/ncgp-out/ob32/test-brms-likelihood-t4.txt dev/ncgp-out/ob32/test-brms-likelihood.txt
rm -rf dev/ciharden-out
echo OB4-DONE > dev/ncgp-out/ob4.done
