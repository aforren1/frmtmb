#!/usr/bin/env bash
# Lane nlband: test-fit-restart.R and the ulp draws (dev/nlband-ulp.R)
# under both OpenBLAS copies with LAPACK (dev/ciharden-openblas.sh <ver>
# lapack), OPENBLAS_CORETYPE=Haswell (this CPU has no AVX-512) and 4
# threads, as the runners have 4 cores. One process at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
export NLBAND_ARM=after
for V in 0.3.26 0.3.32; do
  OB="dev/ciharden-out/Rob$V-lapack/bin/x64/Rscript.exe"
  "$OB" dev/nlstop-obrun.R frmtmb tests/testthat/test-fit-restart.R \
    > "dev/nlband-out/ob-$V-fit-restart.txt" 2>&1
  "$OB" -e 'cat("blas:", extSoftVersion()[["BLAS"]], "|", La_library(), "\n"); source("dev/nlband-ulp.R")' \
    > "dev/nlband-out/ob-$V-ulp.txt" 2>&1
done
echo "OB DONE"
