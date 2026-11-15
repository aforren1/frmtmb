#!/usr/bin/env bash
# Lane nlstop2, fix round: test-fit-restart.R and the explin draws
# (dev/nlstop2-explintest.R) under both OpenBLAS copies with LAPACK
# (dev/ciharden-openblas.sh <ver> lapack), OPENBLAS_CORETYPE=Haswell as
# the review ran them. One process at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop2-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlstop2/dev/stan-cache
export NLSTOP2_ARM=after
for V in 0.3.26 0.3.32; do
  OB="dev/ciharden-out/Rob$V-lapack/bin/x64/Rscript.exe"
  "$OB" dev/nlstop-obrun.R frmtmb tests/testthat/test-fit-restart.R \
    > "dev/nlstop2-out/ob2-$V-fit-restart.txt" 2>&1
  "$OB" dev/nlstop2-explintest.R > "dev/nlstop2-out/ob2-$V-explintest.txt" 2>&1
done
echo "OB2 DONE"
