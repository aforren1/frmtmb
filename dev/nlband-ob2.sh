#!/usr/bin/env bash
# Lane nlband, fix round 2: test-fit-restart.R and the round's 21
# roundings (dev/nlband-fix2ulp.R) under OpenBLAS 0.3.26 with LAPACK,
# OPENBLAS_CORETYPE=Haswell, 4 threads; one process at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
export NLBAND_ARM=after
OB="dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
"$OB" dev/nlstop-obrun.R frmtmb tests/testthat/test-fit-restart.R \
  > dev/nlband-out/ob2-0326-fit-restart.txt 2>&1
"$OB" dev/nlband-fix2ulp.R > dev/nlband-out/fix2ulp-ob0326.txt 2>&1
echo "OB2 DONE"
