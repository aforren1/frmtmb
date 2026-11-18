#!/usr/bin/env bash
# Lane nlband, fix round 1: the nonlinear test files (dev/nlband-nl.jobs)
# under the OpenBLAS 0.3.26 copy with LAPACK, one process at a time, brms
# fit tests on.
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
OB="dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
mkdir -p dev/nlband-out/obnl
while read -r p f; do
  "$OB" dev/nlstop-obrun.R "$p" "$f" > "dev/nlband-out/obnl/$(basename "$f" .R).txt" 2>&1
done < dev/nlband-nl.jobs
echo "OBNL DONE"
