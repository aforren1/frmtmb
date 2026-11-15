#!/usr/bin/env bash
# Lane nlstop2: the lane's test files under an OpenBLAS copy of R
# (dev/ciharden-openblas.sh <ver> lapack), N at a time, one file per
# process. Usage: bash dev/nlstop2-ob.sh <ver> [N]
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
V="$1"
N="${2:-1}"
OB="dev/ciharden-out/Rob$V-lapack/bin/x64/Rscript.exe"
O="dev/nlstop2-out/ob-$V"
mkdir -p "$O"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop2-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlstop2/dev/stan-cache
files="test-fit-restart.R test-nl-start.R test-nl.R test-nl-family.R test-nlf.R test-families.R test-numerical-robustness.R test-simulate-density.R test-update-pool.R"
n=0
for f in $files; do
  ( "$OB" dev/nlstop-obrun.R frmtmb "tests/testthat/$f" > "$O/${f%.R}.txt" 2>&1 ) &
  n=$((n + 1))
  if [ $n -ge "$N" ]; then wait -n; n=$((n - 1)); fi
done
wait
echo "OB DONE $V"
