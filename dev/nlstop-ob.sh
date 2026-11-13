#!/usr/bin/env bash
# Lane nlstop: the lane's nonlinear and prior test files under an
# OpenBLAS copy of R, three at a time. Usage: bash dev/nlstop-ob.sh <ver>
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop || exit 1
V="$1"
OB="dev/ciharden-out/Rob$V-lapack/bin/x64/Rscript.exe"
O="dev/nlstop-out/ob-$V"
mkdir -p "$O"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
files="test-fit-restart.R test-nl-start.R test-nl.R test-nl-family.R test-nlf.R test-map.R $(cd tests/testthat && ls test-prior-*.R test-priors-*.R)"
n=0
for f in $files; do
  ( "$OB" dev/nlstop-obrun.R frmtmb "tests/testthat/$f" > "$O/${f%.R}.txt" 2>&1 ) &
  n=$((n + 1))
  if [ $n -ge 3 ]; then wait -n; n=$((n - 1)); fi
done
wait
echo "OB DONE $V"
