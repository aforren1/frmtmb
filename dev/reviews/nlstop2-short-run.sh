#!/usr/bin/env bash
# Short check of lane nlstop2's fix round. Usage: nlstop2-short-run.sh <set>
#   a: test-fit-restart.R listed on reference BLAS (lane), 0.3.32 Haswell
#      (lane) and reference BLAS (rellib-r21), 3 processes at once
#   b: test-bracket-access.R and test-families.R on the lane, 2 at once
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlstop2/dev/stan-cache
O=dev/reviews/nlstop2-short-out
mkdir -p $O
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
S=dev/reviews/nlstop2-short-list.R
T=tests/testthat
case "$1" in
a)
  "$RS" $S $T/test-fit-restart.R > $O/fit-restart-ref.txt 2>&1 &
  OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell \
    dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe $S \
    $T/test-fit-restart.R > $O/fit-restart-ob0332.txt 2>&1 &
  NLSTOP2_SHORT_ARM=before "$RS" $S $T/test-fit-restart.R \
    > $O/fit-restart-before.txt 2>&1 &
  wait ;;
b)
  "$RS" $S $T/test-bracket-access.R > $O/bracket-access-ref.txt 2>&1 &
  "$RS" $S $T/test-families.R > $O/families-ref.txt 2>&1 &
  wait ;;
esac
echo "SET $1 DONE"
