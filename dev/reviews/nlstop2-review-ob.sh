#!/usr/bin/env bash
# Review of lane nlstop2: one test file, or one review script, under an
# OpenBLAS copy of R with a given core type.
# Usage: bash dev/reviews/nlstop2-review-ob.sh <ver> <coretype|-> <file> <log>
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
V="$1"; CT="$2"; F="$3"; LOG="$4"
OB="dev/ciharden-out/Rob$V-lapack/bin/x64/Rscript.exe"
mkdir -p dev/reviews/nlstop2-review-out
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
[ "$CT" != "-" ] && export OPENBLAS_CORETYPE="$CT"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop2-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
case "$F" in
  tests/*) "$OB" dev/nlstop-obrun.R frmtmb "$F" > "dev/reviews/nlstop2-review-out/$LOG" 2>&1 ;;
  *) "$OB" "$F" > "dev/reviews/nlstop2-review-out/$LOG" 2>&1 ;;
esac
echo "EXIT $? $V $CT $F" >> "dev/reviews/nlstop2-review-out/$LOG"
