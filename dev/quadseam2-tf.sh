#!/usr/bin/env bash
# One test file per process: bash dev/quadseam2-tf.sh <base|lib> <pkg> <file> <tag>
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
LIB=C:/Users/adf44/source/r/wt-quadseam2-lib
arg=$1; [ "$1" = lib ] && arg=$LIB
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=${FRMTMB_BRMS_FIT_TESTS:-true}
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-quadseam2/dev/stan-cache
export R_PARALLEL_PORT=$((11000 + RANDOM % 5000))
mkdir -p "$WT/dev/quadseam2-log/tests"
cd "$WT"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/quadseam2-test.R "$arg" "$2" "$3" \
  > "dev/quadseam2-log/tests/$4.txt" 2>&1
grep -E "^(RESULT|FAIL|lib:)" "dev/quadseam2-log/tests/$4.txt" || tail -5 "dev/quadseam2-log/tests/$4.txt"
