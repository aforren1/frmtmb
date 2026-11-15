#!/usr/bin/env bash
# Short check of lane quadseam2's fix round, item 6: the changed and
# neighboring test files on the lane library, reference BLAS, one file
# per R process, two at a time. Logs
# dev/reviews/quadseam2-short-log/test-<file>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
cd "$WT"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-quadseam2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/reviews/quadseam2-short-log
T=tests/testthat
run() {
  R_PARALLEL_PORT=$((11000 + RANDOM % 5000)) "$RS" dev/release/run-tests.R \
    "$1" "$2" > "$L/test-$(basename "$2" .R).txt" 2>&1
}
run frmtmb $T/test-quadrature-defects.R &
run frmtmb $T/test-se-tier.R &
wait
run frmtmb $T/test-latres.R &
run frmtmb $T/test-quad-escape-record.R &
wait
run frmtmb $T/test-stationary-passes.R &
run frmtmb $T/test-latres-quadratic.R &
wait
run frmtmb.latent extensions/frmtmb.latent/tests/testthat/test-hmm-stationary.R
grep -h "^RESULT\|^lib" $L/test-*.txt | cut -c1-220
echo "TESTS DONE"
