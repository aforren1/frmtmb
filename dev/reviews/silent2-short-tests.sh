#!/usr/bin/env bash
# Short check of lane silent2's fix round: six test files on the lane
# library, three processes at a time.
#   bash dev/reviews/silent2-short-tests.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-silent2
cd "$ROOT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
O=dev/reviews/silent2-short-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
run() {
  "$R" dev/reviews/silent2-short-runtest.R "$1" "$2" > "$O/t-$(basename "$2" .R).txt" 2>&1
}
export -f run
export R O
printf '%s\n' \
  "frmtmb tests/testthat/test-ordinal-mixture.R" \
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-laplace-draws.R" \
  "frmtmb tests/testthat/test-predict-newdata.R" \
  "frmtmb tests/testthat/test-newdata-aterms.R" \
  "frmtmb tests/testthat/test-subset-rate.R" \
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-brms-draws-order.R" |
  xargs -P 3 -L 1 bash -c 'run "$0" "$1"'
grep -h "^RESULT\|^lib:" $O/t-*.txt
echo SHORT-TESTS-DONE
