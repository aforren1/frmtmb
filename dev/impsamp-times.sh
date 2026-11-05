#!/usr/bin/env bash
# Run time before and after on the other test files that use importance =,
# each pair at once so they share the load. A file whose test the lane
# rewrote runs its original on the base build.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-impsamp-lib
O=dev/impsamp-out/oldtests
export NOT_CRAN=true
pair() {
  p="$1"; before="$2"; after="$3"; lab=$(basename "$after" .R)
  "$R" dev/impsamp-testfile.R "" "$p" "$before" > "dev/impsamp-out/tb-$lab.txt" 2>&1 &
  "$R" dev/impsamp-testfile.R "$L" "$p" "$after" > "dev/impsamp-out/ta-$lab.txt" 2>&1 &
  wait
  b=$(grep -h "^RESULT" "dev/impsamp-out/tb-$lab.txt" | sed 's/.*secs=//')
  a=$(grep -h "^RESULT" "dev/impsamp-out/ta-$lab.txt" | sed 's/.*secs=//')
  echo "$lab before $b after $a"
}
pair frmtmb $O/test-grad-verdict.R tests/testthat/test-grad-verdict.R
pair frmtmb $O/test-se-exact-hessian.R tests/testthat/test-se-exact-hessian.R
pair frmtmb $O/test-rl-example.R tests/testthat/test-rl-example.R
T=extensions/frmtmb.learn/tests/testthat
pair frmtmb.learn $T/test-rlddm-ndt.R $T/test-rlddm-ndt.R
pair frmtmb.learn $T/test-session.R $T/test-session.R
pair frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-car.R \
  extensions/frmtmb.spatial/tests/testthat/test-car.R
