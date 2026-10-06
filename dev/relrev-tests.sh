#!/usr/bin/env bash
# Reviewer, item 9: rerun a sample of test files on rellib-r6, gated,
# one R process per file, all at once. Logs dev/relrev-log/t-*.txt.
set -u
rel=/c/Users/adf44/source/r/frmtmb-wt-release
S=/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$S/relrev-stan-cache"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP unset}"
T=$rel/tests/testthat
X=$rel/extensions
O=$rel/dev/relrev-log
run() { tag=$1; shift; ( "$RS" $rel/dev/relrev-run1.R "$@" > "$O/t-$tag.txt" 2>&1 ) & }
for f in test-cumulative-cs test-ce-predict-median test-ordinal-thres-names \
         test-compat test-ordinal-mixture test-gp-by; do
  run "$f" frmtmb "$T/$f.R"
done
run test-gp-by-globalf frmtmb "$T/test-gp-by.R" globalf
run sample-test-incl-thres-draws frmtmb.sample \
  "$X/frmtmb.sample/tests/testthat/test-incl-thres-draws.R"
run spline-test-curve frmtmb.spline "$X/frmtmb.spline/tests/testthat/test-curve.R"
wait
grep -h "^RESULT" $O/t-*.txt
