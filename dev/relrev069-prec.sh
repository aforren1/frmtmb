#!/usr/bin/env bash
# Reviewer of the 0.69.0 consolidation: rerun the files holding the 8
# fits the precedence run (dev/rel069-log/prof-sum.txt) found with two
# kinds of report, with dev/relrev069-profile.R checking each report.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$ROOT"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/relrev069-stan-cache
export R_PROFILE_USER="$ROOT/dev/relrev069-profile.R"
O=$ROOT/dev/relrev069-log/prec; mkdir -p "$O"
run() {
  p="$1"; f="$2"; t="$3"; nm="$p--$(basename "$f" .R)"
  rm -f "$O/$nm.txt"
  RELREV_TARGETS="$t" RELREV_OUT="$O/$nm.txt" \
    "$RS" dev/ciharden-run1.R base "$p" "$f" > "$O/$nm.log" 2>&1 &
}
[ -n "${ONLY:-}" ] || run frmtmb tests/testthat/test-brms-suite-emmeans.R 2
run frmtmb tests/testthat/test-brms-suite-methods.R 2,10,15
run frmtmb tests/testthat/test-diagnostics-ux.R 58
[ -n "${ONLY:-}" ] && [ -z "${FUZZ:-}" ] || run frmtmb tests/testthat/test-fuzz.R 231,232
[ -n "${ONLY:-}" ] || run frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-brms-suite-methods.R 3
wait
echo PREC DONE
