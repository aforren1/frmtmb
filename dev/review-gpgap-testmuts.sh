#!/usr/bin/env bash
# item 7: run the two new test files on each mutant, 4 at once.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/review-gpgap-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
LIB="C:/Users/adf44/source/r/rev-gpgap-lib;C:/Users/adf44/source/r/rellib-r16"
O=dev/review-gpgap-out
for m in none nomsg thr_1e3 thr_1e4 thr_05 thr_2 no_dpar all_types no_se no_z sig_max noise_wrong class_drop; do
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 1; done
  REV_MUT=$m "$RS" dev/review-gpgap-testmut.R test-gp-nugget-floor.R "$LIB" > $O/tm-floor-$m.txt 2>&1 &
done
for m in none kern_102 nug_2e6 nug_1e7 ord_none ord_direct ord_col1 ord_always; do
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 1; done
  REV_MUT=$m "$RS" dev/review-gpgap-testmut.R test-brms-lik-gp-exact.R "$LIB" > $O/tm-lik-$m.txt 2>&1 &
done
wait
# on the base build, the lane's tests as they stand
"$RS" dev/review-gpgap-testmut.R test-gp-nugget-floor.R "C:/Users/adf44/source/r/rellib-r16" "tests/testthat" > $O/tm-floor-r16.txt 2>&1
