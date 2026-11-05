#!/usr/bin/env bash
# Fix round 1's final evidence, on the final build. Library paths are absolute in each
# script; read `lib:` in each log. Nothing in R/, man/, NEWS or the tests changes after
# this starts.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-ncgp-lib
O=dev/ncgp-out
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"

rm -rf $O/gated
FRMTMB_BRMS_FIT_TESTS=true bash dev/ncgp-par.sh $O/gated frmtmb frmtmb.sample frmtmb.spatial frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.survival
grep -h '^RESULT' $O/gated/*.txt > $O/gated.txt
(
  export NOT_CRAN=true
  export NCGP_TLIB="C:/Users/adf44/source/r/rellib-r16"
  "$RS" dev/ncgp-runtest.R frmtmb tests/testthat/test-gp-start.R > $O/seen-fail-start.txt 2>&1
  "$RS" dev/ncgp-runtest.R frmtmb tests/testthat/test-gp-noncentered.R > $O/seen-fail.txt 2>&1
)
echo STAGE-GATED

bash dev/ncgp-l1-rest.sh > $O/l1-rest.log 2>&1
bash dev/ncgp-l1eq.sh > $O/l1eq.log 2>&1
bash dev/ncgp-l1ng.sh > $O/l1ng.log 2>&1
"$RS" dev/ncgp-row10.R > $O/row10.txt 2>&1 &
"$RS" dev/ncgp-readers.R > $O/readers.txt 2>&1 &
"$RS" dev/ncgp-grcheck.R > $O/grcheck.txt 2>&1 &
"$RS" dev/ncgp-prior.R > $O/prior.txt 2>&1 &
"$RS" dev/ncgp-cholscan.R > $O/cholscan.txt 2>&1 &
wait
"$RS" dev/ncgp-fcwarn.R > $O/fcwarn.txt 2>&1
"$RS" dev/review-ncgp-modes2.R lane > $O/modes2-lane.txt 2>&1 &
"$RS" dev/ncgp-truth.R lane $L > $O/truth-lane.txt 2>&1 &
"$RS" dev/ncgp-equiv.R lane > $O/equiv-lane.txt 2>&1 &
"$RS" dev/ncgp-bigunit.R $L > $O/bigunit-lane.txt 2>&1 &
"$RS" dev/ncgp-nuggetrate.R > $O/nuggetrate.txt 2>&1 &
wait
"$RS" dev/ncgp-nuggetrate-sum.R > $O/nuggetrate-sum.txt 2>&1
echo STAGE-EVIDENCE

(
  export NOT_CRAN=true
  export FRMTMB_RELLIB=$L
  export TMP="${TMP:-$TEMP}"
  rm -rf dev/release/ncgpsess-sessions
  bash dev/release/run-sessions.sh ncgpsess 6 > $O/sess-driver.txt 2>&1
  cp dev/release/ncgpsess.log $O/ncgpsess.log 2>/dev/null
)
bash dev/ncgp-ob.sh > $O/ob.log 2>&1
rm -rf dev/ciharden-out
bash dev/ncgp-mut.sh > $O/mut.log 2>&1
echo STAGE-SESS

bash dev/ncgp-check.sh extensions/frmtmb.sample > $O/check-sample.txt 2>&1
bash dev/ncgp-check.sh . > $O/check-core.txt 2>&1
echo CHECK-DONE
