#!/usr/bin/env bash
# Everything again on the final build, after the factor moved to a farthest-first
# order and the sampler's columns went back to the template's order. Absolute library
# paths are in each script; read `lib:` in each log.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-ncgp-lib
O=dev/ncgp-out
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"

cp $O/truth-lane.txt $O/truth-lane-sorted.txt
"$RS" dev/ncgp-truth.R lane $L > $O/truth-lane.txt 2>&1 &
"$RS" dev/ncgp-cholscan.R > $O/cholscan.txt 2>&1 &
"$RS" dev/ncgp-krig.R > $O/krig.txt 2>&1 &
"$RS" dev/ncgp-eps.R > $O/eps.txt 2>&1 &
"$RS" dev/ncgp-grcheck.R > $O/grcheck.txt 2>&1 &
"$RS" dev/ncgp-readers.R > $O/readers.txt 2>&1 &
wait
"$RS" dev/ncgp-row10.R > $O/row10.txt 2>&1 &
"$RS" dev/ncgp-prior.R > $O/prior.txt 2>&1 &
"$RS" dev/ncgp-noisefree2.R > $O/noisefree2.txt 2>&1 &
wait
echo STAGE1-DONE

bash dev/ncgp-l1-rest.sh > $O/l1-rest.log 2>&1
: > $O/l1eq.txt
for k in p022 p023 p094 p109 p136; do
  "$RS" dev/ncgp-l1eq.R $k 2>/dev/null | grep -a NCGPE >> $O/l1eq.txt
done
echo STAGE2-DONE

rm -rf $O/gated
FRMTMB_BRMS_FIT_TESTS=true bash dev/ncgp-par.sh $O/gated frmtmb frmtmb.sample frmtmb.spatial frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.survival
grep -h '^RESULT' $O/gated/*.txt > $O/gated.txt
(
  export NOT_CRAN=true
  NCGP_TLIB="C:/Users/adf44/source/r/rellib-r16" "$RS" dev/ncgp-runtest.R frmtmb tests/testthat/test-gp-noncentered.R > $O/seen-fail.txt 2>&1
)
echo STAGE3-DONE

(
  export NOT_CRAN=true
  export FRMTMB_RELLIB=$L
  export TMP="${TMP:-$TEMP}"
  rm -rf dev/release/ncgpsess-sessions dev/release/ncgpsess2-sessions
  bash dev/release/run-sessions.sh ncgpsess 6 > $O/sess-driver.txt 2>&1
  cp dev/release/ncgpsess.log $O/ncgpsess.log 2>/dev/null
)
bash dev/ncgp-ob.sh > $O/ob.log 2>&1
rm -rf dev/ciharden-out
echo STAGE4-DONE

bash dev/ncgp-mut.sh > $O/mut.log 2>&1
"$RS" dev/ncgp-sampler.R 1e-2 30 1000 0.99 > $O/post/sampler-1e-2.txt 2>&1
echo STAGE5-DONE

bash dev/ncgp-perf.sh > $O/perf.log 2>&1
echo STAGE6-DONE

bash dev/ncgp-check.sh extensions/frmtmb.sample > $O/check-sample.txt 2>&1
bash dev/ncgp-check.sh . > $O/check-core.txt 2>&1
echo CHECK-DONE
