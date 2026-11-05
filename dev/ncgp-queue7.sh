#!/usr/bin/env bash
# Fix round 2's final evidence on the frozen build. Library paths are absolute in each
# script; read `lib:` in each log. R CMD check runs after this, from dev/ncgp-queue8.sh.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-ncgp-lib
O=dev/ncgp-out
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
mkdir -p $O/r2 $O/m1 $O/rs

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

# F1, F2, F3, the review's rows and fresh seeds
"$RS" dev/review-ncgp-r2-edge2.R lane > $O/r2/edge2.txt 2>&1 &
"$RS" dev/review-ncgp-r2-ng.R lane "1" > $O/r2/ng.txt 2>&1 &
"$RS" dev/review-ncgp-r2-fresh.R lane "1:2" > $O/r2/fresh-review.txt 2>&1 &
"$RS" dev/ncgp-r2fresh.R lane "11:12" > $O/r2/fresh-new.txt 2>&1 &
"$RS" dev/ncgp-r2fc.R "1:3" > $O/r2/fc.txt 2>&1 &
"$RS" dev/review-ncgp-modes2.R lane > $O/modes2-lane.txt 2>&1 &
wait
"$RS" dev/ncgp-truth.R lane $L > $O/truth-lane.txt 2>&1 &
"$RS" dev/ncgp-equiv.R lane > $O/equiv-lane.txt 2>&1 &
"$RS" dev/ncgp-bigunit.R $L > $O/bigunit-lane.txt 2>&1 &
"$RS" dev/ncgp-bigunit.R base > $O/bigunit-base.txt 2>&1 &
"$RS" dev/ncgp-nuggetrate.R > $O/nuggetrate.txt 2>&1 &
"$RS" dev/ncgp-fcwarn.R > $O/fcwarn.txt 2>&1 &
"$RS" dev/ncgp-keep.R > $O/keep.txt 2>&1 &
wait
"$RS" dev/ncgp-nuggetrate-sum.R > $O/nuggetrate-sum.txt 2>&1
echo STAGE-ROWS

# the noise floors, noise-free, gpgap's reviewer's floor
"$RS" dev/review-ncgp-floor.R lane 40 > $O/m1/floor-40.txt 2>&1 &
"$RS" dev/review-ncgp-floor.R lane 100 > $O/m1/floor-100.txt 2>&1 &
"$RS" dev/review-ncgp-nfpred.R > $O/m1/nfpred.txt 2>&1 &
"$RS" dev/review-ncgp-mpfr.R > $O/m1/mpfr.txt 2>&1 &
"$RS" dev/review-ncgp-lanerows.R > $O/m1/lanerows.txt 2>&1 &
"$RS" dev/review-ncgp-stationary.R > $O/m1/stationary.txt 2>&1 &
"$RS" ../frmtmb-wt-gpgap/dev/review-gpgap-ncgpfloor.R > $O/m1/ncgpfloor.txt 2>&1 &
wait
"$RS" dev/ncgp-floorsum.R > $O/floorsum.txt 2>&1
echo STAGE-FLOORS

bash dev/ncgp-l1-rest.sh > $O/l1-rest.log 2>&1
bash dev/ncgp-l1eq.sh > $O/l1eq.log 2>&1
bash dev/ncgp-l1ng.sh > $O/l1ng.log 2>&1
"$RS" dev/ncgp-row10.R > $O/row10.txt 2>&1 &
"$RS" dev/ncgp-readers.R > $O/readers.txt 2>&1 &
"$RS" dev/ncgp-readers.R base > $O/readers-base.txt 2>&1 &
"$RS" dev/ncgp-grcheck.R > $O/grcheck.txt 2>&1 &
"$RS" dev/ncgp-prior.R > $O/prior.txt 2>&1 &
"$RS" dev/ncgp-cholscan.R > $O/cholscan.txt 2>&1 &
wait
echo STAGE-L1

# the start rule, memory, the sampler's read: one chain set at a time, own port each
"$RS" dev/ncgp-randstart.R 1e-6 default > $O/rs/1e-6-default.txt 2>&1 &
"$RS" dev/ncgp-randstart.R 1e-6 mode > $O/rs/1e-6-mode.txt 2>&1 &
"$RS" dev/ncgp-randstart.R 1e-3 default > $O/rs/1e-3-default.txt 2>&1 &
wait
R_PARALLEL_PORT=11841 "$RS" dev/review-ncgp-r2-start.R 7701 > $O/r2/start.txt 2>&1
R_PARALLEL_PORT=11843 "$RS" dev/review-ncgp-mem.R lane 300 2 > $O/mem-lane-300.txt 2>&1
R_PARALLEL_PORT=11847 "$RS" dev/ncgp-sampler.R 1e-2 30 1000 0.99 > $O/post/sampler-1e-2.txt 2>&1
echo STAGE-SAMPLER

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
echo QUEUE7-DONE
