#!/usr/bin/env bash
# Fix round 2: every run the report rests on, in order, on the final
# build in wt-cdfs-lib, at most 8 R processes at once. The pre-fix arm
# (dev/release/cdfs-pre2) ran on its own.
#   bash dev/cdfs-fix2-runs.sh
cd "$(dirname "$0")/.."
L=C:/Users/adf44/source/r/wt-cdfs-lib
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
SP=C:/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
echo "start $(date +%H:%M)"
# 1. every file of the eight suites, one per process, NOT_CRAN only
rm -rf dev/release/cdfs-fx-all-files dev/release/cdfs-fx-all.log
( export NOT_CRAN=true FRMTMB_RELLIB=$L
  unset FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS FRMTMB_FUZZ FRMTMB_SCALE_TESTS
  bash dev/release/run-par.sh cdfs-fx-all dev/cdfs-fx-all.jobs 8 > /dev/null 2>&1 )
echo "ungated $(date +%H:%M): $(tail -1 dev/release/cdfs-fx-all.log)"
# 2. every file with every gate on, reference BLAS
bash dev/ciharden-scan.sh cdfs-fx-ref ref $L 8 > /dev/null 2>&1
echo "gated $(date +%H:%M): $(tail -1 dev/ciharden-log/cdfs-fx-ref.sum)"
# 3. the numerical files under OpenBLAS 0.3.26 and 0.3.32, lapack, 4 threads
for v in 0.3.26 0.3.32; do
  [ -x dev/ciharden-out/Rob$v-lapack/bin/x64/Rscript.exe ] ||
    bash dev/ciharden-openblas.sh $v lapack > dev/cdfs-out/ob-build-$v.txt 2>&1
done
bash dev/ciharden-scan.sh cdfs-fx-ob326 ob0.3.26-lapack-t4 $L 4 dev/cdfs-fx-ob.jobs > /dev/null 2>&1 &
bash dev/ciharden-scan.sh cdfs-fx-ob332 ob0.3.32-lapack-t4 $L 4 dev/cdfs-fx-ob.jobs > /dev/null 2>&1 &
wait
echo "openblas $(date +%H:%M)"
# 4. the one-session runs, NOT_CRAN only, beside R CMD check of core
( export NOT_CRAN=true FRMTMB_RELLIB=$L
  unset FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS FRMTMB_FUZZ FRMTMB_SCALE_TESTS
  bash dev/release/run-sessions.sh cdfs-fx-sess 6 > dev/cdfs-out/fx-sess-driver.txt 2>&1 ) &
bash dev/cdfs-check.sh > dev/cdfs-out/check-driver.txt 2>&1 &
wait
echo "sessions and check $(date +%H:%M)"
# 5. the measurements, the review's scripts, bin-2, the brms report and
#    the new tests on the pre-fix build
( for s in acc grad5 grad4; do "$R" $SP/rev/$s.R > dev/cdfs-out/rerev/$s.txt 2>&1; done ) &
( for s in zitrunc4 trlb fuzzweib exgdens2; do
    ARM=lane "$R" $SP/rev/$s.R > dev/cdfs-out/rerev/$s.txt 2>&1; done ) &
( "$R" dev/cdfs-fix1-igamma.R $L > dev/cdfs-out/fix1-igamma.txt 2>&1
  "$R" dev/cdfs-fix1-conv.R $L > dev/cdfs-out/fix1-conv.txt 2>&1
  "$R" dev/cdfs-fix1-exgfit.R $L > dev/cdfs-out/fix1-exgfit-lane.txt 2>&1
  "$R" dev/cdfs-fix1-tpois.R $L > dev/cdfs-out/fix1-tpois.txt 2>&1
  "$R" dev/cdfs-fix1-cost.R $L > dev/cdfs-out/fix1-cost.txt 2>&1
  "$R" dev/cdfs-fix1-imp.R $L > dev/cdfs-out/fix1-imp-wt-cdfs-lib.txt 2>&1
  "$R" dev/cdfs-fix2-simcost.R $L > dev/cdfs-out/fix2-simcost.txt 2>&1
  for s in sumtail sumfit compois win; do
    "$R" $SP/rev/$s.R > dev/cdfs-out/rerev/$s.txt 2>&1
  done
  "$R" dev/cdfs-fuzz-base.R $L 20322452 20446531 20502220 \
    > dev/cdfs-out/fuzz-base-mine.txt 2>&1 ) &
( export FRMTMB_PORT_LIB=$L
  for t in distributions log_lik; do
    "$R" dev/brmsport-bin2-run.R frmtmb \
      dev/cdfs-out/bin2/frmtmb/test-brms-suite-bin2-$t.R \
      dev/cdfs-out/bin2-rec/rec-frmtmb-$t.tsv \
      > dev/cdfs-out/bin2-rec/run-frmtmb-$t.txt 2>&1
  done
  "$R" dev/cdfs-bin2-flips.R > dev/cdfs-out/bin2-flips.txt 2>&1
  "$R" dev/cdfs-brms-lik-report.R > dev/cdfs-out/brms-lik-report.txt 2>&1 ) &
wait
echo "done $(date +%H:%M)"
