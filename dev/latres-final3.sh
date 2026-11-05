#!/usr/bin/env bash
# Fix round 2, on the final build, in order. Library paths are absolute in
# latres-env.sh; every test log's lib: line reads the lane's.
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
o=dev/latres-out
bash dev/latres-install.sh frmtmb.sample frmtmb.spatial > $o/instf.log 2>&1
tail -1 $o/instf.log
rm -f $o/f2core.out $o/f2gated.out
bash dev/latres-suite.sh f2core 7 core > $o/f2core.out 2>&1
echo phase-core
FRMTMB_BRMS_FIT_TESTS=true bash dev/latres-suite.sh f2gated 7 core frmtmb.sample \
  frmtmb.spatial frmtmb.survival frmtmb.latent frmtmb.eam frmtmb.learn \
  frmtmb.coupling frmtmb.ode > $o/f2gated.out 2>&1
echo phase-gated

rm -rf $o/r3
bash dev/latres-r3run.sh > $o/r3.out 2>&1
"$RSCRIPT" dev/latres-r3sum.R > $o/r3sum.txt 2>&1
"$RSCRIPT" dev/review-latres-unb.R > $o/unb.log 2>&1
"$RSCRIPT" dev/latres-r2.R > $o/r2.log 2>&1
"$RSCRIPT" dev/review-latres-m7c.R > $o/m7c.log 2>&1
"$RSCRIPT" dev/latres-r7.R > $o/r7.log 2>&1
"$RSCRIPT" dev/latres-m5.R > $o/m5.log 2>&1
echo phase-measure

rm -f $o/l1-w*.tsv $o/l1x-w*.tsv
for w in 1 2 3 4; do
  "$RSCRIPT" dev/latres-l1.R $w 4 $o/l1-w$w.tsv > $o/l1-w$w.log 2>&1 &
done
for w in 1 2 3; do
  "$RSCRIPT" dev/latres-l1x.R $w 3 $o/l1x-w$w.tsv > $o/l1x-w$w.log 2>&1 &
done
wait
"$RSCRIPT" dev/latres-report.R > $o/report-f2.txt 2>&1
echo phase-l1

rm -rf $o/grid-after; mkdir -p $o/grid-after
for w in 1 2 3 4 5 6 7; do
  "$RSCRIPT" dev/latres-grid.R $w 7 "$LATRES_LIB;$BASE_LIB" $o/grid-after \
    > $o/grid-after/w$w.log 2>&1 &
done
wait
"$RSCRIPT" dev/latres-grid-summary.R $o/grid-after > $o/grid-after-summary.txt 2>&1
"$RSCRIPT" dev/latres-gridcount.R > $o/gridcount.txt 2>&1
echo phase-grid

for p in p167 p156 p164 p154 p112 p028 p128 p132 p158; do
  rm -f dev/latres-equiv-out/$p.tsv
done
FRMTMB_EQUIV_LIB="$LATRES_LIB;$BASE_LIB" bash dev/equiv-run.sh 7 dev/latres-equiv-out \
  > $o/equiv-run.log 2>&1
echo phase-equiv

FRMTMB_RELLIB="$LATRES_LIB" NOT_CRAN=true bash dev/release/run-sessions.sh latres-sess 6 \
  > $o/sessions.log 2>&1
echo phase-sessions

for v in 0.3.26 0.3.32; do
  bash dev/ciharden-openblas.sh $v lapack > $o/ob-$v-build.log 2>&1
  R="dev/ciharden-out/Rob$v-lapack/bin/Rscript.exe"
  for f in latres latres-brms autocor autocor-cond fcor; do
    OPENBLAS_NUM_THREADS=4 FRMTMB_BRMS_FIT_TESTS=true "$R" dev/latres-test1.R \
      tests/testthat/test-$f.R > $o/ob-$v-$f.log 2>&1 &
    while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  done
  wait
  grep -h "^RESULT" $o/ob-$v-*.log
done
rm -rf dev/ciharden-out
echo phase-openblas

bash dev/latres-check.sh . frmtmb > $o/check-core.out 2>&1
bash dev/latres-check.sh extensions/frmtmb.sample frmtmb.sample > $o/check-sample.out 2>&1
echo final3-done
