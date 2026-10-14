#!/usr/bin/env bash
# Lane illse: the SE-check batteries of the 0.70.0 consolidation
# (dev/rel070-seflags-run.sh) on the lane library ("lane") and on
# rellib-r8 ("base"). TAG names the log directory, CFGS is "ref", "ob"
# or "ref ob", OBV the OpenBLAS release of the ob arm.
#   TAG=r1 CFGS=ref bash dev/illse-sf-run.sh
# At most 6 R processes at once. Waits for every job.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-illse
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob${OBV:-0.3.26}-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-illse/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
L=dev/illse-log/sf-$TAG
mkdir -p $L/falarm
N=6
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { # <cfg> <logname> <script> args...
  local c=$1 lg=$2; shift 2
  local R="$RS"
  if [ $c = ob ]; then R="$OB"; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  jobs_wait
  echo "$(date +%T) START $c $lg :: $*" >> $L/driver.txt
  "$R" "$@" < /dev/null > $L/$lg.txt 2>&1 &
}
for c in $CFGS; do
  for g in nl lin poly mixed bounded; do
    run $c falarm/$g-$c-lane dev/illse-sf-rev3-falarm.R lane $g $L/falarm/$g-$c-lane.tsv "1:20"
    run $c falarm/$g-$c-base dev/illse-sf-rev3-falarm.R base $g $L/falarm/$g-$c-base.tsv "1:20"
  done
  for a in lane base; do
    t=lane; [ $a = base ] && t=base
    for s in whole cum illpred routes; do
      run $c $s-$c-$t dev/illse-sf-rev4-$s.R $a
    done
    run $c illid-$c-$t dev/illse-sf-rev4-illid.R $a $L/illid-$c-$t.tsv
    run $c glm-$c-$t dev/illse-sf-rev4-glm.R $a $L/glm-$c-$t.tsv
    run $c bitwise-$c-$t dev/illse-sf-rev4-bitwise.R $a $L/bitwise-$c-$t.rds
    run $c deg7-$c-$t dev/illse-sf-rev2-deg7.R $a
    run $c deg7exact-$c-$t dev/illse-sf-rev3-deg7exact.R $a
    run $c deg7pred-$c-$t dev/illse-sf-rev3-deg7pred.R $a
  done
  for a in lane base; do
    run $c gain-$c-$a dev/illse-sf-rev4-gain.R $a $L/gain-$c-$a.tsv
    run $c join-$c-$a dev/illse-sf-rev2-join.R $a
    run $c sd2-$c-$a dev/illse-sf-rev2-sd2.R $a
    run $c pred3x-$c-$a dev/illse-sf-rev3-pred3x.R $a
  done
  run $c b1ref-$c-lane dev/illse-sf-b1ref.R
done
wait
for c in $CFGS; do
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-lane.rds $L/bitwise-$c-base.rds > $L/bitcmp-$c-lane-vs-base.txt 2>&1
done
echo "$(date +%T) all done" >> $L/driver.txt
