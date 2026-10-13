#!/usr/bin/env bash
# rel070 seflags verification driver: the 2c check's scripts on rellib-r8
# ("rel", via dev/rel070-seflags-* copies whose lane arm is rellib-r8) and
# on rellib-r7 ("base", 0.69.1), reference BLAS and OpenBLAS 0.3.26 LAPACK.
# At most 6 R processes at once. Waits for every job.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
L=dev/rel070-log/studies/seflags
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
for c in ref ob; do
  for g in nl lin poly mixed bounded; do
    run $c falarm/$g-$c-rel dev/rel070-seflags-rev3-falarm.R lane $g $L/falarm/$g-$c-rel.tsv "1:20"
    run $c falarm/$g-$c-base dev/rel070-seflags-rev3-falarm.R base $g $L/falarm/$g-$c-base.tsv "1:20"
  done
  for a in lane base; do
    t=rel; [ $a = base ] && t=base
    for s in whole cum illpred routes; do
      run $c $s-$c-$t dev/rel070-seflags-rev4-$s.R $a
    done
    run $c illid-$c-$t dev/rel070-seflags-rev4-illid.R $a $L/illid-$c-$t.tsv
    run $c glm-$c-$t dev/rel070-seflags-rev4-glm.R $a $L/glm-$c-$t.tsv
    run $c bitwise-$c-$t dev/rel070-seflags-rev4-bitwise.R $a $L/bitwise-$c-$t.rds
    run $c deg7-$c-$t dev/rel070-seflags-rev2-deg7.R $a
    run $c deg7exact-$c-$t dev/rel070-seflags-rev3-deg7exact.R $a
    run $c deg7pred-$c-$t dev/rel070-seflags-rev3-deg7pred.R $a
  done
  run $c gain-$c-rel dev/rel070-seflags-rev4-gain.R lane $L/gain-$c-rel.tsv
  run $c join-$c-rel dev/rel070-seflags-rev2-join.R lane
  run $c sd2-$c-rel dev/rel070-seflags-rev2-sd2.R lane
  run $c pred3x-$c-rel dev/rel070-seflags-rev3-pred3x.R lane
  run $c b1ref-$c-rel dev/rel070-seflags-b1ref.R
done
wait
for c in ref ob; do
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-rel.rds $L/bitwise-$c-base.rds > $L/bitcmp-$c-rel-vs-base.txt 2>&1
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-rel.rds C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-log/bitwise-$c-base.rds > $L/bitcmp-$c-rel-vs-rev4base.txt 2>&1
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-base.rds C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-log/bitwise-$c-base.rds > $L/bitcmp-$c-base-vs-rev4base.txt 2>&1
done
echo "$(date +%T) all done" >> $L/driver.txt
