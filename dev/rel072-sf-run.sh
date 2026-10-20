#!/usr/bin/env bash
# Release 0.72.0 (copy of dev/adhess-sf-run.sh): the SE false-alarm battery (5 model groups, seeds 1:20)
# and the bitwise set of lane illse (dev/illse-sf-run.sh), on the lane
# library ("lane", here wt-release-lib) and on rellib-r9 ("base"). TAG names the log
# directory, CFGS is "ref", "ob" or "ref ob", OBV the OpenBLAS release
# of the ob arm. At most N R processes at once.
#   TAG=a1 CFGS="ref ob" OBV=0.3.32 bash dev/rel072-sf-run.sh
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob${OBV:-0.3.32}-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
L=dev/rel072-log/sf-$TAG
mkdir -p $L/falarm
N=${N:-8}
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
    for a in lane base; do
      run $c falarm/$g-$c-$a dev/rel072-sf-rev3-falarm.R $a $g $L/falarm/$g-$c-$a.tsv "1:20"
    done
  done
  for a in lane base; do
    run $c bitwise-$c-$a dev/rel072-sf-rev4-bitwise.R $a $L/bitwise-$c-$a.rds
  done
done
wait
for c in $CFGS; do
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-lane.rds $L/bitwise-$c-base.rds > $L/bitcmp-$c-lane-vs-base.txt 2>&1
done
ILLSE_CFGS="$CFGS" "$RS" dev/rel072-sf-falarm-cmp.R $L/falarm lane $L/falarm base > $L/falarm-cmp.txt 2>&1
echo "$(date +%T) all done" >> $L/driver.txt
