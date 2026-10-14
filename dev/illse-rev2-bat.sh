#!/usr/bin/env bash
# Reviewer of lane illse, round 2: rerun the lane's false-alarm battery and the
# 30-fit bitwise set (the lane's own scripts, read-only) under the
# reference BLAS and OpenBLAS 0.3.32 LAPACK, lane vs base.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-illse
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
L=dev/illse-rev2-log/bat
N=5
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { local c=$1 lg=$2; shift 2; local R="$RS"
  if [ $c = ob ]; then R="$OB"; export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
  else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  jobs_wait
  echo "$(date +%T) START $c $lg" >> $L/driver.txt
  "$R" "$@" < /dev/null > $L/$lg.txt 2>&1 & }
for c in ref ob; do
  for g in nl lin poly mixed bounded; do
    for a in lane base; do
      run $c falarm/$g-$c-$a dev/illse-sf-rev3-falarm.R $a $g $L/falarm/$g-$c-$a.tsv "1:20"
    done
  done
  for a in lane base; do
    run $c bitwise-$c-$a dev/illse-sf-rev4-bitwise.R $a $L/bitwise-$c-$a.rds
  done
done
wait
for c in ref ob; do
  "$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-$c-lane.rds $L/bitwise-$c-base.rds > $L/bitcmp-$c.txt 2>&1
done
ILLSE_CFGS="ref ob" "$RS" dev/illse-sf-falarm-cmp.R $L/falarm lane $L/falarm base > $L/falarm-cmp.txt 2>&1
echo "$(date +%T) all done" >> $L/driver.txt
