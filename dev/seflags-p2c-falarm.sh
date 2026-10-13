#!/usr/bin/env bash
# Lane seflags, punch round 2c: the reviewer's false-alarm battery
# (dev/seflags-rev3-falarm.R) on the lane build, both BLAS builds, one
# process per group, into dev/seflags-log/p2c/falarm. The base arm is
# the reviewer's own run of rellib-r7 (dev/seflags-rev3-log/falarm),
# copied beside it, since the base did not change. Waits for every job.
#   bash dev/seflags-p2c-falarm.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
L=dev/seflags-log/p2c/falarm
mkdir -p $L
cp dev/seflags-rev3-log/falarm/*-base.tsv $L/
for g in nl lin poly mixed bounded; do
  "$RS" dev/seflags-rev3-falarm.R lane $g $L/$g-ref-lane.tsv "1:20" > $L/$g-ref-lane.txt 2>&1 &
  "$OB" dev/seflags-rev3-falarm.R lane $g $L/$g-ob-lane.tsv "1:20" > $L/$g-ob-lane.txt 2>&1 &
done
wait
"$RS" dev/seflags-p2c-falarm-sum.R > $L/../falarm-sum.txt 2>&1
