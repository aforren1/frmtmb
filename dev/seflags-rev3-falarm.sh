#!/usr/bin/env bash
# Reviewer of lane seflags (final check): the false-alarm battery, lane and base, both
# BLAS builds, one process per group.
#   bash dev/seflags-rev3-falarm.sh <seeds expr>
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
L=dev/seflags-rev3-log/falarm
mkdir -p $L
for g in nl lin poly mixed bounded; do
  for a in lane base; do
    "$RS" dev/seflags-rev3-falarm.R $a $g $L/$g-ref-$a.tsv "$1" > $L/$g-ref-$a.txt 2>&1 &
    "$OB" dev/seflags-rev3-falarm.R $a $g $L/$g-ob-$a.tsv "$1" > $L/$g-ob-$a.txt 2>&1 &
  done
done
wait
echo all done
