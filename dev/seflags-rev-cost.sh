#!/usr/bin/env bash
# Reviewer of lane seflags: dev/seflags-rev-cost.R, lane and base of each
# design run side by side (the same load on both arms).
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/seflags-rev-log/cost
mkdir -p $L
for dz in "$@"; do
  for a in lane base; do
    "$RS" dev/seflags-rev-cost.R $a $dz 3 > $L/$dz-$a.txt 2>&1 &
  done
  wait
done
echo done
