#!/usr/bin/env bash
# Lane seflags, punch round 1: the review's cost designs
# (dev/seflags-rev-cost2.R) and the 2051-coefficient bernoulli fit
# (dev/seflags-bigp.R), lane and base side by side per design, into
# dev/seflags-log/cost-p1/.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/seflags-log/cost-p1
mkdir -p $L
for dz in glmm_1e5_year glmm_1e5 lmm_400 glm_400; do
  for a in lane base; do
    "$RS" dev/seflags-rev-cost2.R $a $dz 3 > $L/$dz-$a.txt 2>&1 &
  done
  wait
done
for a in lane base; do
  "$RS" dev/seflags-bigp.R $a > $L/bigp-$a.txt 2>&1 &
done
wait
echo COST DONE
