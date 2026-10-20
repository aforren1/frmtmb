#!/usr/bin/env bash
# Lane illse2, fix round 1: every accuracy script on the final build,
# P at a time. Writes the job list and runs it.
#   bash dev/illse2-f2-jobs.sh [P]
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-illse2
P="${1:-6}"
R=/c/PROGRA~1/R/R-46~1.1/bin/Rscript.exe
L=dev/illse2-log/f2
mkdir -p $L dev/illse2-log/fx
J=$L/jobs1.txt
: > $J
for c in "0 5" "0 6" "0 7" "0 8" "1 5" "1 6" "1 7" "1 8"; do
  set -- $c
  echo "ILLSE2_ARM=lane $R dev/illse2-p1.R $1 $2 1:3 > dev/illse2-log/p1-lane-$1-$2.txt 2>&1" >> $J
done
echo "ILLSE2_ARM=lane $R dev/illse2-p2.R 7 4 > dev/illse2-log/p2-lane-7-4.txt 2>&1" >> $J
for dg in 5 6 7 8; do
  echo "ILLSE2_ARM=lane $R dev/illse2-edge4.R $dg > $L/edge4-ref-$dg.txt 2>&1" >> $J
  echo "ILLSE2_ARM=lane $R dev/illse2-p3reml.R $dg 1:3 > dev/illse2-log/p3reml-lane-$dg.txt 2>&1" >> $J
done
for fam in bernoulli poisson; do
  echo "ILLSE2_ARM=lane $R dev/illse2-p3remlglm.R $fam 7 1:3 > dev/illse2-log/p3remlglm-lane-$fam.txt 2>&1" >> $J
done
for c in poly sin rat; do
  echo "ILLSE2_ARM=lane $R dev/illse2-p3nl.R $c 1:3 > dev/illse2-log/p3nl-lane-$c.txt 2>&1" >> $J
done
for job in "gaussian 0" "gaussian 1" "bernoulli 0" "poisson 0"; do
  set -- $job
  echo "ILLSE2_ARM=lane $R dev/illse2-acc.R $1 $2 6:8 1:3 > $L/acc-$1-$2.txt 2>&1" >> $J
done
echo "ILLSE2_ARM=lane $R dev/illse2-p4bound.R 5:8 1:5 > dev/illse2-log/p4bound-lane.txt 2>&1" >> $J
for K in 150 292; do
  echo "ILLSE2_ARM=lane $R dev/illse2-p4wide.R $K 2000 1:2 > dev/illse2-log/p4wide-$K-final.txt 2>&1" >> $J
done
echo "ILLSE2_ARM=lane $R dev/illse2-asym2.R > dev/illse2-log/asym2.txt 2>&1" >> $J
echo "ILLSE2_ARM=lane $R dev/illse2-fx-remlll.R > dev/illse2-log/fx/remlll-lane.txt 2>&1" >> $J
for A in 1 0.2 0.05; do
  echo "ILLSE2_ARM=lane $R dev/illse2-rev-asym.R $A 1:2 > dev/illse2-log/fx/asym-$A-lane.txt 2>&1" >> $J
done
echo "ILLSE2_ARM=lane $R dev/illse2-rev-esc2.R > dev/illse2-log/fx/esc2-lane.txt 2>&1" >> $J
echo "ILLSE2_ARM=lane $R dev/illse2-p2tier.R > dev/illse2-log/p2tier-lane.txt 2>&1" >> $J
echo "ILLSE2_ARM=lane $R dev/illse2-p3ratdbg.R 3 > dev/illse2-log/p3ratdbg-lane-3.txt 2>&1" >> $J
for h in "poisson 6 1e5" "bernoulli 7 1e5" "poisson 8 1e4"; do
  set -- $h
  echo "ILLSE2_ARM=lane $R dev/illse2-rev-hess.R $1 $2 $3 1:2 > dev/illse2-log/fx/hess-$1$2-lane.txt 2>&1" >> $J
  echo "ILLSE2_ARM=base $R dev/illse2-rev-hess.R $1 $2 $3 1:2 > dev/illse2-log/fx/hess-$1$2-base.txt 2>&1" >> $J
done
wc -l $J
xargs -P "$P" -I{} bash -c "{}" < $J
echo "ALL DONE" >> $L/done.txt
