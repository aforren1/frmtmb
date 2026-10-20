#!/usr/bin/env bash
# Lane illse2, fix round 1 resumed: every accuracy script of
# dev/illse2-f2-jobs.sh again, on the final install, with every log in
# dev/illse2-log/f4/ so that the install-20 logs stay for comparison.
# Adds the seen-to-fail runs of dev/illse2-fx-mock.R.
#   bash dev/illse2-f4-jobs.sh [P]
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-illse2
P="${1:-6}"
R=/c/PROGRA~1/R/R-46~1.1/bin/Rscript.exe
L=dev/illse2-log/f4
mkdir -p $L
J=$L/jobs.txt
: > $J
add() { echo "ILLSE2_ARM=lane $R $1 > $L/$2 2>&1" >> $J; }
for c in "0 5" "0 6" "0 7" "0 8" "1 5" "1 6" "1 7" "1 8"; do
  set -- $c
  add "dev/illse2-p1.R $1 $2 1:3" "p1-lane-$1-$2.txt"
done
add "dev/illse2-p2.R 7 4" "p2-lane-7-4.txt"
for dg in 5 6 7 8; do
  add "dev/illse2-edge4.R $dg" "edge4-ref-$dg.txt"
  add "dev/illse2-p3reml.R $dg 1:3" "p3reml-lane-$dg.txt"
done
for fam in bernoulli poisson; do
  add "dev/illse2-p3remlglm.R $fam 7 1:3" "p3remlglm-lane-$fam.txt"
done
for c in poly sin rat; do
  add "dev/illse2-p3nl.R $c 1:3" "p3nl-lane-$c.txt"
done
for job in "gaussian 0" "gaussian 1" "bernoulli 0" "poisson 0"; do
  set -- $job
  add "dev/illse2-acc.R $1 $2 6:8 1:3" "acc-$1-$2.txt"
done
add "dev/illse2-p4bound.R 5:8 1:5" "p4bound-lane.txt"
for K in 150 292; do
  add "dev/illse2-p4wide.R $K 2000 1:2" "p4wide-$K.txt"
done
add "dev/illse2-asym2.R" "asym2.txt"
add "dev/illse2-fx-remlll.R" "remlll-lane.txt"
for A in 1 0.2 0.05; do
  add "dev/illse2-rev-asym.R $A 1:2" "asym-$A-lane.txt"
done
add "dev/illse2-rev-esc2.R" "esc2-lane.txt"
add "dev/illse2-p2tier.R" "p2tier-lane.txt"
add "dev/illse2-p3ratdbg.R 3" "p3ratdbg-lane-3.txt"
for h in "poisson 6 1e5" "bernoulli 7 1e5" "poisson 8 1e4"; do
  set -- $h
  add "dev/illse2-rev-hess.R $1 $2 $3 1:2" "hess-$1$2-lane.txt"
done
for m in esc17 schemeE noguard allguard; do
  add "dev/illse2-fx-mock.R $m" "mock-$m.txt"
done
add "dev/illse2-fx-nojac.R" "t-reml-nojac.txt"
wc -l $J
xargs -P "$P" -I{} bash -c "{}" < $J
echo "ALL DONE" >> $L/done.txt
