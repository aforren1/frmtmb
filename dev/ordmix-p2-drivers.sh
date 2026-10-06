#!/usr/bin/env bash
# Punch round 2: every warning driver on the final lane build, in
# turn. The review's margin script (dev/ordmix-rev2-margin.R, its ten
# designs) to dev/ordmix-p2-margin-log/, its B1 edge script
# (dev/ordmix-rev2-b1.R, its seven cases) to dev/ordmix-p2-b1-log/, the
# round-1 false-alarm driver (dev/ordmix-p1-fa.sh) to
# dev/ordmix-p1-fa-log/ (round 1's copy: dev/ordmix-p1-fa-log-round1/),
# and the criterion measurements (dev/ordmix-p2-crit.sh).
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$ROOT/dev/ordmix-p2-margin-log" "$ROOT/dev/ordmix-p2-b1-log"
for c in "cauchit 1 500" "cauchit 5 500" "cauchit 10 500" "cloglog 1 500" \
         "cloglog 5 500" "logit 1 5000" "logit 5 500" "logit 10 500" \
         "logit 10 2000" "probit 5 500"; do
  set -- $c
  ( "$RS" "$ROOT/dev/ordmix-rev2-margin.R" $1 $2 $3 \
      > "$ROOT/dev/ordmix-p2-margin-log/$1-$2-$3.txt" 2>&1 ) &
done
for k in hu1_only thres_only beta11 none nopred_mu1 disc_differs disc_same; do
  ( "$RS" "$ROOT/dev/ordmix-rev2-b1.R" $k 10 \
      > "$ROOT/dev/ordmix-p2-b1-log/$k.txt" 2>&1 ) &
done
wait
bash "$ROOT/dev/ordmix-p1-fa.sh"
bash "$ROOT/dev/ordmix-p2-crit.sh"
echo DRIVERS DONE
