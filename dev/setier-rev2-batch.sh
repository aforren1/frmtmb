#!/usr/bin/env bash
# Reviewer of lane setier, re-check after punch round 1: rerun the
# round-1 review scripts on the current lane build into
# dev/setier-rev2-log/ (reference BLAS unless named -ob).
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-setier
cd "$W"
L=dev/setier-rev2-log
mkdir -p $L/mo
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/setier-out/Rob/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4
LANE=C:/Users/adf44/source/r/wt-setier-lib
bash dev/setier-rev2-attacks.sh > $L/attacks.driver 2>&1 &
for c in 1 2 3 4; do a=$(( (c-1)*50+1 )); b=$(( c*50 ))
  "$R" dev/setier-rev-mo.R $LANE "$a:$b" $L/mo/mo-lane-$c.tsv > $L/mo/mo-lane-$c.txt 2>&1 &
done
"$R" dev/setier-rev-probe.R lane > $L/probe-lane-ref.txt 2>&1 &
"$OB" dev/setier-rev-probe.R lane > $L/probe-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-c0k.R lane > $L/c0k-lane.txt 2>&1 &
"$R" dev/setier-rev-c0k-sweep.R lane > $L/c0ks-lane-ref.txt 2>&1 &
"$OB" dev/setier-rev-c0k-sweep.R lane > $L/c0ks-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-aexpb.R lane > $L/aexpb-lane-ref.txt 2>&1 &
"$OB" dev/setier-rev-aexpb.R lane > $L/aexpb-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-deg7.R lane 7 > $L/deg7-lane.txt 2>&1 &
"$R" dev/setier-rev-bnd.R lane > $L/bnd-lane.txt 2>&1 &
"$R" dev/setier-rev-sep2.R lane > $L/sep2-lane.txt 2>&1 &
"$R" dev/setier-rev-probit.R $LANE > $L/probit-lane.txt 2>&1 &
"$R" dev/setier-rev-bound.R lane > $L/bound-lane.txt 2>&1 &
"$R" dev/setier-rev-smooth2.R lane > $L/smooth2-lane.txt 2>&1 &
"$R" dev/setier-rev-condvar.R lane > $L/condvar-lane.txt 2>&1 &
"$R" dev/setier-rev-adapt.R lane > $L/adapt-lane.txt 2>&1 &
"$R" dev/setier-rev-miss.R lane > $L/miss-lane.txt 2>&1 &
"$R" dev/setier-rev-rs20.R $LANE > $L/rs20-lane.txt 2>&1 &
"$R" dev/setier-rev-fixes-table.R merge > $L/fixes-lane.txt 2>&1 &
"$R" dev/setier-singular.R $LANE $L/sing-lane-ref.tsv > $L/sing-lane-ref.txt 2>&1 &
"$OB" dev/setier-singular.R $LANE $L/sing-lane-ob.tsv > $L/sing-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-smallsd.R $LANE $L/smallsd-lane-ref.tsv 40 > $L/smallsd-lane-ref.txt 2>&1 &
"$R" dev/setier-sep.R $LANE $L/sep-lane-ref.tsv > $L/sep-lane-ref.txt 2>&1 &
wait
"$R" dev/setier-rev-sepcost.R lane > $L/sepcost-lane.txt 2>&1
echo done
