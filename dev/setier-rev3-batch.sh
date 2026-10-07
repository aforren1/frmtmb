#!/usr/bin/env bash
# Reviewer of lane setier, final check after punch round 2: rerun the
# review scripts on the current lane build into dev/setier-rev3-log/.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-setier
cd "$W"
L=dev/setier-rev3-log
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/setier-out/Rob/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4
LANE=C:/Users/adf44/source/r/wt-setier-lib
"$R" dev/setier-rev2-window.R lane > $L/window-lane.txt 2>&1 &
"$OB" dev/setier-rev2-window.R lane > $L/window-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-probe.R lane poly curved c0k expb bound > $L/probe-lane.txt 2>&1 &
"$R" dev/setier-rev-c0k.R lane > $L/c0k-lane.txt 2>&1 &
"$R" dev/setier-rev-c0k-sweep.R lane > $L/c0ks-lane.txt 2>&1 &
"$R" dev/setier-rev-aexpb.R lane > $L/aexpb-lane.txt 2>&1 &
"$R" dev/setier-rev-deg7.R lane 7 > $L/deg7-lane.txt 2>&1 &
"$R" dev/setier-rev2-trap.R $LANE > $L/trap-lane.txt 2>&1 &
"$R" dev/setier-rev2-short.R $LANE > $L/short-lane.txt 2>&1 &
"$R" dev/setier-rev2-scale2.R $LANE > $L/scale2-lane.txt 2>&1 &
"$R" dev/setier-rev2-item8.R $LANE > $L/item8-lane.txt 2>&1 &
"$R" dev/setier-rev2-sepmiss.R $LANE budget few > $L/sepmiss-lane.txt 2>&1 &
"$R" dev/setier-rev2-sepmiss.R $LANE bigp > $L/sepmiss-bigp-lane.txt 2>&1 &
"$R" dev/setier-rev-rs20.R $LANE > $L/rs20-lane.txt 2>&1 &
"$R" dev/setier-singular.R $LANE $L/sing-lane-ref.tsv > $L/sing-lane-ref.txt 2>&1 &
"$OB" dev/setier-singular.R $LANE $L/sing-lane-ob.tsv > $L/sing-lane-ob.txt 2>&1 &
"$R" dev/setier-rev-smallsd.R $LANE $L/smallsd-lane.tsv 40 > $L/smallsd-lane.txt 2>&1 &
"$R" dev/setier-sep.R $LANE $L/sep-lane-ref.tsv > $L/sep-lane-ref.txt 2>&1 &
"$R" dev/setier-rev2-bnd.R lane > $L/bnd-lane.txt 2>&1 &
"$R" dev/setier-rev-smooth2.R lane > $L/smooth2-lane.txt 2>&1 &
"$R" dev/setier-rev-fixes-table.R merge > $L/fixes-lane.txt 2>&1 &
wait
"$R" dev/setier-rev-sepcost.R lane > $L/sepcost-lane.txt 2>&1
echo done
