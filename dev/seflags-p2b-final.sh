#!/usr/bin/env bash
# Lane seflags, punch round 2b: every rerun on the final build, in turn,
# then the two scans of all eight suites. Waits for every job.
#   bash dev/seflags-p2b-final.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
P=dev/seflags-log/p2b
echo "start $(date)"
rm -rf $P/run-lane $P/run-m9
bash dev/seflags-p2b-run.sh lane $P/run-lane &
bash dev/seflags-p2b-run.sh m9 $P/run-m9 &
wait
echo "p2b runs done $(date)"
bash dev/seflags-p2b-batch.sh
echo "batch done $(date)"
bash dev/ciharden-scan.sh seflags-ref8 ref C:/Users/adf44/source/r/wt-seflags-lib 12 > $P/scan-ref8-driver.txt 2>&1 &
bash dev/ciharden-scan.sh seflags-ob8 ob0.3.26-lapack-t2 C:/Users/adf44/source/r/wt-seflags-lib 8 > $P/scan-ob8-driver.txt 2>&1 &
wait
echo "scans done $(date)"
