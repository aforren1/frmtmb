#!/usr/bin/env bash
# Lane seflags, punch round 2c: the false-alarm battery, then the core,
# frmtmb.sample, frmtmb.spline and frmtmb.eam suites under both BLAS
# builds, one R process per file. Waits for every job.
#   bash dev/seflags-p2c-final.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
P=dev/seflags-log/p2c
LIST=$W/dev/ciharden-log/seflags-p2c.list
echo "start $(date)"
bash dev/seflags-p2c-falarm.sh
echo "falarm done $(date)"
bash dev/ciharden-scan.sh seflags-ref9 ref C:/Users/adf44/source/r/wt-seflags-lib 12 "$LIST" > $P/scan-ref9-driver.txt 2>&1 &
bash dev/ciharden-scan.sh seflags-ob9 ob0.3.26-lapack-t2 C:/Users/adf44/source/r/wt-seflags-lib 8 "$LIST" > $P/scan-ob9-driver.txt 2>&1 &
wait
echo "scans done $(date)"
