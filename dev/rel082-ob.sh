#!/usr/bin/env bash
# OpenBLAS with LAPACK on the merged build, as the ubuntu runners round:
# 0.3.26 on the eleven files of dev/rel082-log/ob26.jobs and 0.3.32 on the three
# of ob32.jobs (the brief's list). 4 threads, OB_N files at once (3 by
# default); every gate on. The R copies are deleted after.
#   bash dev/rel082-ob.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/rel082-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r21
L=dev/rel082-log
bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build.txt 2>&1; echo "build 0.3.26 exit $?"
bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob32-build.txt 2>&1; echo "build 0.3.32 exit $?"
bash dev/ciharden-scan.sh r082-ob26 ob0.3.26-lapack-t4 base ${OB_N:-3} $L/ob26.jobs
bash dev/ciharden-scan.sh r082-ob32 ob0.3.32-lapack-t4 base ${OB_N:-3} $L/ob32.jobs
cat dev/ciharden-log/r082-ob26.sum dev/ciharden-log/r082-ob32.sum
grep -h '^blas:' dev/ciharden-log/r082-ob26-suite/*.txt dev/ciharden-log/r082-ob32-suite/*.txt
rm -rf dev/ciharden-out
echo "OB END"
