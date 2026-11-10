#!/usr/bin/env bash
# OpenBLAS with LAPACK on the merged build, as the ubuntu runners round:
# 0.3.26 on test-se-tier.R, test-se-check.R, test-families-extreme.R and
# test-cens-trunc.R; 0.3.32 on test-se-tier.R (the brief's list). 4
# threads, 3 files at once; every gate on. The R copies are deleted after.
#   bash dev/rel081-ob.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r20
L=dev/rel081-log
bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build.txt 2>&1; echo "build 0.3.26 exit $?"
bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob32-build.txt 2>&1; echo "build 0.3.32 exit $?"
bash dev/ciharden-scan.sh r081-ob26 ob0.3.26-lapack-t4 base 3 $L/ob26.jobs
bash dev/ciharden-scan.sh r081-ob32 ob0.3.32-lapack-t4 base 3 $L/ob32.jobs
cat dev/ciharden-log/r081-ob26.sum dev/ciharden-log/r081-ob32.sum
grep -h '^blas:' dev/ciharden-log/r081-ob26-suite/*.txt dev/ciharden-log/r081-ob32-suite/*.txt
rm -rf dev/ciharden-out
echo "OB END"
