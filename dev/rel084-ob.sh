#!/usr/bin/env bash
# OpenBLAS with LAPACK on the merged build, as the ubuntu runners round:
# 0.3.26 and 0.3.32, each with OPENBLAS_CORETYPE=Haswell and unset, on the
# 23 files of dev/rel084-log/ob.jobs (the brief's list). 4 threads, OB_N
# files at once (3 by default); every gate on. The R copies are built
# first, alone, and deleted after. Adapted from dev/rel083-ob.sh.
#   bash dev/rel084-ob.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel084 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r23
L=dev/rel084-log
bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build.txt 2>&1; echo "build 0.3.26 exit $?"
bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob32-build.txt 2>&1; echo "build 0.3.32 exit $?"
for v in 26 32; do
  unset OPENBLAS_CORETYPE
  bash dev/ciharden-scan.sh r084-ob$v ob0.3.$v-lapack-t4 base ${OB_N:-3} $L/ob.jobs
  OPENBLAS_CORETYPE=Haswell bash dev/ciharden-scan.sh r084-ob$v-hsw ob0.3.$v-lapack-t4 base ${OB_N:-3} $L/ob.jobs
done
for r in ob26 ob26-hsw ob32 ob32-hsw; do
  echo "== $r"; cat dev/ciharden-log/r084-$r.sum
  grep -h '^blas:' dev/ciharden-log/r084-$r-suite/*.txt | sort | uniq -c | head -3
  grep -h '^lib:' dev/ciharden-log/r084-$r-suite/*.txt | awk '{print $7}' | sort | uniq -c
done
rm -rf dev/ciharden-out
echo "OB END"
