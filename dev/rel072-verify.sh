#!/usr/bin/env bash
# The verification runs of the 0.72.0 consolidation on the merged build
# (wt-release-lib), one after another, 8 R processes at a time:
# all eight suites one file per process under OpenBLAS 0.3.26 and
# 0.3.32 (BLAS and LAPACK, 4 threads) and the reference BLAS, every gate
# on (dev/ciharden-scan.sh); the escaped-warning scan
# (dev/rel072-warnscan.sh); the one-session runs, forward, reverse and
# preload (dev/release/run-sessions.sh, NOT_CRAN only, as coverage).
#   bash dev/rel072-verify.sh
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/3 TEMP=C:/Users/adf44/AppData/Local/Temp/3
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
LIB=C:/Users/adf44/source/r/wt-release-lib
cd "$ROOT"
V=dev/rel072-log/verify.txt
echo "$(date '+%F %T') start" >> $V
bash dev/ciharden-scan.sh r072-ob26 ob0.3.26-lapack-t4 $LIB 8 >> $V 2>&1
echo "$(date '+%F %T') ob26 done" >> $V
bash dev/ciharden-scan.sh r072-ob32 ob0.3.32-lapack-t4 $LIB 8 >> $V 2>&1
echo "$(date '+%F %T') ob32 done" >> $V
bash dev/ciharden-scan.sh r072b-ref ref $LIB 8 >> $V 2>&1
echo "$(date '+%F %T') ref done" >> $V
FRMTMB_RELLIB=$LIB bash dev/rel072-warnscan.sh dev/ciharden-log/r072b-ref.jobs 8 >> $V 2>&1
echo "$(date '+%F %T') warnscan done" >> $V
NOT_CRAN=true FRMTMB_RELLIB=$LIB bash dev/release/run-sessions.sh rel0.72.0 8 >> $V 2>&1
echo "$(date '+%F %T') sessions done" >> $V
