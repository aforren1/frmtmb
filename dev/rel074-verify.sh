#!/usr/bin/env bash
# The verification runs of the 0.74.0 consolidation on rellib-r12, one
# after another, 8 R processes at a time (adapted from
# dev/rel073-verify.sh): every file of the nine suites one per process
# with every gate on, under the reference BLAS (r074-ref); the files
# the round changed or added (dev/rel074-log/ob.jobs) under OpenBLAS
# 0.3.26 and 0.3.32, BLAS and LAPACK, 4 threads; the escaped-warning
# scan over every file. The one-session runs ran before it
# (dev/release/rel0.74.0.log, rel0.74.0g.log).
#   bash dev/rel074-verify.sh
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/3 TEMP=C:/Users/adf44/AppData/Local/Temp/3
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r12
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$ROOT"
V=dev/rel074-log/verify.txt
bash dev/ciharden-scan.sh r074-ref ref base 8 >> $V 2>&1
echo "$(date '+%F %T') ref done" >> $V
bash dev/ciharden-scan.sh r074-ob26 ob0.3.26-lapack-t4 base 8 dev/rel074-log/ob.jobs >> $V 2>&1
echo "$(date '+%F %T') ob26 done" >> $V
bash dev/ciharden-scan.sh r074-ob32 ob0.3.32-lapack-t4 base 8 dev/rel074-log/ob.jobs >> $V 2>&1
echo "$(date '+%F %T') ob32 done" >> $V
FRMTMB_RELLIB=C:/Users/adf44/source/r/rellib-r12 bash dev/rel074-warnscan.sh dev/ciharden-log/r074-ref.jobs 8 >> $V 2>&1
echo "$(date '+%F %T') warnscan done" >> $V
