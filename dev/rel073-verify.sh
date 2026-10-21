#!/usr/bin/env bash
# The verification runs of the 0.73.0 consolidation on the merged build
# (wt-release-lib), one after another, 8 R processes at a time, after
# the reference-BLAS scan (r073-ref): all eight suites one file per
# process under OpenBLAS 0.3.26 and 0.3.32 (BLAS and LAPACK, 4 threads),
# every gate on (dev/ciharden-scan.sh); the escaped-warning scan
# (dev/rel073-warnscan.sh); the one-session runs, forward, reverse and
# preload, ungated (NOT_CRAN only) and gated (every gate the scans set).
# ROUND_BASE is the 0.72.0 library: rellib-r11 does not exist yet.
#   bash dev/rel073-verify.sh
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/3 TEMP=C:/Users/adf44/AppData/Local/Temp/3
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r10
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
LIB=C:/Users/adf44/source/r/wt-release-lib
cd "$ROOT"
V=dev/rel073-log/verify.txt
bash dev/ciharden-scan.sh r073-ob26 ob0.3.26-lapack-t4 $LIB 8 >> $V 2>&1
echo "$(date '+%F %T') ob26 done" >> $V
bash dev/ciharden-scan.sh r073-ob32 ob0.3.32-lapack-t4 $LIB 8 >> $V 2>&1
echo "$(date '+%F %T') ob32 done" >> $V
FRMTMB_RELLIB=$LIB bash dev/rel073-warnscan.sh dev/ciharden-log/r073-ref.jobs 8 >> $V 2>&1
echo "$(date '+%F %T') warnscan done" >> $V
( export NOT_CRAN=true FRMTMB_RELLIB=$LIB
  bash dev/release/run-sessions.sh rel0.73.0 8 >> $V 2>&1 )
echo "$(date '+%F %T') sessions ungated done" >> $V
( export NOT_CRAN=true FRMTMB_RELLIB=$LIB FRMTMB_BRMS_FIT_TESTS=true \
    FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
  bash dev/release/run-sessions.sh rel0.73.0g 8 >> $V 2>&1 )
echo "$(date '+%F %T') sessions gated done" >> $V
