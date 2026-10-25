#!/usr/bin/env bash
# The verification runs of the 0.75.0 consolidation on rellib-r13,
# adapted from dev/rel074-verify.sh, one stage per call so that the
# R CMD checks (dev/rel075-check.sh, one at a time) can run beside a
# stage with fewer processes:
#   plain  every file of the nine suites, one per process, NOT_CRAN only
#   ref    the same with every gate on (brms, drmTMB, fuzz, scale)
#   ob     the round's numeric files (dev/rel075-log/ob.jobs) under
#          OpenBLAS 0.3.26 and 0.3.32, BLAS and LAPACK, 4 threads
#   warn   the escaped-warning scan over every file, every gate on
#
#   bash dev/rel075-verify.sh <stage> [N]
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r13
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel075
cd "$ROOT" || exit 1
N="${2:-8}"
V=dev/rel075-log/verify.txt
echo "$(date '+%F %T') $1 start N=$N" >> $V
case "$1" in
  plain) CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r075-plain ref base "$N" >> $V 2>&1;;
  ref) bash dev/ciharden-scan.sh r075-ref ref base "$N" >> $V 2>&1;;
  ob)
    bash dev/ciharden-scan.sh r075-ob26 ob0.3.26-lapack-t4 base "$N" dev/rel075-log/ob.jobs >> $V 2>&1
    bash dev/ciharden-scan.sh r075-ob32 ob0.3.32-lapack-t4 base "$N" dev/rel075-log/ob.jobs >> $V 2>&1;;
  warn) FRMTMB_RELLIB=C:/Users/adf44/source/r/rellib-r13 bash dev/rel075-warnscan.sh dev/ciharden-log/r075-ref.jobs "$N" >> $V 2>&1;;
  *) echo "unknown stage $1"; exit 1;;
esac
echo "$(date '+%F %T') $1 done" >> $V
