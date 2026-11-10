#!/usr/bin/env bash
# The verification runs of the 0.81.0 consolidation on rellib-r20,
# adapted from dev/rel080-verify.sh, one stage per call so that the
# R CMD checks (dev/rel081-check.sh, one at a time) can run beside a
# stage with fewer processes:
#   focus  the round's changed test files (dev/rel081-log/focus.jobs),
#          every gate on
#   base   the same files on rellib-r19 (0.80.0), every gate on: the
#          round's new tests seen to fail before
#   plain  every file of the nine suites, one per process, NOT_CRAN only
#   ref    the same with every gate on (brms, drmTMB, fuzz, scale)
#   ob     the round's numeric files (dev/rel081-log/ob.jobs) under
#          OpenBLAS 0.3.26 and 0.3.32, BLAS and LAPACK, 4 threads
#
#   bash dev/rel081-verify.sh <stage> [N]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r20
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel081
cd "$ROOT" || exit 1
N="${2:-8}"
V=dev/rel081-log/verify.txt
echo "$(date '+%F %T') $1 start N=$N" >> $V
case "$1" in
  focus) bash dev/ciharden-scan.sh r081-focus ref base "$N" dev/rel081-log/focus.jobs >> $V 2>&1;;
  base) ROUND_BASE=C:/Users/adf44/source/r/rellib-r19 bash dev/ciharden-scan.sh r081-focus-r19 ref base "$N" dev/rel081-log/focus.jobs >> $V 2>&1;;
  plain) CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r081-plain ref base "$N" >> $V 2>&1;;
  ref) bash dev/ciharden-scan.sh r081-ref ref base "$N" >> $V 2>&1;;
  ob)
    bash dev/ciharden-scan.sh r081-ob26 ob0.3.26-lapack-t4 base "$N" dev/rel081-log/ob.jobs >> $V 2>&1
    bash dev/ciharden-scan.sh r081-ob32 ob0.3.32-lapack-t4 base "$N" dev/rel081-log/ob.jobs >> $V 2>&1;;
  *) echo "unknown stage $1"; exit 1;;
esac
echo "$(date '+%F %T') $1 done" >> $V
