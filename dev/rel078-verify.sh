#!/usr/bin/env bash
# The verification runs of the 0.78.0 consolidation on rellib-r16,
# adapted from dev/rel077-verify.sh, one stage per call so that the
# R CMD checks (dev/rel078-check.sh, one at a time) can run beside a
# stage with fewer processes:
#   focus  the round's changed test files (dev/rel078-log/focus.jobs),
#          every gate on
#   base   the same files on rellib-r15 (0.77.0), every gate on: the
#          round's new tests seen to fail before
#   plain  every file of the nine suites, one per process, NOT_CRAN only
#   ref    the same with every gate on (brms, drmTMB, fuzz, scale)
#   ob     the round's numeric files (dev/rel078-log/ob.jobs) under
#          OpenBLAS 0.3.26 and 0.3.32, BLAS and LAPACK, 4 threads
#
#   bash dev/rel078-verify.sh <stage> [N]
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r16
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel078
cd "$ROOT" || exit 1
N="${2:-8}"
V=dev/rel078-log/verify.txt
echo "$(date '+%F %T') $1 start N=$N" >> $V
case "$1" in
  focus) bash dev/ciharden-scan.sh r078-focus ref base "$N" dev/rel078-log/focus.jobs >> $V 2>&1;;
  base) ROUND_BASE=C:/Users/adf44/source/r/rellib-r15 bash dev/ciharden-scan.sh r078-focus-r15 ref base "$N" dev/rel078-log/focus.jobs >> $V 2>&1;;
  plain) CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r078-plain ref base "$N" >> $V 2>&1;;
  ref) bash dev/ciharden-scan.sh r078-ref ref base "$N" >> $V 2>&1;;
  ob)
    bash dev/ciharden-scan.sh r078-ob26 ob0.3.26-lapack-t4 base "$N" dev/rel078-log/ob.jobs >> $V 2>&1
    bash dev/ciharden-scan.sh r078-ob32 ob0.3.32-lapack-t4 base "$N" dev/rel078-log/ob.jobs >> $V 2>&1;;
  *) echo "unknown stage $1"; exit 1;;
esac
echo "$(date '+%F %T') $1 done" >> $V
