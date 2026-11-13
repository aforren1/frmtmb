#!/usr/bin/env bash
# The reruns after the release review's fixes and lane boundfam's
# consolidation fixes, at most three R processes at once:
#   A. every file of the nine suites, one per process, NOT_CRAN only (3);
#   B. the same with every gate on (3);
#   C. R CMD check of frmtmb, frmtmb.eam and frmtmb.ode (2: a check and its
#      test process) beside the docs check (1);
#   D. OpenBLAS 0.3.26 with LAPACK on the changed numeric files: the R copy
#      is built alone, then 3 files at a time;
#   E. frmtmb's three one-session runs, plain, then every gate on (3).
#
#   bash dev/rel082-chain3.sh [stage ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/rel082-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r21
L=dev/rel082-log
note() { echo "$(date '+%F %T') CHAIN3 $*" >> $L/chain.txt; }
for s in ${*:-A B C D E}; do
  note "stage $s start"
  case "$s" in
    A) CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r082-plain2 ref base 3 > $L/plain2.txt 2>&1;;
    B) bash dev/ciharden-scan.sh r082-ref2 ref base 3 > $L/ref2.txt 2>&1;;
    C)
      ( bash dev/rel082-check.sh frmtmb frmtmb.eam frmtmb.ode > $L/check2.txt 2>&1 ) &
      ( bash dev/rel082-docs-check.sh > $L/docs-check2.txt 2>&1 ) &
      wait;;
    D)
      bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build2.txt 2>&1
      echo "build exit $?" >> $L/ob2.txt
      bash dev/ciharden-scan.sh r082-ob26fix ob0.3.26-lapack-t4 base 3 $L/ob-fix.jobs >> $L/ob2.txt 2>&1
      rm -rf dev/ciharden-out;;
    E)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.82.0c 3 frmtmb > $L/sessions-core-plain.txt 2>&1
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.82.0cg 3 frmtmb > $L/sessions-core-gated.txt 2>&1;;
    *) note "unknown stage $s";;
  esac
  note "stage $s end"
done
note END
