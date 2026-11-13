#!/usr/bin/env bash
# The reruns after lane boundfam's OpenBLAS fix (its findings,
# "Consolidation fixes", item 3), at most three R processes at once:
#   1. test-se-limit.R under the reference BLAS and under OpenBLAS 0.3.26
#      and 0.3.32 with LAPACK (each copy built alone first);
#   2. every file of the nine suites, one per process, NOT_CRAN only (3);
#   3. the same with every gate on (3);
#   4. R CMD check of frmtmb (a check and its test process);
#   5. frmtmb's three one-session runs, plain, then every gate on (3).
#
#   bash dev/rel082-chain4.sh [stage ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/rel082-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r21
L=dev/rel082-log
note() { echo "$(date '+%F %T') CHAIN4 $*" >> $L/chain.txt; }
for s in ${*:-1 2 3 4 5}; do
  note "stage $s start"
  case "$s" in
    1)
      bash dev/ciharden-scan.sh r082-sl-ref ref base 1 $L/selimit.jobs > $L/sl.txt 2>&1
      for v in 0.3.26 0.3.32; do
        bash dev/ciharden-openblas.sh $v lapack > $L/ob-build-sl-$v.txt 2>&1
        echo "build $v exit $?" >> $L/sl.txt
        bash dev/ciharden-scan.sh r082-sl-ob${v##*.} ob$v-lapack-t4 base 1 $L/selimit.jobs >> $L/sl.txt 2>&1
      done
      rm -rf dev/ciharden-out;;
    2) CIHARDEN_UNGATED=true bash dev/ciharden-scan.sh r082-plain3 ref base 3 > $L/plain3.txt 2>&1;;
    3) bash dev/ciharden-scan.sh r082-ref3 ref base 3 > $L/ref3.txt 2>&1;;
    4) bash dev/rel082-check.sh frmtmb > $L/check3.txt 2>&1;;
    5)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.82.0d 3 frmtmb > $L/sessions-core-plain2.txt 2>&1
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.82.0dg 3 frmtmb > $L/sessions-core-gated2.txt 2>&1;;
    *) note "unknown stage $s";;
  esac
  note "stage $s end"
done
note END
