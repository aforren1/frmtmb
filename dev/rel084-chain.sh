#!/usr/bin/env bash
# The 0.84.0 consolidation's runs on rellib-r23 after the focus tier, in
# phases that keep at most TWO R processes at once, so that the lanes'
# key tables can run beside them in a third:
#   1. every file plain, one per process (2);
#   2. every file with every gate on (2);
#   3. the ported brms suite: bin 1 (2 at once), the ledger, gen and
#      guards, bin 2 (1), then the docs check (1);
#   4. the compat frames and the recorded survey (2: one pair at a time);
#   5. the floors (2);
#   6. the translator audit (2);
#   7. R CMD check --as-cran of frmtmb, frmtmb.eam, frmtmb.survival and
#      frmtmb.coupling, one at a time (a check and its test process);
#   8. OpenBLAS 0.3.26 and 0.3.32 with LAPACK, CORETYPE Haswell and unset (2);
#   9. the 27 one-session runs, plain (2);
#  10. the 27 one-session runs, every gate on (2).
# Progress lines go to dev/rel084-log/chain.txt ("CHAIN phase <n> end").
#
#   bash dev/rel084-chain.sh [phase ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel084 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r23
L=dev/rel084-log
note() { echo "$(date '+%F %T') CHAIN $*" >> $L/chain.txt; }
for ph in ${*:-1 2 3 4 5 6 7 8 9 10}; do
  note "phase $ph start"
  case "$ph" in
    1) bash dev/rel084-verify.sh plain 2;;
    2) bash dev/rel084-verify.sh ref 2;;
    3)
      { REC_JOBS=2 bash dev/rel084-ledger.sh bin1; bash dev/rel084-ledger.sh ledger
        bash dev/rel084-ledger.sh gen; bash dev/rel084-bin2.sh; } > $L/ledger-run.txt 2>&1
      bash dev/rel084-docs-check.sh > $L/docs-check.txt 2>&1;;
    4) bash dev/rel084-survey.sh all > $L/survey-run.txt 2>&1;;
    5) bash dev/rel084-floor.sh 2 > $L/floor.txt 2>&1;;
    6) bash dev/rel084-audit.sh > $L/audit.txt 2>&1;;
    7) bash dev/rel084-check.sh frmtmb frmtmb.eam frmtmb.survival frmtmb.coupling > $L/check.txt 2>&1;;
    8) OB_N=2 bash dev/rel084-ob.sh > $L/ob.txt 2>&1;;
    9)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.84.0 2 > $L/sessions-plain.txt 2>&1
      bash dev/rel084-sessions-sum.sh rel0.84.0 > $L/sessions-rel0.84.0-sum.txt 2>&1;;
    10)
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.84.0g 2 > $L/sessions-gated.txt 2>&1
      bash dev/rel084-sessions-sum.sh rel0.84.0g > $L/sessions-rel0.84.0g-sum.txt 2>&1;;
    *) note "unknown phase $ph";;
  esac
  note "phase $ph end"
done
note END
