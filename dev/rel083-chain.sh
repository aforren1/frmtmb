#!/usr/bin/env bash
# The 0.83.0 consolidation's runs on rellib-r22 after the focus tier, in
# phases that keep at most three R processes at once:
#   1. every file plain, one per process (2) beside the ported brms
#      suite: bin 1, the ledger, gen and guards, bin 2 (1);
#   2. every file with every gate on (2) beside the docs check (1);
#   3. the compat frames and the recorded survey (2) beside the floors (1);
#   4. the translator audit (2);
#   5. R CMD check --as-cran of frmtmb, frmtmb.sample, frmtmb.spatial and
#      frmtmb.latent, one at a time (a check and its test process);
#   6. OpenBLAS 0.3.26 and 0.3.32 with LAPACK, CORETYPE Haswell and unset (3);
#   7. the 27 one-session runs, plain (3);
#   8. the 27 one-session runs, every gate on (3).
# Progress lines go to dev/rel083-log/chain.txt ("CHAIN phase <n> end").
#
#   bash dev/rel083-chain.sh [phase ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel083 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r22
L=dev/rel083-log
note() { echo "$(date '+%F %T') CHAIN $*" >> $L/chain.txt; }
for ph in ${*:-1 2 3 4 5 6 7 8}; do
  note "phase $ph start"
  case "$ph" in
    1)
      ( bash dev/rel083-verify.sh plain 2 ) &
      ( { REC_JOBS=1 bash dev/rel083-ledger.sh bin1; bash dev/rel083-ledger.sh ledger
          bash dev/rel083-ledger.sh gen; bash dev/rel083-bin2.sh; } > $L/ledger-run.txt 2>&1 ) &
      wait;;
    2)
      ( bash dev/rel083-verify.sh ref 2 ) &
      ( bash dev/rel083-docs-check.sh > $L/docs-check.txt 2>&1 ) &
      wait;;
    3)
      ( bash dev/rel083-survey.sh all > $L/survey-run.txt 2>&1 ) &
      ( bash dev/rel083-floor.sh 1 > $L/floor.txt 2>&1 ) &
      wait;;
    4) bash dev/rel083-audit.sh > $L/audit.txt 2>&1;;
    5) bash dev/rel083-check.sh frmtmb frmtmb.sample frmtmb.spatial frmtmb.latent > $L/check.txt 2>&1;;
    6) bash dev/rel083-ob.sh > $L/ob.txt 2>&1;;
    7)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.83.0 3 > $L/sessions-plain.txt 2>&1
      bash dev/rel083-sessions-sum.sh rel0.83.0 > $L/sessions-rel0.83.0-sum.txt 2>&1;;
    8)
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.83.0g 3 > $L/sessions-gated.txt 2>&1
      bash dev/rel083-sessions-sum.sh rel0.83.0g > $L/sessions-rel0.83.0g-sum.txt 2>&1;;
    *) note "unknown phase $ph";;
  esac
  note "phase $ph end"
done
note END
