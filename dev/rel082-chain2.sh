#!/usr/bin/env bash
# The 0.82.0 consolidation after the survey and the ledger, in phases that
# keep at most three R processes at once (the round's limit):
#   1. the translator audit (2 jobs) beside the docs check (1);
#   2. R CMD check, one package at a time (a check and its test process),
#      beside the floors (1 process at a time);
#   3. OpenBLAS: the two R copies are built first, alone, then the files
#      run 3 at a time;
#   4. the 27 one-session runs, plain tier (3);
#   5. the 27 one-session runs, every gate on (3).
#
#   bash dev/rel082-chain2.sh [phase ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/rel082-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r21
L=dev/rel082-log
note() { echo "$(date '+%F %T') CHAIN2 $*" >> $L/chain.txt; }
for ph in ${*:-1 2 3 4 5}; do
  note "phase $ph start"
  case "$ph" in
    1)
      ( bash dev/rel082-audit.sh > $L/audit.txt 2>&1 ) &
      ( bash dev/rel082-docs-check.sh > $L/docs-check.txt 2>&1 ) &
      wait;;
    2)
      ( bash dev/rel082-check.sh frmtmb frmtmb.sample frmtmb.spatial frmtmb.survival \
          frmtmb.latent > $L/check.txt 2>&1 ) &
      ( bash dev/rel082-floor.sh 1 > $L/floor.txt 2>&1 ) &
      wait;;
    3)
      bash dev/rel082-ob.sh > $L/ob.txt 2>&1;;
    4)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.82.0 3 > $L/sessions-plain.txt 2>&1
      bash dev/rel082-sessions-sum.sh rel0.82.0 > $L/sessions-rel0.82.0-sum.txt 2>&1;;
    5)
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.82.0g 3 > $L/sessions-gated.txt 2>&1
      bash dev/rel082-sessions-sum.sh rel0.82.0g > $L/sessions-rel0.82.0g-sum.txt 2>&1;;
    *) note "unknown phase $ph";;
  esac
  note "phase $ph end"
done
note END
