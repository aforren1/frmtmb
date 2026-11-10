#!/usr/bin/env bash
# The 0.81.0 consolidation after the test tiers, in phases that keep at
# most three R processes at once (the round's limit):
#   1. the survey (2: one pair and its live-row child) beside the frames
#      and the lanes' key tables (1, one script at a time);
#   2. the translator audit beside the record (2 jobs) beside the floors
#      (1 process at a time);
#   3. the ported brms suite: bin 1 at 3; then the ledger and bin 2 (1)
#      beside the docs check (1); then the regenerated test-brms-suite
#      files, every gate on (3);
#   4. R CMD check, one package at a time (a check and its test process);
#   5. the 27 one-session runs, plain tier (3).
#
#   bash dev/rel081-chain2.sh [phase ...]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r20
L=dev/rel081-log
note() { echo "$(date '+%F %T') CHAIN2 $*" >> $L/chain.txt; }
for ph in ${*:-1 2 3 4 5}; do
  note "phase $ph start"
  case "$ph" in
    1)
      bash dev/rel081-survey.sh frames > $L/frames.txt 2>&1
      ( bash dev/rel081-survey.sh survey ) &
      ( bash dev/rel081-keytables.sh all > $L/keytables.txt 2>&1 ) &
      wait;;
    2)
      ( bash dev/rel081-audit.sh > $L/audit.txt 2>&1 ) &
      ( bash dev/rel081-floor.sh 1 > $L/floor.txt 2>&1 ) &
      wait;;
    3)
      REC_JOBS=3 bash dev/rel081-ledger.sh bin1 > $L/ledger-bin1.txt 2>&1
      ( bash dev/rel081-ledger.sh ledger > $L/ledger-ledger.txt 2>&1
        bash dev/rel081-ledger.sh bin2 > $L/ledger-bin2.txt 2>&1 ) &
      ( bash dev/rel081-docs-check.sh > $L/docs-check.txt 2>&1 ) &
      wait
      : > $L/suite.jobs
      for f in tests/testthat/test-brms-suite-*.R; do echo "frmtmb $f" >> $L/suite.jobs; done
      for f in extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
        echo "frmtmb.sample $f" >> $L/suite.jobs; done
      bash dev/ciharden-scan.sh r081-suite ref base 3 $L/suite.jobs > $L/suite.txt 2>&1;;
    4)
      bash dev/rel081-check.sh frmtmb frmtmb.sample frmtmb.survival > $L/check.txt 2>&1;;
    5)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.81.0 3 > $L/sessions-plain.txt 2>&1
      bash dev/rel081-sessions-sum.sh rel0.81.0 > $L/sessions-rel0.81.0-sum.txt 2>&1;;
    *) note "unknown phase $ph";;
  esac
  note "phase $ph end"
done
note END
