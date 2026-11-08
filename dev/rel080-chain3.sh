#!/usr/bin/env bash
# The rest of the evidence on the final build (rellib-r19: core 23:32:17,
# extensions 21:37), one stage after another, at most three R processes.
# Run it after dev/rel080-chain2.sh and key-table batch b4 have ended.
#
#   bash dev/rel080-chain3.sh [stage ...]
#
# Stages: obbuild ob check sessions docs floor survey ledger (default: all,
# in that order). Each stage's start and end go to dev/rel080-log/chain.txt.
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/rel080-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r19
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel080-log
C=$L/chain.txt
note() { echo "$(date '+%F %T') CHAIN3 $*" >> "$C"; }
stages="${*:-obbuild ob check sessions docs floor survey ledger}"
for s in $stages; do
  note "$s start"
  case "$s" in
    obbuild)
      [ -x dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe ] ||
        bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/ob26-build2.txt 2>&1
      bash dev/ciharden-openblas.sh 0.3.32 lapack > $L/ob32-build.txt 2>&1
      note "obbuild exit $?";;
    ob) bash dev/rel080-verify.sh ob 3;;
    check)
      bash dev/rel080-check.sh frmtmb frmtmb.sample frmtmb.latent frmtmb.ode \
        frmtmb.survival > $L/check.txt 2>&1;;
    sessions)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.80.0 3 > $L/sessions-plain.txt 2>&1
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.80.0g 3 > $L/sessions-gated.txt 2>&1
      bash dev/rel080-sessions-sum.sh rel0.80.0 > $L/sessions-rel0.80.0-sum.txt 2>&1
      bash dev/rel080-sessions-sum.sh rel0.80.0g > $L/sessions-rel0.80.0g-sum.txt 2>&1;;
    docs) bash dev/rel080-docs-check.sh > $L/docs-check.txt 2>&1;;
    floor) bash dev/rel080-floor.sh 3 > $L/floor.txt 2>&1;;
    survey) bash dev/rel080-survey.sh;;
    ledger)
      bash dev/rel080-ledger.sh > $L/ledger-run.txt 2>&1
      "$RS" dev/rel080-ledger-diff.R > $L/ledger-diff.txt 2>&1
      "$RS" dev/rel080-bin2-diff.R > $L/bin2-diff.txt 2>&1;;
    *) note "unknown stage $s";;
  esac
  note "$s end"
done
note "END"
