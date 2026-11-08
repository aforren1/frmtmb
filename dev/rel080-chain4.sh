#!/usr/bin/env bash
# The release review's fixes (dev/rel080-release-fix3.R; core reinstalled
# 2026-11-07 08:59:26): the tiers and sessions again, 2 R processes.
#   bash dev/rel080-chain4.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/rel080-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
L=dev/rel080-log
note() { echo "$(date '+%F %T') CHAIN4 $*" >> $L/chain.txt; }
for s in ref plain; do
  note "$s start"; bash dev/rel080-verify.sh "$s" 2; note "$s end"
done
note "sessions start"
NOT_CRAN=true bash dev/release/run-sessions.sh rel0.80.0 2 > $L/sessions-plain.txt 2>&1
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
  FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.80.0g 2 > $L/sessions-gated.txt 2>&1
bash dev/rel080-sessions-sum.sh rel0.80.0 > $L/sessions-rel0.80.0-sum.txt 2>&1
bash dev/rel080-sessions-sum.sh rel0.80.0g > $L/sessions-rel0.80.0g-sum.txt 2>&1
note "sessions end"
note END
