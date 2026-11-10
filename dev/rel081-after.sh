#!/usr/bin/env bash
# What follows chain2 (dev/rel081-chain2.sh), started by the second
# consolidator: it waits for chain2's outer shell to exit, so the round
# stays at three R processes, then runs OpenBLAS (dev/rel081-ob.sh, 3 at
# a time) and the 27 one-session runs of the gated tier (3 at a time).
# Progress lines go to dev/rel081-log/chain.txt.
#
#   bash dev/rel081-after.sh <chain2 pid, as tasklist shows it>
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
L=dev/rel081-log
note() { echo "$(date '+%F %T') AFTER $*" >> $L/chain.txt; }
pid="${1:?chain2 pid}"
note "waiting for chain2 (pid $pid)"
while tasklist //FI "PID eq $pid" 2>/dev/null | grep -q " $pid "; do sleep 30; done
note "chain2 gone; ob start"
bash dev/rel081-ob.sh > $L/ob.txt 2>&1
note "ob end"
note "gated sessions start"
NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
  FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.81.0g 3 > $L/sessions-gated.txt 2>&1
bash dev/rel081-sessions-sum.sh rel0.81.0g > $L/sessions-rel0.81.0g-sum.txt 2>&1
note "gated sessions end"
note END
