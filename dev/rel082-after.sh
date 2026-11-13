#!/usr/bin/env bash
# After the audit and the floors: the docs check beside R CMD check (1 + 2
# R processes), then dev/rel082-chain2.sh phases 3 to 5 (OpenBLAS, the
# plain and the gated sessions), so the round stays at three R processes.
#
#   bash dev/rel082-after.sh
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
L=dev/rel082-log
note() { echo "$(date '+%F %T') AFTER $*" >> $L/chain.txt; }
note "waiting for the audit and the floors"
until grep -q "AUDIT END" $L/audit.txt 2>/dev/null && grep -q "FLOOR DONE" $L/floor.txt 2>/dev/null; do sleep 30; done
note "docs check and R CMD check start"
( bash dev/rel082-docs-check.sh > $L/docs-check.txt 2>&1 ) &
( bash dev/rel082-check.sh frmtmb frmtmb.sample frmtmb.spatial frmtmb.survival frmtmb.latent > $L/check.txt 2>&1 ) &
wait
note "docs check and R CMD check end"
bash dev/rel082-chain2.sh 3 4 5
note END
