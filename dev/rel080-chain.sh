#!/usr/bin/env bash
# The long verification runs of the 0.80.0 consolidation, one after the
# other, 2 R processes each, so that one slot stays free under the round's
# limit of three. Each stage writes dev/ciharden-log/<run>.sum.
#   bash dev/rel080-chain.sh
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
for s in base ref plain; do
  echo "$(date '+%F %T') CHAIN $s start" >> dev/rel080-log/chain.txt
  bash dev/rel080-verify.sh "$s" 2
  echo "$(date '+%F %T') CHAIN $s end" >> dev/rel080-log/chain.txt
done
echo "$(date '+%F %T') CHAIN END" >> dev/rel080-log/chain.txt
