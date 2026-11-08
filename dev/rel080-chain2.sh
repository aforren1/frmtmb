#!/usr/bin/env bash
# The two full tiers again on the final build (core reinstalled 23:32:17 after
# the bigx x nlstart fix), 2 R processes each, so that one slot stays free.
#   bash dev/rel080-chain2.sh
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
for s in ref plain; do
  echo "$(date '+%F %T') CHAIN2 $s start" >> dev/rel080-log/chain.txt
  bash dev/rel080-verify.sh "$s" 2
  echo "$(date '+%F %T') CHAIN2 $s end" >> dev/rel080-log/chain.txt
done
echo "$(date '+%F %T') CHAIN2 END" >> dev/rel080-log/chain.txt
