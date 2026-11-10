#!/usr/bin/env bash
# The test tiers of the 0.81.0 consolidation, one after the other, 3 R
# processes each (the round's limit): the changed files on rellib-r20 and
# on rellib-r19, then every file plain and gated.
#   bash dev/rel081-chain.sh [stages]
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
for s in ${*:-focus base plain ref}; do
  echo "$(date '+%F %T') CHAIN $s start" >> dev/rel081-log/chain.txt
  bash dev/rel081-verify.sh "$s" 3
  echo "$(date '+%F %T') CHAIN $s end" >> dev/rel081-log/chain.txt
done
echo "$(date '+%F %T') CHAIN END" >> dev/rel081-log/chain.txt
