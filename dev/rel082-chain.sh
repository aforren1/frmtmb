#!/usr/bin/env bash
# The test tiers of the 0.82.0 consolidation, one after the other, 3 R
# processes each (the round's limit): the changed files on rellib-r21 and
# on rellib-r20, then every file plain and gated.
#   bash dev/rel082-chain.sh [stages]
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
for s in ${*:-focus base plain ref}; do
  echo "$(date '+%F %T') CHAIN $s start" >> dev/rel082-log/chain.txt
  bash dev/rel082-verify.sh "$s" "${CHAIN_N:-3}"
  echo "$(date '+%F %T') CHAIN $s end" >> dev/rel082-log/chain.txt
done
echo "$(date '+%F %T') CHAIN END" >> dev/rel082-log/chain.txt
