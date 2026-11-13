#!/usr/bin/env bash
# Runs test files one per process in sequence. usage: chain.sh <tag> <pkg>:<file> ...
tag=$1; shift
cd C:/Users/adf44/source/r/frmtmb-wt-rel082
for j in "$@"; do
  p=${j%%:*}; f=${j#*:}
  b=$(basename "$f" .R)
  bash dev/reviews/rel082-review-r.sh dev/reviews/rel082-review-test1.R \
    dev/reviews/rel082-review-out/t-$tag-$b.log "$p" "$f"
  grep -h "^RESULT\|^lib:" dev/reviews/rel082-review-out/t-$tag-$b.log >> dev/reviews/rel082-review-out/chain-$tag.sum
done
echo END >> dev/reviews/rel082-review-out/chain-$tag.sum
