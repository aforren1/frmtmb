#!/usr/bin/env bash
# Short check of 0.82.0: does the release tree hold every lane's change to
# a path that more than one lane takes? For each lane and each such path,
# a three-way merge of the lane's copy into the release copy over the base
# (e8abed0a) must give the release copy back. Read-only: git show only.
set -u
export GIT_OPTIONAL_LOCKS=0
R=C:/Users/adf44/source/r
REL=$R/frmtmb-wt-rel082
O=$REL/dev/reviews/rel082-review-out/multilane
mkdir -p $O/tmp
LANES="boundfam quadseam ppdraw wcheck priortab nlstop gpstop ordstart silent arspat names sampfront"
cat $REL/dev/rel082-log/commit-dry-{boundfam,quadseam,ppdraw,wcheck,priortab,nlstop,gpstop,ordstart,silent,arspat,names,sampfront}.txt |
  sed 's/^ *//' | sort | uniq -c | awk '$1 > 1 {print $2}' > $O/multi.txt
: > $O/result.txt
for l in $LANES; do
  sed 's/^ *//' $REL/dev/rel082-log/commit-dry-$l.txt | grep -Fxf $O/multi.txt | while read -r p; do
    a=$R/frmtmb-wt-$l/$p; b=$REL/$p
    if [ ! -e "$a" ]; then
      [ -e "$b" ] && echo "LANE-DELETED-REL-HAS $l $p" >> $O/result.txt
      continue
    fi
    tr -d '\r' < "$a" > $O/tmp/lane
    tr -d '\r' < "$b" > $O/tmp/rel 2>/dev/null || : > $O/tmp/rel
    git -C $REL show "e8abed0a:$p" 2>/dev/null | tr -d '\r' > $O/tmp/base
    git merge-file -p $O/tmp/rel $O/tmp/base $O/tmp/lane > $O/tmp/out 2>/dev/null
    rc=$?
    if [ $rc -ne 0 ]; then echo "CONFLICT($rc) $l $p" >> $O/result.txt
    elif ! cmp -s $O/tmp/out $O/tmp/rel; then echo "MISSING $l $p" >> $O/result.txt
    else echo "ok $l $p" >> $O/result.txt; fi
  done
done
grep -vc "^ok" $O/result.txt
grep -c "^ok" $O/result.txt
