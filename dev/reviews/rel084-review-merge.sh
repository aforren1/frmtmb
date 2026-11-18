#!/usr/bin/env bash
# rel084 release review, merge correctness of the multi-lane files.
# For each file: which lanes changed it against b255350e (LF-normalized),
# whether each lane's copy merged three-way into the release copy leaves
# it unchanged (its hunks are present), and a replay of the merges from
# b255350e in the consolidation's order, compared with the release copy
# (what differs is a hand resolution or a consolidation edit).
# Read-only: git show in the release tree only; lane files are read from
# disk; scratch in $SCR.
#   bash dev/reviews/rel084-review-merge.sh <scratch dir>
set -u
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel084
SCR=$1
mkdir -p "$SCR"
LANES="sparsese gpselimit slopetrig alapse nlband quadesc nbrobust squant drawdir ordstarts cdfname samplers codesize cran"
FILES="R/se-check.R R/fit.R R/families.R R/predict.R R/predict-brms.R R/compat.R R/par-template.R R/fit-restart.R R/confint.R NEWS.md dev/test-backlog.md dev/upstream-bugs.md"
for f in $FILES; do
  k=$(echo "$f" | tr '/.' '__')
  git -C "$REL" show "b255350e:$f" | tr -d '\r' > "$SCR/$k.base"
  tr -d '\r' < "$REL/$f" > "$SCR/$k.rel"
  cp "$SCR/$k.base" "$SCR/$k.replay"
  union=""
  case "$f" in NEWS.md|dev/test-backlog.md|dev/upstream-bugs.md) union="--union";; esac
  touched=""
  for l in $LANES; do
    lf=$R/frmtmb-wt-$l/$f
    [ -f "$lf" ] || continue
    tr -d '\r' < "$lf" > "$SCR/$k.$l"
    if cmp -s "$SCR/$k.$l" "$SCR/$k.base"; then rm "$SCR/$k.$l"; continue; fi
    touched="$touched $l"
    # lane into the release copy
    git merge-file -p $union "$SCR/$k.rel" "$SCR/$k.base" "$SCR/$k.$l" > "$SCR/$k.$l.into" 2>/dev/null
    rc=$?
    if cmp -s "$SCR/$k.$l.into" "$SCR/$k.rel"; then st="present"; else
      st="DIFFERS(rc=$rc,$(diff "$SCR/$k.rel" "$SCR/$k.$l.into" | grep -c '^[<>]') lines)"; fi
    echo "INTO $f $l $st"
    git merge-file $union "$SCR/$k.replay" "$SCR/$k.base" "$SCR/$k.$l" 2>/dev/null
    echo "REPLAY-STEP $f $l rc=$?"
  done
  nd=$(diff "$SCR/$k.replay" "$SCR/$k.rel" | grep -c '^[<>]')
  echo "FILE $f lanes:$touched | replay vs release: $nd differing lines"
  diff "$SCR/$k.replay" "$SCR/$k.rel" > "$SCR/$k.replay.diff"
done
echo DONE
