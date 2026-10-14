#!/usr/bin/env bash
# Merge each file that two or more sources changed (ci701 in HEAD
# counts as one) in the reverse lane order onto its HEAD blob and
# compare with the release tree: an identical result means no hunk of
# one lane depends on where another lane's landed. NEWS files and
# dev/test-backlog.md are rebuilt by scripts and are left out. A file
# the resolve edits after the merge differs by exactly those edits, so
# the diff is printed for reading.
#
#   bash dev/rel071-revorder.sh [illse]
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/rel071-log/revorder
BASE=d2ef6312
if [ "${1:-}" = illse ]; then
  LANES="remeasure illse optstart postfix famuni famsimplex dblbar"
else
  LANES="remeasure optstart postfix famuni famsimplex dblbar"
fi
cd "$REL"; rm -rf "$W"; mkdir -p "$W"
ci701=$(git diff --name-only $BASE HEAD)
: > "$W/all"
for l in $LANES; do
  git -C ../frmtmb-wt-$l diff --name-only $BASE 2>/dev/null >> "$W/all"
done
echo "$ci701" >> "$W/all"
sort "$W/all" | uniq -d | grep -vE 'NEWS\.md$|test-backlog\.md$|man/|NAMESPACE$|links-brms\.R$' > "$W/files"
echo "files changed by two or more sources: $(wc -l < "$W/files")"
while IFS= read -r f; do
  m="$W/$(echo "$f" | tr '/' '_')"
  git show "HEAD:$f" | tr -d '\r' > "$m"
  who=""
  echo "$ci701" | grep -qx "$f" && who="ci701(HEAD)"
  for l in $LANES; do
    src=/c/Users/adf44/source/r/frmtmb-wt-$l
    git -C "$src" diff --quiet $BASE -- "$f" 2>/dev/null && continue
    git show "$BASE:$f" | tr -d '\r' > "$W/base"
    tr -d '\r' < "$src/$f" > "$W/theirs"
    git merge-file -L release -L base -L "$l" "$m" "$W/base" "$W/theirs"
    rc=$?
    who="$who $l"; [ $rc -ne 0 ] && who="$who(conflicts=$rc)"
  done
  tr -d '\r' < "$f" > "$m.rel"
  if cmp -s "$m" "$m.rel"; then echo "SAME: $f [$who ]"
  else
    echo "DIFFERENT: $f [$who ] ($(diff "$m" "$m.rel" | grep -c '^[<>]') lines)"
    diff "$m" "$m.rel" | sed 's/^/    /' | head -60
  fi
done < "$W/files"
