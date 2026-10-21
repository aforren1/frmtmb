#!/usr/bin/env bash
# Merge each file that two or more lanes changed in the reverse lane
# order onto its 4854e425 blob and compare with the release tree: an
# identical result means no hunk of one lane depends on where another
# lane's landed. NEWS files and generated files are left out. A file
# the resolve edits after the merge differs by exactly those edits, so
# the diff is printed for reading. Adapted from dev/rel072-revorder.sh.
#
#   bash dev/rel073-revorder.sh
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/rel073-log/revorder
BASE=4854e425
LANES="fuzzfix nb2dm cdfs loodraws"
cd "$REL"; rm -rf "$W"; mkdir -p "$W"
: > "$W/all"
for l in $LANES; do
  git -C ../frmtmb-wt-$l diff --name-only $BASE 2>/dev/null >> "$W/all"
done
sort "$W/all" | uniq -d | grep -vE 'NEWS\.md$|man/|NAMESPACE$|links-brms\.R$' > "$W/files"
echo "files changed by two or more lanes: $(wc -l < "$W/files")"
while IFS= read -r f; do
  m="$W/$(echo "$f" | tr '/' '_')"
  git show "$BASE:$f" | tr -d '\r' > "$m"
  who=""
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
    diff "$m" "$m.rel" | sed 's/^/    /' | head -80
  fi
done < "$W/files"
