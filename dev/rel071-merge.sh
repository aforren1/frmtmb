#!/usr/bin/env bash
# Merge the lanes of the round opened 2026-10-13 onto the release tree,
# file by file with git merge-file, in the order dblbar, famsimplex,
# famuni, postfix, optstart, illse, remeasure. Every lane is on d2ef6312
# (frmtmb 0.70.0); the release tree is on 253b3ef7 (0.70.1), which
# differs from d2ef6312 only by lane ci701's commit. So each file starts
# from its 253b3ef7 blob (HEAD) and each lane's change is merged in with
# the d2ef6312 blob as the base. It touches no index and makes no git
# write: blobs are read with `git show`.
#
# Lane illse was still in its fix round when the six-lane tree was
# built, so it is merged only when asked:
#
#   bash dev/rel071-merge.sh           # six lanes, without illse
#   bash dev/rel071-merge.sh illse     # all seven
#
# Rerunnable, as dev/rel070-merge.sh: every run first restores each file
# any lane changed (of all seven, so a run without illse also removes
# what a run with it merged) to its HEAD blob and removes the untracked
# lane files an earlier run copied ($LOG/copied-<lane>.lst), then
# merges the lanes in order and copies their untracked files again.
# The hand resolutions are applied by dev/rel071-resolve.R and the NEWS
# files are written by dev/rel071-news-merge.R, both run at the end.
#
# Line endings are normalized to LF first, as git stores them.
# Generated files (man/*.Rd, NAMESPACE, frmtmb.sample's copy of
# helper-brms-suite.R, and R/links-brms.R, which dev/famlink-gen-links.R
# writes) are left out; roxygen and the generator write them. Untracked
# lane files are copied, except generated Rd files, Stan caches and
# OpenBLAS builds of R.
set -u
BASE=d2ef6312
HEADC=253b3ef7
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LOG=$REL/dev/rel071-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ALL="dblbar famsimplex famuni postfix optstart illse remeasure"
if [ "${1:-}" = illse ]; then
  LANES="$ALL"
else
  LANES="dblbar famsimplex famuni postfix optstart remeasure"
fi
export REL071_LANES="$LANES"
cd "$REL" || exit 1
[ "$(git rev-parse --short=8 HEAD)" = "$HEADC" ] ||
  { echo "release tree is not on $HEADC"; exit 1; }
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|R/links-brms\.R|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
# dev/ciharden-out/ holds a lane's OpenBLAS copies of R (925 MB in lane
# illse's tree); it is gitignored, and listed here so that no change to
# the ignore file can let it in
uskip='^(man/.*\.Rd|extensions/[^/]*/man/.*\.Rd|.*stan-cache.*|dev/ciharden-out/.*)$'
echo "lanes: $LANES"

# every tracked file any lane of the round changed, illse included
: > "$W/all"
for l in $ALL; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  [ "$(git -C "$src" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "lane $l is not on $BASE"; exit 1; }
  git -C "$src" diff --name-status "$BASE" 2>/dev/null > "$W/$l.status"
  cut -f2 "$W/$l.status" >> "$W/all"
done
sort -u "$W/all" | grep -vE "$excl" > "$W/files"
echo "files changed by some lane (without generated ones): $(wc -l < "$W/files")"

# restore the HEAD blob of each, and remove what an earlier run copied
while IFS= read -r f; do
  git show "HEAD:$f" | tr -d '\r' > "$f"
done < "$W/files"
git show "HEAD:R/links-brms.R" | tr -d '\r' > R/links-brms.R
for l in $ALL; do
  m="$LOG/copied-$l.lst"
  [ -f "$m" ] || continue
  while IFS= read -r f; do
    # never remove a tracked file
    git ls-files --error-unmatch "$f" > /dev/null 2>&1 && continue
    rm -f "$f"
  done < "$m"
  rm -f "$m"
done

for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  echo "== lane $l"
  nc=0; nk=0
  while IFS=$'\t' read -r st f; do
    echo "$f" | grep -qE "$excl" && continue
    case "$st" in
      M)
        git show "$BASE:$f" | tr -d '\r' > "$W/base"
        tr -d '\r' < "$src/$f" > "$W/theirs"
        tr -d '\r' < "$f" > "$W/ours"
        k="$W/keep/$l.$(echo "$f" | tr '/' '_')"
        cp "$W/ours" "$k.ours"; cp "$W/base" "$k.base"
        cp "$W/theirs" "$k.theirs"
        if git merge-file -L release -L base -L "$l" \
             "$W/ours" "$W/base" "$W/theirs"; then
          rm -f "$k.ours" "$k.base" "$k.theirs"; nc=$((nc + 1))
        else
          echo "   CONFLICT($?) $f"; nk=$((nk + 1))
        fi
        cp "$W/ours" "$f";;
      D) echo "   UNHANDLED deletion $f";;
      *) echo "   UNHANDLED $st $f";;
    esac
  done < "$W/$l.status"
  echo "   merged clean: $nc, with conflicts: $nk"
  # untracked lane files
  lst="$W/$l.untracked"; : > "$lst"; clash=0; skip=0
  while IFS= read -r f; do
    if echo "$f" | grep -qE "$uskip"; then skip=$((skip + 1)); continue; fi
    if [ -e "$f" ]; then
      echo "   CLASH untracked: $f"; clash=$((clash + 1)); continue
    fi
    printf '%s\n' "$f" >> "$lst"
  done < <(git -C "$src" ls-files --others --exclude-standard)
  ( cd "$src" && tar -cf - -T "$lst" ) | tar -xpf -
  cp "$lst" "$LOG/copied-$l.lst"
  echo "   untracked copied: $(grep -c '' "$lst"), skipped (Rd, caches): $skip," \
       "clashes: $clash"
done
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"
grep -ln '^<<<<<<< ' $(cat "$W/files") 2>/dev/null |
  sed 's/^/   conflicted after the textual merge: /'
[ "${REL071_TEXTUAL_ONLY:-}" = true ] && exit 0

# the consolidation's edits, then the NEWS files
"$RS" dev/rel071-resolve.R || exit 1
"$RS" dev/rel071-news-merge.R || exit 1
grep -ln '^<<<<<<< \|^>>>>>>> ' $(cat "$W/files") R/links-brms.R && {
  echo "CONFLICT MARKERS LEFT"; exit 1; }
echo "MERGE, RESOLVE AND NEWS DONE"
