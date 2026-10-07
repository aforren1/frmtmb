#!/usr/bin/env bash
# Merge the four lanes of the round of 2026-10-06 onto the release tree
# in the order ciharden, surface, setier, optima, file by file with
# git merge-file against the base blob (4f5ea39f, which every lane is
# on). It touches no index and makes no git write: the base blob is
# read with `git show`.
#
# Rerunnable, in the manner of dev/rel068-nanse-merge.sh: every run
# first restores each file any lane changed to its base blob, then
# merges the lanes in order, so a replacement file from a lane is
# merged again by rerunning. The hand resolutions are applied by
# dev/rel069-resolve.R and the NEWS files are written by
# dev/rel069-news-merge.R, both run at the end.
#
# Line endings are normalized to LF first, as git stores them.
# Generated files (man/*.Rd, NAMESPACE, frmtmb.sample's copy of
# helper-brms-suite.R) are left out; roxygen writes them. Untracked
# lane files are copied without overwriting, except generated Rd files
# and Stan caches.
#
#   bash dev/rel069-merge.sh
set -u
BASE=4f5ea39f
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LOG=$REL/dev/rel069-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LANES="ciharden surface setier optima"
cd "$REL" || exit 1
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'

# every tracked file any lane changed
: > "$W/all"
for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  [ "$(git -C "$src" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "lane $l is not on $BASE"; exit 1; }
  git -C "$src" diff --name-status "$BASE" 2>/dev/null > "$W/$l.status"
  cut -f2 "$W/$l.status" >> "$W/all"
done
sort -u "$W/all" | grep -vE "$excl" > "$W/files"
echo "files changed by some lane (without generated ones): $(wc -l < "$W/files")"

# restore the base blob of each
while IFS= read -r f; do
  git show "$BASE:$f" | tr -d '\r' > "$f"
done < "$W/files"

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
  # untracked lane files: copy without overwriting
  lst="$W/$l.untracked"; : > "$lst"; clash=0; skip=0; same=0
  while IFS= read -r f; do
    case "$f" in
      man/*.Rd|extensions/*/man/*.Rd|*stan-cache*/*) skip=$((skip + 1))
        continue;;
    esac
    if [ -e "$f" ]; then
      if cmp -s "$src/$f" "$f"; then same=$((same + 1)); else
        echo "   CLASH untracked: $f"; clash=$((clash + 1)); fi
      continue
    fi
    printf '%s\n' "$f" >> "$lst"
  done < <(git -C "$src" ls-files --others --exclude-standard)
  ( cd "$src" && tar -cf - -T "$lst" ) | tar -xpf -
  echo "   untracked copied: $(grep -c '' "$lst"), already here and" \
       "identical: $same, skipped (Rd, caches): $skip, clashes: $clash"
done
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"

# the consolidation's edits, then the NEWS files
"$RS" dev/rel069-resolve.R || exit 1
"$RS" dev/rel069-news-merge.R || exit 1
grep -ln '^<<<<<<< \|^>>>>>>> ' $(cat "$W/files") && {
  echo "CONFLICT MARKERS LEFT"; exit 1; }
echo "MERGE, RESOLVE AND NEWS DONE"
