#!/usr/bin/env bash
# Merge the three lanes of the round opened 2026-10-07 onto the release
# tree in the order moparity, optrap, seflags, file by file with
# git merge-file. Every lane is on 48c1e6c9 (frmtmb 0.69.0); the release
# tree is on 2884ccbd (0.69.1), which differs from 48c1e6c9 only by lane
# ci691's commit. So each file starts from its 2884ccbd blob (HEAD) and
# each lane's change is merged in with the 48c1e6c9 blob as the base.
# It touches no index and makes no git write: blobs are read with
# `git show`.
#
# Rerunnable, in the manner of dev/rel069-merge.sh: every run first
# restores each file any lane changed to its HEAD blob and removes the
# untracked lane files an earlier run copied (listed in
# $LOG/copied-<lane>.lst), then merges the lanes in order and copies
# their untracked files again. So a replacement file in a lane's tree is
# merged by rerunning. The hand resolutions are applied by
# dev/rel070-resolve.R and the NEWS files are written by
# dev/rel070-news-merge.R, both run at the end.
#
# Line endings are normalized to LF first, as git stores them.
# Generated files (man/*.Rd, NAMESPACE, frmtmb.sample's copy of
# helper-brms-suite.R) are left out; roxygen writes them. Untracked
# lane files are copied, except generated Rd files, Stan caches and the
# trial-merge copies of the whole package (the optrap reviewer's
# dev/optrap-rev*-merge*/, 35 MB, and lane optrap's scratch copy of
# this release tree for its consolidation fixes, dev/optrap-relmerge/,
# 41 MB).
#
#   bash dev/rel070-merge.sh
set -u
BASE=48c1e6c9
HEADC=2884ccbd
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LOG=$REL/dev/rel070-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LANES="moparity optrap seflags"
cd "$REL" || exit 1
[ "$(git rev-parse --short=8 HEAD)" = "$HEADC" ] ||
  { echo "release tree is not on $HEADC"; exit 1; }
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
uskip='^(man/.*\.Rd|extensions/[^/]*/man/.*\.Rd|.*stan-cache.*|dev/optrap-rev[0-9]*-merge[0-9]*/.*|dev/optrap-relmerge/.*)$'

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

# restore the HEAD blob of each, and remove what an earlier run copied
while IFS= read -r f; do
  git show "HEAD:$f" | tr -d '\r' > "$f"
done < "$W/files"
for l in $LANES; do
  m="$LOG/copied-$l.lst"
  [ -f "$m" ] || continue
  while IFS= read -r f; do
    # never remove a tracked file
    git ls-files --error-unmatch "$f" > /dev/null 2>&1 && continue
    rm -f "$f"
  done < "$m"
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
  echo "   untracked copied: $(grep -c '' "$lst"), skipped (Rd, caches," \
       "trial merges): $skip, clashes: $clash"
done
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"

# the consolidation's edits, then the NEWS files
"$RS" dev/rel070-resolve.R || exit 1
"$RS" dev/rel070-news-merge.R || exit 1
grep -ln '^<<<<<<< \|^>>>>>>> ' $(cat "$W/files") && {
  echo "CONFLICT MARKERS LEFT"; exit 1; }
echo "MERGE, RESOLVE AND NEWS DONE"
