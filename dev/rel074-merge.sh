#!/usr/bin/env bash
# Merge the three lanes of the round opened 2026-10-23 onto the release
# tree, file by file with git merge-file, in the order splitsurv,
# mvnlatent, spatial. Every lane and the release tree are on 19b162b9
# (frmtmb 0.73.1), so each file starts from its 19b162b9 blob and each
# lane's change is merged in with the same blob as the base. It touches
# no index and makes no git write: blobs are read with `git show`.
# Adapted from dev/rel073-merge.sh.
#
#   bash dev/rel074-merge.sh
#
# Rerunnable: every run first restores each file any lane changed or
# deleted to its HEAD blob and removes the untracked lane files an
# earlier run copied ($LOG/copied-<lane>.lst), then merges the lanes in
# order and copies their untracked files again. The hand resolutions are
# applied by dev/rel074-resolve.R and the NEWS files are written by
# dev/rel074-news-merge.R, both run at the end.
#
# This round is the first that deletes tracked files (frmtmb.spline is
# retired, three test files moved to extensions). A deletion is applied
# when the file is still the base blob; a deletion of a file an earlier
# lane changed, or a change to a file an earlier lane deleted, is
# reported as a conflict. Deletions apply to generated files too, since
# a retired package's man/ has no generator left.
#
# Line endings are normalized to LF first, as git stores them. Generated
# files are left out of the textual merge, because their generators run
# after it: man/*.Rd and NAMESPACE (roxygen), frmtmb.sample's copy of
# helper-brms-suite.R, and R/links-brms.R, which dev/famlink-gen-links.R
# writes from brms (lane splitsurv changed the generator: cox() left
# core). Untracked lane files are copied, except generated Rd files, Stan
# caches, the brms suite, OpenBLAS builds of R, the reviewers' scripts
# (dev/review-*) and the round's shared dev/core-boundary.md, which the
# release tree carries itself.
set -u
BASE=19b162b9
REL=${REL074_ROOT:-C:/Users/adf44/source/r/frmtmb-wt-release}
# the checkout whose blobs are read: the release tree itself, or it for
# a scratch copy made with git archive (REL074_ROOT), which has no git
G=${REL074_GIT:-$REL}
export REL074_ROOT="$REL" REL074_GIT="$G"
LOG=$REL/dev/rel074-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LANES="splitsurv mvnlatent spatial"
export REL074_LANES="$LANES"
cd "$REL" || exit 1
[ "$(git -C "$G" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "release tree is not on $BASE"; exit 1; }
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|R/links-brms\.R|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
uskip='^(man/.*\.Rd|extensions/[^/]*/man/.*\.Rd|.*stan-cache.*|dev/ciharden-out/.*|dev/brms-suite/.*|dev/review-.*|dev/core-boundary\.md)$'
echo "lanes: $LANES"

: > "$W/all"
for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  [ "$(git -C "$src" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "lane $l is not on $BASE"; exit 1; }
  git -C "$src" diff --name-status "$BASE" 2>/dev/null > "$W/$l.status"
  cut -f2 "$W/$l.status" >> "$W/all"
done
sort -u "$W/all" > "$W/all-files"
grep -vE "$excl" "$W/all-files" > "$W/files"
echo "files changed or deleted by some lane: $(wc -l < "$W/all-files")," \
     "without generated ones: $(wc -l < "$W/files")"

# restore every file any lane changed or deleted, generated ones too, so
# that a rerun starts from the base
while IFS= read -r f; do
  mkdir -p "$(dirname "$f")"
  git -C "$G" show "HEAD:$f" | tr -d '\r' > "$f"
done < "$W/all-files"
for l in $LANES; do
  m="$LOG/copied-$l.lst"
  [ -f "$m" ] || continue
  while IFS= read -r f; do
    git -C "$G" ls-files --error-unmatch "$f" > /dev/null 2>&1 && continue
    rm -f "$f"
  done < "$m"
  rm -f "$m"
done

: > "$W/deleted"
for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  echo "== lane $l"
  nc=0; nk=0; nd=0
  while IFS=$'\t' read -r st f; do
    case "$st" in
      D)
        git -C "$G" show "$BASE:$f" | tr -d '\r' > "$W/base"
        if [ ! -e "$f" ]; then
          nd=$((nd + 1))                 # an earlier lane deleted it too
        elif tr -d '\r' < "$f" | cmp -s - "$W/base"; then
          rm -f "$f"; nd=$((nd + 1)); printf '%s\n' "$f" >> "$W/deleted"
        else
          echo "   CONFLICT deletion of a file an earlier lane changed: $f"
          nk=$((nk + 1))
        fi
        continue;;
    esac
    echo "$f" | grep -qE "$excl" && continue
    case "$st" in
      M)
        if [ ! -e "$f" ]; then
          echo "   CONFLICT change to a file an earlier lane deleted: $f"
          nk=$((nk + 1)); continue
        fi
        git -C "$G" show "$BASE:$f" | tr -d '\r' > "$W/base"
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
      *) echo "   UNHANDLED $st $f";;
    esac
  done < "$W/$l.status"
  echo "   merged clean: $nc, with conflicts: $nk, deleted: $nd"
  lst="$W/$l.untracked"; : > "$lst"; clash=0; skip=0
  while IFS= read -r f; do
    if echo "$f" | grep -qE "$uskip"; then skip=$((skip + 1)); continue; fi
    if [ -e "$f" ]; then
      echo "   CLASH untracked: $f"; clash=$((clash + 1)); continue
    fi
    printf '%s\n' "$f" >> "$lst"
  done < <(git -C "$src" ls-files --others --exclude-standard)
  ( cd "$src" && tar -cf - -T "$lst" ) | tar -xpf -
  # untracked copies keep the lane's line endings; store them as LF
  while IFS= read -r f; do
    case "$f" in *.R|*.r|*.md|*.Rmd|*.sh|*.tsv|*.yaml|*.yml|*.txt|*.jobs|*.ps1|*.patch|*.sed|DESCRIPTION|NAMESPACE|.Rbuildignore)
      if grep -q $'\r' "$f"; then sed -i 's/\r$//' "$f"; fi;;
    esac
  done < "$lst"
  cp "$lst" "$LOG/copied-$l.lst"
  echo "   untracked copied: $(grep -c '' "$lst"), skipped: $skip," \
       "clashes: $clash"
done
# a retired package leaves its directory tree behind; remove what is
# empty, never a file
if [ -d extensions/frmtmb.spline ]; then
  find extensions/frmtmb.spline -depth -type d -empty -delete
  [ -d extensions/frmtmb.spline ] &&
    { echo "   LEFT in extensions/frmtmb.spline:";
      find extensions/frmtmb.spline -type f | sed 's/^/     /'; }
fi
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"
grep -ln '^<<<<<<< ' $(grep -vxFf "$W/deleted" "$W/files" |
  while read -r f; do [ -e "$f" ] && echo "$f"; done) 2>/dev/null |
  sed 's/^/   conflicted after the textual merge: /'
[ "${REL074_TEXTUAL_ONLY:-}" = true ] && exit 0

# R/links-brms.R from the merged generator; it reads brms from the user
# library
"$RS" dev/famlink-gen-links.R || exit 1
"$RS" dev/rel074-resolve.R || exit 1
"$RS" dev/rel074-news-merge.R || exit 1
left=$(for f in $(cat "$W/files") $(cat "$LOG"/copied-*.lst); do
  [ -f "$f" ] && grep -l '^<<<<<<< \|^>>>>>>> \|^=======$' "$f"; done)
[ -n "$left" ] && { echo "CONFLICT MARKERS LEFT:"; echo "$left"; exit 1; }
echo "MERGE, RESOLVE AND NEWS DONE"
