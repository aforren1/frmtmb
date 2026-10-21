#!/usr/bin/env bash
# Merge the four lanes of the round opened 2026-10-20 onto the release
# tree, file by file with git merge-file, in the order loodraws, cdfs,
# nb2dm, fuzzfix. Every lane and the release tree are on 4854e425, so
# each file starts from its 4854e425 blob and each lane's change is
# merged in with the same blob as the base. It touches no index and
# makes no git write: blobs are read with `git show`.
#
#   bash dev/rel073-merge.sh
#
# Rerunnable, as dev/rel072-merge.sh: every run first restores each
# file any lane changed to its HEAD blob and removes the untracked lane
# files an earlier run copied ($LOG/copied-<lane>.lst), then merges the
# lanes in order and copies their untracked files again. The hand
# resolutions are applied by dev/rel073-resolve.R and the NEWS files
# are written by dev/rel073-news-merge.R, both run at the end.
#
# Line endings are normalized to LF first, as git stores them: the lanes
# left LF and CRLF working copies, and merging those against LF blobs
# gives a conflict on every line. Generated files are left out, because
# their generators run after the merge: man/*.Rd and NAMESPACE
# (roxygen), frmtmb.sample's copy of helper-brms-suite.R, and
# R/links-brms.R, which dev/famlink-gen-links.R writes from brms (run
# here, after the merge, since lane nb2dm changed the generator).
# Untracked lane files are copied, except generated Rd files, Stan
# caches, the brms suite and OpenBLAS builds of R (gitignored anyway;
# named so that no change to the ignore file lets them in).
set -u
BASE=4854e425
REL=${REL073_ROOT:-C:/Users/adf44/source/r/frmtmb-wt-release}
# the checkout whose blobs are read: the release tree itself, or it for
# a scratch copy made with git archive (REL073_ROOT), which has no git
G=${REL073_GIT:-$REL}
export REL073_ROOT="$REL" REL073_GIT="$G"
LOG=$REL/dev/rel073-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LANES="loodraws cdfs nb2dm fuzzfix"
export REL073_LANES="$LANES"
cd "$REL" || exit 1
[ "$(git -C "$G" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "release tree is not on $BASE"; exit 1; }
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|R/links-brms\.R|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
uskip='^(man/.*\.Rd|extensions/[^/]*/man/.*\.Rd|.*stan-cache.*|dev/ciharden-out/.*|dev/brms-suite/.*)$'
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
echo "files changed by some lane (without generated ones): $(wc -l < "$W/files")"

# restore every file any lane changed, generated ones too, so that a
# rerun starts from the base
while IFS= read -r f; do
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

for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  echo "== lane $l"
  nc=0; nk=0
  while IFS=$'\t' read -r st f; do
    echo "$f" | grep -qE "$excl" && continue
    case "$st" in
      M)
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
      D) echo "   UNHANDLED deletion $f";;
      *) echo "   UNHANDLED $st $f";;
    esac
  done < "$W/$l.status"
  echo "   merged clean: $nc, with conflicts: $nk"
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
    case "$f" in *.R|*.r|*.md|*.sh|*.tsv|*.Rmd|*.yaml|*.yml|*.txt|*.jobs|*.ps1|*.patch)
      if grep -q $'\r' "$f"; then sed -i 's/\r$//' "$f"; fi;;
    esac
  done < "$lst"
  cp "$lst" "$LOG/copied-$l.lst"
  echo "   untracked copied: $(grep -c '' "$lst"), skipped: $skip," \
       "clashes: $clash"
done
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"
grep -ln '^<<<<<<< ' $(cat "$W/files") 2>/dev/null |
  sed 's/^/   conflicted after the textual merge: /'
[ "${REL073_TEXTUAL_ONLY:-}" = true ] && exit 0

# R/links-brms.R from the merged generator (lane nb2dm added two
# families to it); it reads brms from the user library
"$RS" dev/famlink-gen-links.R || exit 1
"$RS" dev/rel073-resolve.R || exit 1
"$RS" dev/rel073-news-merge.R || exit 1
grep -ln '^<<<<<<< \|^>>>>>>> \|^=======$' $(cat "$W/files") && {
  echo "CONFLICT MARKERS LEFT"; exit 1; }
echo "MERGE, RESOLVE AND NEWS DONE"
