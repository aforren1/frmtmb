#!/usr/bin/env bash
# Merge the fourteen lane trees of the round of 2026-11-13 into the release
# tree, file by file, three-way: base = the file at ec8178d3 (the lanes'
# base), current = the release tree (8159568e, frmtmb 0.82.1), other = the
# lane's file. Run once, from a clean release tree. It writes only inside
# the release tree and reads the lanes. Adapted from dev/rel082-merge.sh.
#
#   bash dev/rel083-merge.sh
#
# A lane's paths are read by dev/rel083-list.sh (git ls-files plus a
# content compare), never by git status or git diff, which refresh a
# lane's index. Every path is taken but machine output that dev/.gitignore
# lets through: reviewers' output folders, a lane's equivalence-survey
# copy, trial-merge folders, *.progress and *.done, and this round's two:
# a Python byte-code cache a reviewer's script left in .github/scripts/
# and lane hygiene's reviewer copy of the base codemeta.json.
# Order: selimit2, quadseam2, slopescale first (all three edit
# R/se-check.R; the reviewers merged them in this order with 0
# conflicts), then the brief's order.
# Text is merged LF-normalized and written back with the line endings the
# release checkout has for that file. NEWS files, dev/test-backlog.md and
# dev/upstream-bugs.md are merged with --union (every lane adds its own
# section) and rebuilt after.
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel083
BASE=ec8178d3
RELBASE=8159568e
ORDER="selimit2 quadseam2 slopescale nlstop2 priors2 famend silent2 arspat2 refits parity2 wcheck2 gpfield predraw2 hygiene"
OUT=$REL/dev/rel083-log/merge
SKIP='^(dev/(review-[^/]*-out[^/]*/|[^/]*-equiv-out[^/]*/|review-[^/]*-merge/|compat-m-race[^/]*/|[^/]*\.progress$|[^/]*\.done$|reviews/hygiene-review-cm-base\.json$)|\.github/scripts/__pycache__/)'
UNION='^(NEWS\.md|extensions/[^/]*/NEWS\.md|dev/test-backlog\.md|dev/upstream-bugs\.md)$'
[ "$(git -C "$REL" rev-parse --short=8 HEAD)" = "$RELBASE" ] || { echo "release tree not on $RELBASE"; exit 1; }
pre=$(bash "$REL/dev/rel083-list.sh" "$REL" | grep -v '^?? ' || true)
[ -z "$pre" ] || { echo "release tree has tracked changes; merge only once"; echo "$pre"; exit 1; }
mkdir -p "$OUT/tmp"
: > "$OUT/conflicts.txt"; : > "$OUT/taken.txt"; : > "$OUT/skipped.txt"
lf() { tr -d '\r' < "$1"; }
# a last line with no newline keeps none (and gets no lone CR)
crlf() { sed 's/$/\r/' | perl -0pe 's/\r\z//'; }
is_text() {  # text as git sees it (binary files are -text)
  local a
  a=$(cd "$REL" && git ls-files --eol -- "$1" | awk '{print $1}')
  [ -z "$a" ] || [ "$a" = "i/lf" ] || [ "$a" = "i/none" ]
}
for l in $ORDER; do
  T=$R/frmtmb-wt-$l
  [ "$(git -C "$T" rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "$l not on $BASE"; exit 1; }
  bash "$REL/dev/rel083-list.sh" "$T" > "$OUT/tmp/list-$l"
  n=0; c=0
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    st=${line:0:2}; p=${line:3}
    if [[ "$p" =~ $SKIP ]]; then echo "$l $p" >> "$OUT/skipped.txt"; continue; fi
    case "$st" in
      "??")
        if [ -e "$REL/$p" ]; then
          cmp -s "$T/$p" "$REL/$p" ||
            { echo "$l ADDED-DIFFERS $p" | tee -a "$OUT/conflicts.txt"; c=$((c+1)); }
        else
          mkdir -p "$(dirname "$REL/$p")"; cp -p "$T/$p" "$REL/$p"
        fi ;;
      " M")
        if ! is_text "$p"; then
          if cmp -s <(git -C "$REL" show "$BASE:$p") "$REL/$p"; then cp -p "$T/$p" "$REL/$p"
          else echo "$l BINARY-BOTH $p" | tee -a "$OUT/conflicts.txt"; c=$((c+1)); fi
          echo "$l $st $p" >> "$OUT/taken.txt"; n=$((n+1)); continue
        fi
        b=$OUT/tmp/base; o=$OUT/tmp/ours; t=$OUT/tmp/theirs
        git -C "$REL" show "$BASE:$p" | tr -d '\r' > "$b"
        if [ "$(tr -dc '\r' < "$REL/$p" | head -c 1 | wc -c)" -gt 0 ]; then eol=crlf; else eol=lf; fi
        lf "$REL/$p" > "$o"; lf "$T/$p" > "$t"
        cp "$o" "$OUT/tmp/m"
        if [[ "$p" =~ $UNION ]]; then
          git merge-file --union -L rel -L base -L "$l" "$OUT/tmp/m" "$b" "$t"; rc=0
        else
          set +e; git merge-file -L rel -L base -L "$l" "$OUT/tmp/m" "$b" "$t"; rc=$?; set -e
          [ "$rc" -ge 0 ] || { echo "merge-file failed on $p"; exit 1; }
        fi
        if [ "$eol" = crlf ]; then crlf < "$OUT/tmp/m" > "$REL/$p"; else cp "$OUT/tmp/m" "$REL/$p"; fi
        if [ "$rc" -ne 0 ]; then
          k=$(grep -c '^<<<<<<< ' "$OUT/tmp/m" || true)
          echo "$l CONFLICT $p ($k)" | tee -a "$OUT/conflicts.txt"; c=$((c+1))
        fi ;;
      *) echo "$l UNHANDLED '$st' $p" | tee -a "$OUT/conflicts.txt"; c=$((c+1)) ;;
    esac
    echo "$l $st $p" >> "$OUT/taken.txt"; n=$((n+1))
  done < "$OUT/tmp/list-$l"
  echo "MERGED $l: $n paths, $c needing a hand"
done
left=$(grep -rlE '^(<<<<<<<|>>>>>>>) ' "$REL/R" "$REL/tests" "$REL/man" \
         "$REL/extensions" "$REL/NEWS.md" "$REL/dev" --include='*.R' --include='*.md' \
         --include='*.Rd' --include='*.tsv' --include='*.sh' 2>/dev/null || true)
[ -z "$left" ] || { echo "markers left:"; echo "$left"; }
echo "MERGE DONE"
