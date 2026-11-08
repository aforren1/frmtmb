#!/usr/bin/env bash
# Merge the twelve lane trees of the round of 2026-11-06 into the release
# tree, file by file, three-way against 3cd2a0c8. Run once, from a clean
# release tree. It writes only inside the release tree and reads the lanes.
#
#   bash dev/rel080-merge.sh
#
# A lane's whole tree is taken, as its commit will hold it: every modified
# or untracked file that is not ignored, except machine output that
# dev/.gitignore lets through (reviewers' output folders, a lane's copy of
# the equivalence survey's output, a reviewer's trial-merge folder).
# Text is merged LF-normalized and written back CRLF, as core.autocrlf
# checks it out. R/gp-start.R is merged with diff3 -E, fitedge first, so
# that the written resolutions (dev/rel080-resolve.R) apply as reviewed.
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel080
BASE=3cd2a0c8
ORDER="ci latperf mfxse loglik osa nestquad arstep fitedge bigx nlstart readers superassign"
OUT=$REL/dev/rel080-log/merge
SKIP='^dev/(review-[^/]*-out[^/]*/|[^/]*-equiv-out[^/]*/|review-[^/]*-merge/|[^/]*\.progress$|[^/]*\.done$)'
UNION='^(NEWS\.md|extensions/[^/]*/NEWS\.md|dev/test-backlog\.md|dev/upstream-bugs\.md)$'
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
[ "$(git -C "$REL" rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "release tree not on $BASE"; exit 1; }
[ -z "$(git -C "$REL" status --porcelain --untracked-files=no)" ] ||
  { echo "release tree has tracked changes; merge only once"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT/tmp"
: > "$OUT/conflicts.txt"; : > "$OUT/taken.txt"; : > "$OUT/skipped.txt"
lf() { tr -d '\r' < "$1"; }
# a last line with no newline keeps none (and gets no lone CR)
crlf() { sed 's/$/\r/' | perl -0pe 's/\r\z//'; }
is_text() {  # text as git sees it at the base (binary files are -text)
  local a
  a=$(cd "$REL" && git ls-files --eol -- "$1" | awk '{print $1}')
  [ -z "$a" ] || [ "$a" = "i/lf" ] || [ "$a" = "i/none" ]
}
for l in $ORDER; do
  T=$R/frmtmb-wt-$l
  [ "$(git -C "$T" rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "$l not on $BASE"; exit 1; }
  (cd "$T" && git status --porcelain=v1 -z --untracked-files=all --no-renames) |
    tr '\0' '\n' > "$OUT/tmp/list-$l"
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
        git -C "$REL" show "$BASE:$p" > "$b"
        lf "$REL/$p" > "$o"; lf "$T/$p" > "$t"
        if [[ "$p" =~ $UNION ]]; then
          # every lane adds its own top section: keep all, edit after
          # (dev/rel080-news.R)
          cp "$o" "$OUT/tmp/m"
          git merge-file --union -L rel -L base -L "$l" "$OUT/tmp/m" "$b" "$t"; rc=0
        elif [ "$p" = "R/gp-start.R" ] && [ "$l" = fitedge ]; then
          # fitedge first, as dev/rel080-resolve.R expects
          set +e; diff3 -m -E "$t" "$b" "$o" > "$OUT/tmp/m"; rc=$?; set -e
          [ "$rc" -le 1 ] || { echo "diff3 failed on $p"; exit 1; }
        else
          cp "$o" "$OUT/tmp/m"
          set +e; git merge-file -L rel -L base -L "$l" "$OUT/tmp/m" "$b" "$t"; rc=$?; set -e
          [ "$rc" -ge 0 ] || { echo "merge-file failed on $p"; exit 1; }
        fi
        crlf < "$OUT/tmp/m" > "$REL/$p"
        if [ "$rc" -ne 0 ]; then
          k=$(grep -c '^<<<<<<< ' "$OUT/tmp/m" || true)
          echo "$l CONFLICT $p ($k)" | tee -a "$OUT/conflicts.txt"; c=$((c+1))
        fi ;;
      *) echo "$l UNHANDLED '$st' $p" | tee -a "$OUT/conflicts.txt"; c=$((c+1)) ;;
    esac
    echo "$l $st $p" >> "$OUT/taken.txt"; n=$((n+1))
  done < "$OUT/tmp/list-$l"
  echo "MERGED $l: $n paths, $c needing a hand"
  if [ "$c" -gt 0 ]; then
    "$RS" "$REL/dev/rel080-resolve.R" "$l" | tee -a "$OUT/resolved.txt"
  fi
  left=$(grep -rlE '^(<<<<<<<|>>>>>>>) ' "$REL/R" "$REL/tests" "$REL/man" \
           "$REL/extensions" "$REL/NEWS.md" "$REL/dev/test-backlog.md" || true)
  [ -z "$left" ] || { echo "markers left after $l:"; echo "$left"; exit 1; }
done
echo "MERGE DONE"
