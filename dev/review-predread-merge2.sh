#!/usr/bin/env bash
# Item 8: trial three-way merge of predread, latres and nlfix over base d5d0056d,
# in a source copy outside every worktree. Plain files only: git merge-file works
# on copies and touches no repository; git show reads the base. Line endings are
# normalized first, because the worktrees differ in CRLF and a consolidation
# through git would normalize them.
set -u
R=/c/Users/adf44/source/r
P=$R/frmtmb-wt-predread
L=$R/frmtmb-wt-latres
N=$R/frmtmb-wt-nlfix
M=$R/rev-predread-src/merge
BASE=d5d0056d
rm -rf "$M"; mkdir -p "$M"
LOG=$P/dev/review-predread-out/r1/merge.txt
: > "$LOG"
files() { (cd "$1" && git status --short 2>/dev/null | sed 's/^...//' | grep -E '^(R/|NAMESPACE$|NEWS.md$|extensions/[^/]+/R/|extensions/[^/]+/NEWS.md$|tests/testthat/|extensions/frmtmb.sample/tests/|man/)'); }
files "$P" > "$M.p"; files "$L" > "$M.l"; files "$N" > "$M.n"
cat "$M.p" "$M.l" "$M.n" | sort -u > "$M.all"
echo "files changed by at least one lane: $(wc -l < "$M.all")" >> "$LOG"
cp -r "$P/R" "$P/DESCRIPTION" "$P/NAMESPACE" "$P/inst" "$P/data" "$P/man" "$M/"
mkdir -p "$M/extensions"
cp -r "$P/extensions/frmtmb.sample" "$M/extensions/"
nocr() { sed 's/\r$//' "$1"; }
while read -r f; do
  in_p=0; in_l=0; in_n=0
  grep -qxF "$f" "$M.p" && in_p=1
  grep -qxF "$f" "$M.l" && in_l=1
  grep -qxF "$f" "$M.n" && in_n=1
  count=$((in_p + in_l + in_n))
  mkdir -p "$M/$(dirname "$f")"
  base=$(mktemp); cur=$(mktemp); other=$(mktemp); merged=$(mktemp)
  (cd "$P" && git show "$BASE:$f" 2>/dev/null | sed 's/\r$//' > "$base")
  if [ $in_p = 1 ]; then nocr "$P/$f" > "$cur"
  elif [ $in_l = 1 ]; then nocr "$L/$f" > "$cur"
  else nocr "$N/$f" > "$cur"; fi
  c=0
  if [ $in_p = 1 ] && [ $in_l = 1 ]; then
    nocr "$L/$f" > "$other"
    git merge-file -p --diff3 -L predread -L base -L latres "$cur" "$base" "$other" > "$merged"; r=$?
    c=$((c + r)); cp "$merged" "$cur"
  fi
  if [ $in_n = 1 ] && { [ $in_p = 1 ] || [ $in_l = 1 ]; }; then
    nocr "$N/$f" > "$other"
    git merge-file -p --diff3 -L merged -L base -L nlfix "$cur" "$base" "$other" > "$merged"; r=$?
    c=$((c + r)); cp "$merged" "$cur"
  fi
  cp "$cur" "$M/$f"
  lanes=""; [ $in_p = 1 ] && lanes="$lanes predread"; [ $in_l = 1 ] && lanes="$lanes latres"; [ $in_n = 1 ] && lanes="$lanes nlfix"
  if [ "$count" -gt 1 ]; then echo "$f | lanes:$lanes | conflicts: $c" >> "$LOG"; fi
  rm -f "$base" "$cur" "$other" "$merged"
done < "$M.all"
cat "$LOG"
