#!/usr/bin/env bash
# Review of lane priors2, item 8: for every file priors2 changed, the
# base-line ranges of its hunks against every other r083 lane's hunks on
# the same file. Prints a pair when the ranges overlap or sit within 3
# lines (a three-way merge conflicts there). Read-only: git show to a
# scratch dir and diff, no git diff.
set -u
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
P=$R/frmtmb-wt-priors2
TMP=${TMPDIR:-/tmp}/priors2-overlap.$$
mkdir -p "$TMP"
lanes="ci082"
files="$(cd $P && cat "$1")"
ranges() { # base file, new file -> "start end" per hunk, in base lines
  diff -U0 --strip-trailing-cr "$1" "$2" | sed -n 's/^@@ -\([0-9]*\)\(,\([0-9]*\)\)\? .*/\1 \3/p' |
    awk '{s=$1; n=($2=="")?1:$2; e=(n==0)?s:s+n-1; print s, e}'
}
for f in $files; do
  mkdir -p "$TMP/base/$(dirname $f)"
  (cd $P && git show ec8178d3:$f > "$TMP/base/$f" 2>/dev/null) || continue
  ranges "$TMP/base/$f" "$P/$f" > "$TMP/p.txt"
  for l in $lanes; do
    W=$R/frmtmb-wt-$l
    [ -f "$W/$f" ] || continue
    cmp -s <(tr -d '\r' < "$W/$f") <(tr -d '\r' < "$TMP/base/$f") && continue
    ranges "$TMP/base/$f" "$W/$f" > "$TMP/o.txt"
    while read ps pe; do
      while read os oe; do
        if [ $((os - 3)) -le "$pe" ] && [ $((oe + 3)) -ge "$ps" ]; then
          echo "$f  priors2 $ps-$pe  $l $os-$oe"
        fi
      done < "$TMP/o.txt"
    done < "$TMP/p.txt"
    echo "touched: $f by $l ($(wc -l < "$TMP/o.txt") hunks)" >> "$TMP/touched.txt"
  done
done
echo "--- files priors2 changed that other lanes also changed"
cat "$TMP/touched.txt" 2>/dev/null
rm -rf "$TMP"
