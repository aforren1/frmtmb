#!/usr/bin/env bash
# Merge each file that two sources changed in the reverse lane order
# (seflags, optrap, moparity) onto its HEAD blob and compare with the
# release tree: an identical result means no hunk of one lane depends
# on where another lane's landed. R/fit.R has the mo_search() conflict
# in either order, so it is compared outside mo_search() and roxygen
# (the resolve's two documentation edits), and the reverse merge's
# conflict block is printed for reading.
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/rel070-log/revorder
cd "$REL"; rm -rf "$W"; mkdir -p "$W"
cut_mo() {
  awk '/^mo_search <- function/{skip=1} skip&&/^}$/{skip=0; next} !skip' "$1" |
    grep -v "^#'"
}
for f in R/confint.R R/predict.R R/se-check.R R/fit.R \
         tests/testthat/test-se-check.R; do
  m="$W/$(echo "$f" | tr '/' '_')"
  git show "HEAD:$f" | tr -d '\r' > "$m"
  for l in seflags optrap moparity; do
    src=/c/Users/adf44/source/r/frmtmb-wt-$l
    git -C "$src" diff --quiet 48c1e6c9 -- "$f" 2>/dev/null && continue
    git show "48c1e6c9:$f" | tr -d '\r' > "$W/base"
    tr -d '\r' < "$src/$f" > "$W/theirs"
    git merge-file -L release -L base -L "$l" "$m" "$W/base" "$W/theirs"
    echo "$f $l conflicts=$?"
  done
  tr -d '\r' < "$f" > "$m.rel"
  if [ "$f" = R/fit.R ]; then
    if diff -q <(cut_mo "$m") <(cut_mo "$m.rel") > /dev/null; then
      echo "SAME outside mo_search() and roxygen: $f"
    else echo "DIFFERENT: $f"; fi
    awk '/^mo_search <- function/{p=1} p' "$m" | sed -n '/^<<<<<<</,/^>>>>>>>/p'
  elif cmp -s "$m" "$m.rel"; then echo "SAME: $f"
  else echo "DIFFERENT: $f"; fi
done
