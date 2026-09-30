#!/bin/sh
# Lane wt-defects: a package's whole test suite, ONE FILE PER R PROCESS,
# in $P parallel lanes, against the lane library (DEFECTS_ARM=after) or
# the base build (DEFECTS_ARM=before). Gated tiers run too when
# FRMTMB_BRMS_FIT_TESTS=true is exported by the caller.
#   DEFECTS_ARM=after P=12 sh dev/defects-suite.sh <pkg> <tag>
# Output: dev/defects-log/suite-<tag>/<file>.txt, and a summary of every
# RESULT line in dev/defects-log/suite-<tag>.txt.
#
# Rscript is called by its full path and nothing re-enters a shell by
# name: with Rtools first on PATH, `sh` and `xargs` resolve to Rtools'
# own, which drop TMP (dev/lane-rules.md), and a first spelling of this
# script ran 0 of 183 files that way.
cd "$(dirname "$0")/.." || exit 1
pkg=$1; tag=$2; P=${P:-12}
case "$pkg" in
  frmtmb) dir=tests/testthat ;;
  *) dir=extensions/$pkg/tests/testthat ;;
esac
out=dev/defects-log/suite-$tag
mkdir -p "$out"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export DEFECTS_ARM="${DEFECTS_ARM:-after}"
files=$(ls $dir/test-*.R)
i=0
for k in $(seq 1 "$P"); do
  (
    j=0
    for f in $files; do
      j=$((j + 1))
      [ $(( (j - 1) % P + 1 )) -eq "$k" ] || continue
      b=$(basename "$f" .R)
      "$RS" dev/defects-runtest.R "$pkg" "$f" > "$out/$b.txt" 2>&1
    done
  ) &
done
wait
n=$(echo "$files" | wc -l)
grep -h "^RESULT" $out/*.txt > dev/defects-log/suite-$tag.txt
got=$(wc -l < dev/defects-log/suite-$tag.txt)
echo "SUITE $pkg ran $got of $n files" >> dev/defects-log/suite-$tag.txt
for f in $out/*.txt; do
  grep -q "^RESULT" "$f" || echo "NO RESULT $(basename "$f")" >> dev/defects-log/suite-$tag.txt
done
