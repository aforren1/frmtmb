#!/usr/bin/env bash
# Lane wt-arcovsample: the whole frmtmb.sample suite, ONE test file per R
# process, up to 8 at a time. Counts come from the RESULT lines the
# runner prints, which it emits from the test_file() result and not from
# the launch, so a file that aborts is visible.
#
#   bash dev/arcovsample-suite.sh
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-arcovsample
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$W/dev/arcovsample-log/suite"
mkdir -p "$OUT"
cd "$W"
files=$(ls extensions/frmtmb.sample/tests/testthat/test-*.R)
n=0
for f in $files; do
  b=$(basename "$f" .R)
  ( cd "$W"; "$R" dev/arcovsample-run.R lane frmtmb.sample "$f" \
      > "$OUT/$b.txt" 2>&1 ) &
  n=$((n + 1))
  if [ $((n % 8)) -eq 0 ]; then wait; fi
done
wait
echo "launched $n files"
grep -h "^RESULT" "$OUT"/*.txt | sort
echo "--- files with no RESULT line (aborted):"
for f in $files; do
  b=$(basename "$f" .R)
  grep -q "^RESULT" "$OUT/$b.txt" || echo "  $b"
done
echo "--- totals"
grep -h "^RESULT" "$OUT"/*.txt | sed 's/.*pass=//' | \
  awk -F'[= ]' '{p+=$1; f+=$3; e+=$5; s+=$7} END {
    printf "  files=%d pass=%d fail=%d err=%d skip=%d\n", NR, p, f, e, s}'
echo "--- BAD lines"
grep -h "^  BAD" "$OUT"/*.txt || echo "  (none)"
