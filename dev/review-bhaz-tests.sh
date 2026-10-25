#!/usr/bin/env bash
# Run one queue of test files, one file per R process, against the lane
# library. Usage: review-bhaz-tests.sh <queue-name> <pkg>:<file> ...
WT=/c/Users/adf44/source/r/frmtmb-wt-bhaz
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$WT/dev/review-bhaz-log
mkdir -p "$OUT"
q=$1; shift
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
for spec in "$@"; do
  pkg=${spec%%:*}; f=${spec#*:}
  log=$OUT/t-$(basename "$f" .R)-$pkg.txt
  "$RS" "$WT/dev/bhaz-run-test.R" "$pkg" "$WT/$f" lane > "$log" 2>&1
  echo "$q $(grep -h '^lib:' "$log" | head -1) :: $(grep -h '^RESULT' "$log")" >> "$OUT/summary-$q.txt"
done
echo "$q DONE" >> "$OUT/summary-$q.txt"
