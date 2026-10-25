#!/usr/bin/env bash
# Review short check: four files one per process, then core reverse session.
REL=/c/Users/adf44/source/r/frmtmb-wt-rel075
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$REL/dev/review-rel075-log/fix
mkdir -p "$OUT"
export NOT_CRAN=true
for j in "frmtmb tests/testthat/test-extension-home-mixed.R" "frmtmb tests/testthat/test-autocor-seam.R" "frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-sar.R" "frmtmb tests/testthat/test-guardrail1.R"; do
  set -- $j
  ( "$RS" "$REL/dev/ciharden-run1.R" base "$1" "$REL/$2" > "$OUT/$1--$(basename "$2" .R).log" 2>&1 ) &
done
"$RS" "$REL/dev/release/run-session.R" frmtmb reverse > "$OUT/session-core-reverse.log" 2>&1 &
wait
for f in "$OUT"/*--*.log; do echo "$(grep '^RESULT' "$f") | $(grep '^lib:' "$f")"; done > "$OUT/summary.txt"
grep '^SESSION' "$OUT/session-core-reverse.log" >> "$OUT/summary.txt"
echo ALLDONE >> "$OUT/summary.txt"
