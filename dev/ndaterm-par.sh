#!/usr/bin/env bash
# Lane ndaterm: every test file of the listed packages, ONE R process per
# file, N at a time, through dev/ndaterm-runfile.R on the lane library.
#   bash dev/ndaterm-par.sh <tag> <N> <pkg>...
# pkg is frmtmb or an extension name. Per-file logs in
# dev/ndaterm-suite-<tag>/, the RESULT lines in
# dev/ndaterm-suite-<tag>/RESULTS.txt. Extra env (FRMTMB_BRMS_FIT_TESTS)
# is the caller's.
set -u
tag="$1"; N="$2"; shift 2
ROOT=/c/Users/adf44/source/r/frmtmb-wt-ndaterm
RS="${RS:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
OUT="$ROOT/dev/ndaterm-suite-$tag"
mkdir -p "$OUT"
list="$OUT/jobs.txt"
if [ -n "${JOBS:-}" ]; then cp "$JOBS" "$list"; set --; else : > "$list"; fi
for p in "$@"; do
  if [ "$p" = frmtmb ]; then d="$ROOT/tests/testthat"
  else d="$ROOT/extensions/$p/tests/testthat"; fi
  for f in "$d"/test-*.R; do echo "$p $f" >> "$list"; done
done
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$ROOT/dev/ndaterm-runfile.R" "$p" "$f" mine > "$o" 2>&1 ) &
done < "$list"
wait
res="$OUT/RESULTS.txt"; : > "$res"
while read -r p f; do
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$res"
  else echo "$p NO RESULT $(basename "$f")" >> "$res"; fi
done < "$list"
echo "ran $(grep -c RESULT "$res") of $(wc -l < "$list")" >> "$res"
