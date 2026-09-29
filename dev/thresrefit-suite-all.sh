#!/usr/bin/env bash
# The whole core suite, ONE test file per R process, $JOBS at a time.
# Logs in dev/thresrefit-testlog-all/. Run ONCE, on the final pass.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-thresrefit
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/thresrefit-testlog-all"
JOBS=${JOBS:-7}
mkdir -p "$OUT"
cd "$ROOT/tests/testthat" || exit 1
FILES=$(ls test-*.R)
cd "$ROOT" || exit 1

n=0
for f in $FILES; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  n=$((n + 1))
  ( "$RS" dev/thresrefit-runtest.R "$f" > "$OUT/$f.log" 2>&1 ) &
done
wait
echo "launched $n files"
grep -h "^RESULT" "$OUT"/*.log | sort
echo "--- files with no RESULT line (aborted) ---"
miss=0
for f in $FILES; do
  if ! grep -q "^RESULT" "$OUT/$f.log" 2>/dev/null; then
    echo "$f"; miss=$((miss + 1))
  fi
done
echo "aborted: $miss"
echo "--- totals, counted from RESULT lines only ---"
grep -h "^RESULT" "$OUT"/*.log |
  sed 's/.*pass=\([0-9]*\) fail=\([0-9]*\) err=\([0-9]*\) skip=\([0-9]*\)/\1 \2 \3 \4/' |
  awk '{p+=$1; f+=$2; e+=$3; s+=$4; n++}
       END {print "files="n" pass="p" fail="f" err="e" skip="s}'
echo "--- files with a nonzero fail or err ---"
grep -h "^RESULT" "$OUT"/*.log | grep -v "fail=0 err=0" || echo "none"
echo "--- files with a nonzero skip ---"
grep -h "^RESULT" "$OUT"/*.log | grep -v "skip=0" || echo "none"
