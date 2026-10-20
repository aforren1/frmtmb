#!/usr/bin/env bash
# Lane sampleprior: run a job list of "<arm> <package> <gated> <file>"
# lines, one R process per file, N at a time, through
# dev/sampleprior-runtest.R. Per-file logs go to dev/sampleprior-suite/,
# which the R process writes itself (no watcher tails them).
#   bash dev/sampleprior-par.sh <job list> <out dir name> [N]
set -u
list="$1"; name="$2"; N="${3:-8}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-sampleprior
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/sampleprior-suite/$name"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
export TEMP="$TMP"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r arm p gated f; do
  [ -n "$arm" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( FRMTMB_BRMS_FIT_TESTS="$gated" "$RS" "$ROOT/dev/sampleprior-runtest.R" \
      "$arm" "$p" "$ROOT/$f" > "$o" 2>&1 ) &
done < "$list"
wait
for o in "$OUT"/*.txt; do
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  lib=$(grep -a "^lib:" "$o" | head -1)
  echo "$(basename "$o" .txt) | ${hit:-NO RESULT} | $lib"
done > "$OUT.summary"
echo "done: $(grep -c RESULT "$OUT.summary") of $(wc -l < "$OUT.summary")"
