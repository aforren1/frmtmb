#!/usr/bin/env bash
# One R process per test file, N at a time, through dev/cifix-run1.R.
#   bash dev/cifix-pardrv.sh <run name> <ref|ob> <lib|base> <job list> [N]
# The job list has "<package> <test file>" lines. Per-file logs go to
# dev/cifix-log/<run name>-suite/; the summary to
# dev/cifix-log/<run name>.sum with one RESULT line per file, or NO
# RESULT LINE and the tail of the log.
set -u
run="$1"; kind="$2"; lib="$3"; list="$4"; N="${5:-12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-cifix
cd "$ROOT"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP is unset}"
if [ "$kind" = ob ]; then
  RS="$ROOT/dev/cifix-out/Rob/bin/x64/Rscript.exe"
  export OPENBLAS_NUM_THREADS=4
else
  RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
fi
OUT="$ROOT/dev/cifix-log/$run-suite"
SUM="$ROOT/dev/cifix-log/$run.sum"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$SUM"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/cifix-run1.R "$lib" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
while read -r p f; do
  [ -n "$p" ] || continue
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$SUM"; else
    echo "$p NO RESULT LINE for $f" >> "$SUM"; tail -3 "$o" >> "$SUM"; fi
done < "$list"
echo "ran $(grep -c ' RESULT ' "$SUM") of $(grep -c . "$list")" >> "$SUM"
