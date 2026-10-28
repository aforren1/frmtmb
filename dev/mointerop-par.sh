#!/usr/bin/env bash
# One R process per test file, N at a time, through dev/mointerop-run-test.R
# (the lane's library, or the base with arm "base").
#   bash dev/mointerop-par.sh <tier> <job list> <N> <lane|base>
# The job list has one "<package> <test file>" line per file. The caller
# sets NOT_CRAN and the gates. Output: dev/mointerop-out/<tier>-files/ and
# dev/mointerop-out/<tier>.log with the RESULT lines.
set -u
tier="$1"; list="$2"; N="${3:-8}"; arm="${4:-lane}"
ROOT="/c/Users/adf44/source/r/frmtmb-wt-mointerop"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/mointerop-run-test.R"
OUT="$ROOT/dev/mointerop-out/$tier-files"
LOG="$ROOT/dev/mointerop-out/$tier.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$(dirname "$f")" && "$RS" "$RUN" "$arm" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"
  fi
done < "$list"
echo "$tier ran $ran of $want" >> "$LOG"
