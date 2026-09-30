#!/usr/bin/env bash
# Reviewer, lane ordinal: one R process per test file, N at a time.
#   bash dev/ordinal-rev-run-par.sh <arm> <tier> <job list> [N]
# Job list: "<package> <test file path>" per line.
# Output: dev/ordinal-rev-suite-<tier>/, with <tier>.log of RESULT lines.
set -u
arm="$1"; tier="$2"; list="$3"; N="${4:-20}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordinal"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/ordinal-rev-runtest.R"
OUT="$ROOT/dev/ordinal-rev-suite-$tier"
LOG="$OUT/$tier.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" "$arm" > "$o" 2>&1 ) &
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
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
