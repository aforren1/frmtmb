#!/usr/bin/env bash
# Reviewer driver: one R process per test file, N at a time.
#   bash formrobust-rev-run-par.sh <tier> <job list> [N]
# Per-file output: dev/formrobust-rev-log/<tier>-files/<pkg>--<file>.txt
# Summary: dev/formrobust-rev-log/<tier>.log
set -u
tier="$1"; list="$2"; N="${3:-12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/formrobust-rev-run-tests.R"
OUT="$ROOT/dev/formrobust-rev-log/$tier-files"
LOG="$ROOT/dev/formrobust-rev-log/$tier.log"
[ -s "$list" ] || { echo "empty job list: $list"; exit 1; }
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/formrobust-rev-stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$(dirname "$f")" && "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else echo "$p NORESULT $f" >> "$LOG"; fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
