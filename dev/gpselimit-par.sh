#!/usr/bin/env bash
# Lane gpselimit: run a job list ("<package> <test file>" per line) one
# file per R process, N at a time (default 2), through
# dev/release/run-tests.R with the lane's library first (FRMTMB_RELLIB),
# or the base build alone with BASE=1. Per-file logs go to
# dev/gpselimit-out/<name>/, the RESULT lines to dev/gpselimit-out/<name>.log.
# Usage: bash dev/gpselimit-par.sh <name> <job list> [N]
set -u
name="$1"; list="$2"; N="${3:-2}"
ROOT=C:/Users/adf44/source/r/frmtmb-wt-gpselimit
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/gpselimit-out/$name"
LOG="$ROOT/dev/gpselimit-out/$name.log"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
if [ "${BASE:-0}" = "1" ]; then
  unset FRMTMB_RELLIB
else
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpselimit-lib
fi
: "${TMP:?TMP is unset}"
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
  lib=$(grep -a "^lib:" "$o" | head -1)
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit | $lib" >> "$LOG"
  else
    echo "NO RESULT LINE for $p $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "$name ran $ran of $want" >> "$LOG"
