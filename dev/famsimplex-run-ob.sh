#!/usr/bin/env bash
# Run test files one R process per file, N at a time, through
# dev/release/run-tests.R, with the lane library ahead of rellib-r8
# (FRMTMB_RELLIB), or rellib-r8 alone for a "before" arm.
#
#   bash dev/famsimplex-run.sh <tag> <job list> [N] [before]
#
# Job list lines: "<package> <test file>". Output: one file per test in
# dev/famsimplex-tests/<tag>/ and the RESULT lines in
# dev/famsimplex-tests/<tag>.log. Never tail these while a run is live.
set -u
tag="$1"; list="$2"; N="${3:-8}"; arm="${4:-after}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-famsimplex
RS="${FAMSIMPLEX_RS:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/famsimplex-tests/$tag"
LOG="$ROOT/dev/famsimplex-tests/$tag.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
if [ "$arm" = "before" ]; then unset FRMTMB_RELLIB; else
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-famsimplex-lib; fi
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
  lib=$(grep -a "^lib: " "$o" | tail -1)
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$hit | $lib" >> "$LOG"
  else echo "NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"; fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
