#!/usr/bin/env bash
# Run a list of test files, one R process per file, N at a time, through
# dev/release/run-tests.R, and assemble one tier log in list order.
#
#   bash dev/release/run-par.sh <tier> <job list> [N]
#
# The job list has one "<package> <test file>" line per file. The caller
# sets the tier's environment (NOT_CRAN, FRMTMB_BRMS_FIT_TESTS, ...).
# Output: dev/release/<tier>-files/<package>--<file>.txt per file and
# dev/release/<tier>.log with the RESULT lines, grouped by package as
# run-suite.ps1 writes them, then "<TIER> ran X of Y".
#
# Why a second driver beside the PowerShell ones: those run one file at
# a time, and on the 63 GB machine of 2026-09-28 the whole suite at one
# process per file takes hours that parallel processes cut to minutes.
# Each process still runs ONE file (dev/lane-rules.md). The per-file
# logs are written by the R process's own redirection, never appended
# by a watcher, so the tail -f trap of lane-rules.md cannot reach them;
# still, do not tail them while a run is live.
set -u
tier="$1"; list="$2"; N="${3:-16}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/release/$tier-files"
LOG="$ROOT/dev/release/$tier.log"
[ -f "$RUN" ] || { echo "runner missing: $RUN"; exit 1; }
[ -s "$list" ] || { echo "empty job list: $list"; exit 1; }
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
MK=C:/Users/adf44/Documents/.R/Makevars.win
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-$MK}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0; prev=""
while read -r p f; do
  [ -n "$p" ] || continue
  if [ "$p" != "$prev" ]; then
    [ -n "$prev" ] && echo "== $prev END ==" >> "$LOG"
    echo "== $p $(grep -c "^$p " "$list") files ==" >> "$LOG"
    prev="$p"
  fi
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$hit" >> "$LOG"
  else
    echo "NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
[ -n "$prev" ] && echo "== $prev END ==" >> "$LOG"
up=$(echo "$tier" | tr '[:lower:]' '[:upper:]')
echo "$up ran $ran of $want" >> "$LOG"
echo "$up ran $ran of $want"
