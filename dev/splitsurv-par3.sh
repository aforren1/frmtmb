#!/usr/bin/env bash
# One R process per test file, N at a time, through
# dev/splitsurv-run-tests2.R; as dev/release/run-par.sh, with per-test
# times kept in a TSV beside each log.
#
#   bash dev/splitsurv-par.sh <out dir> <job list> [N]
#
# The caller sets NOT_CRAN and, for the lane's build, SPLITSURV_LIB.
set -u
out="$1"; list="$2"; N="${3:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="${SPLITSURV_RS:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
RUN="$ROOT/dev/splitsurv-run-tests2.R"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-C:/Users/adf44/Documents/.R/Makevars.win}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
: "${TMP:?TMP is unset}"
rm -rf "$out"; mkdir -p "$out"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$out/$p--$(basename "$f" .R)"
  ( "$RS" "$RUN" "$p" "$f" "$b.tsv" > "$b.txt" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  b="$out/$p--$(basename "$f" .R)"
  hit=$(grep -a "^RESULT" "$b.txt" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$p $hit"
  else echo "$p NO RESULT LINE for $f"; tail -3 "$b.txt"; fi
done < "$list" > "$out/summary.log"
echo "ran $ran of $want" >> "$out/summary.log"
tail -1 "$out/summary.log"
