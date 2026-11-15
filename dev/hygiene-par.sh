#!/usr/bin/env bash
# dev/release/run-par.sh for lane hygiene: same contract, but runs
# dev/hygiene-run-tests.R (HYG_ARM picks the library) and writes to
# dev/hygiene-out/<tier>-files/ and dev/hygiene-out/<tier>.log.
#
#   HYG_ARM=mine|base bash dev/hygiene-par.sh <tier> <job list> [N]
set -u
tier="$1"; list="$2"; N="${3:-3}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/hygiene-run-tests.R"
OUT="$ROOT/dev/hygiene-out/$tier-files"
LOG="$ROOT/dev/hygiene-out/$tier.log"
[ -f "$RUN" ] || { echo "runner missing: $RUN"; exit 1; }
[ -s "$list" ] || { echo "empty job list: $list"; exit 1; }
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
MK=C:/Users/adf44/Documents/.R/Makevars.win
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-$MK}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
export NOT_CRAN="${NOT_CRAN:-true}"
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
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "TIER $tier arm=${HYG_ARM:-mine} ran $ran of $want" >> "$LOG"
cat "$LOG"
