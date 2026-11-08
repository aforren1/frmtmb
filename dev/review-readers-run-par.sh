#!/usr/bin/env bash
# Review of lane readers: run a list of test files, one R process per file,
# N at a time (at most 3), through the lane's dev/readers-run-tests.R, which
# puts FRMTMB_RELLIB ahead of rellib-r18. Logs go to
# dev/review-readers-out/<tier>-files/ and dev/review-readers-out/<tier>.log.
#
#   bash dev/review-readers-run-par.sh <tier> <job list> [N] [lib]
set -u
tier="$1"; list="$2"; N="${3:-3}"
[ "$N" -le 3 ] || N=3
export FRMTMB_RELLIB="${4:-C:/Users/adf44/source/r/wt-readers-lib}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-readers"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/readers-run-tests.R"
OUT="$ROOT/dev/review-readers-out/$tier-files"
LOG="$ROOT/dev/review-readers-out/$tier.log"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
cd "$ROOT"
i=0
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  i=$((i + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( R_PARALLEL_PORT=$((21700 + i)) "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
while read -r p f; do
  [ -n "$p" ] || continue
  o="$OUT/$p--$(basename "$f" .R).txt"
  echo "$(grep -a '^lib:' "$o" | head -1)" >> "$LOG"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then echo "$hit" >> "$LOG"; else
    echo "NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"; fi
  grep -a "^DETAIL \(FAIL\|ERROR\|WARN\)" "$o" >> "$LOG"
done < "$list"
echo "done $tier" >> "$LOG"
