#!/usr/bin/env bash
# Lane sparsese: dev/release/run-par.sh with its logs under this lane's
# prefix. One R process per test file, N at a time (default 2), through
# dev/release/run-tests.R.
#
#   bash dev/sparsese-par.sh <tier> <job list> [N]
#
# FRMTMB_RELLIB set: the lane's library ahead of rellib-r22 (the "after"
# arm); unset: rellib-r22 alone (the "before" arm).
# Output: dev/sparsese-suite-<tier>/<package>--<file>.txt per file and
# dev/sparsese-suite-<tier>/RESULTS.log with the RESULT lines in list order.
set -u
tier="$1"; list="$2"; N="${3:-2}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sparsese"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/sparsese-suite-$tier"
LOG="$OUT/RESULTS.log"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-C:/Users/adf44/Documents/.R/Makevars.win}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
export NOT_CRAN=true
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
cd "$ROOT"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do wait -n; done
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
    echo "NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "SUITE ran $ran of $want" >> "$LOG"
awk '/RESULT/ { for (i = 1; i <= NF; i++) { split($i, kv, "=");
  if (kv[1] in s || kv[1] ~ /^(pass|fail|err|skip|warn)$/) s[kv[1]] += kv[2] } }
  END { printf "TOTAL pass=%d fail=%d err=%d skip=%d warn=%d\n",
        s["pass"], s["fail"], s["err"], s["skip"], s["warn"] }' "$LOG" >> "$LOG"
tail -2 "$LOG"
