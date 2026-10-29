#!/usr/bin/env bash
# Lane rescorpp: run a job list ("<package> <test file>" per line), one
# R process per file, N at a time, through dev/rescorpp-runtest.R.
#   bash dev/rescorpp-par.sh <arm> <tier> <job list> [N]
# Output: dev/rescorpp-tests/<tier>/<package>--<file>.txt per file and
# dev/rescorpp-tests/<tier>.log with the RESULT lines in list order.
set -u
arm="$1"; tier="$2"; list="$3"; N="${4:-8}"
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rescorpp
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/rescorpp-tests/$tier"
LOG="$ROOT/dev/rescorpp-tests/$tier.log"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$ROOT/dev/rescorpp-runtest.R" "$arm" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  libl=$(grep -a "^lib:" "$o" | head -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "NO RESULT LINE for $p $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
  echo "$p $libl" >> "$LOG.lib"
done < "$list"
echo "$tier ran $ran of $want" >> "$LOG"
