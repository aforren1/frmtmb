#!/usr/bin/env bash
# Reviewer: one R process per test file, N at a time.
#   bash dev/cifixrev-pardrv.sh <run name> <ref|ob> <libs> <job list> [N]
# <libs> is ';'-separated, user library appended by cifixrev-run1.R.
# Per-file logs in dev/cifixrev-log/<run>-suite/, summary in
# dev/cifixrev-log/<run>.sum. FRMTMB_REV_LOG is set per file when
# REVLOG=1 is in the environment.
set -u
run="$1"; kind="$2"; libs="$3"; list="$4"; N="${5:-16}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-cifix
cd "$ROOT"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP is unset}"
if [ "$kind" = ob ]; then
  RS="$ROOT/dev/cifix-out/Rob/bin/x64/Rscript.exe"
  export OPENBLAS_NUM_THREADS=4
else
  RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
fi
OUT="$ROOT/dev/cifixrev-log/$run-suite"
SUM="$ROOT/dev/cifixrev-log/$run.sum"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$SUM"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$p--$(basename "$f" .R)"
  o="$OUT/$b.txt"
  if [ "${REVLOG:-0}" = 1 ]; then
    ( FRMTMB_REV_LOG="$OUT/$b.rev" "$RS" dev/cifixrev-run1.R "$libs" "$p" "$f" > "$o" 2>&1 ) &
  else
    ( "$RS" dev/cifixrev-run1.R "$libs" "$p" "$f" > "$o" 2>&1 ) &
  fi
done < "$list"
wait
while read -r p f; do
  [ -n "$p" ] || continue
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$SUM"; else
    echo "$p NO RESULT LINE for $f" >> "$SUM"; tail -3 "$o" >> "$SUM"; fi
done < "$list"
echo "ran $(grep -c ' RESULT ' "$SUM") of $(grep -c . "$list")" >> "$SUM"
