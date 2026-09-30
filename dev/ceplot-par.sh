#!/usr/bin/env bash
# Lane ceplot: run a job list through dev/ceplot-runtest.R, one R process
# per file, N at a time. Each job line is "<arm> <package> <test file>".
# Output: dev/ceplot-log/<tag>/<arm>--<package>--<file>.txt, and
# dev/ceplot-log/<tag>.log with every RESULT line in list order.
#   bash dev/ceplot-par.sh <tag> <job list> [N]
set -u
tag="$1"; list="$2"; N="${3:-16}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-ceplot
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/ceplot-log/$tag"
LOG="$ROOT/dev/ceplot-log/$tag.log"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r arm p f; do
  [ -n "$arm" ] || continue
  [ -f "$f" ] || { echo "missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$arm--$p--$(basename "$f" .R).txt"
  ( cd "$(dirname "$f")" && "$RS" "$ROOT/dev/ceplot-runtest.R" "$arm" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r arm p f; do
  [ -n "$arm" ] || continue
  want=$((want + 1))
  o="$OUT/$arm--$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1 | tr -d '\r')
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$hit" >> "$LOG"
    grep -a -e "^  fails:" -e "^  warns:" "$o" | tr -d '\r' >> "$LOG"
  else
    echo "NO RESULT LINE for $arm $p $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "$tag ran $ran of $want" >> "$LOG"
