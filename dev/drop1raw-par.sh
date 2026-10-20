#!/usr/bin/env bash
# Run "<package> <test file>" lines one R process per file, N at a time,
# through dev/drop1raw-runtest.R on the lane library.
#   bash dev/drop1raw-par.sh <tag> <job list> [N] [arm]
# Output: dev/drop1raw-out/<tag>-files/<pkg>--<file>.txt and
# dev/drop1raw-out/<tag>.log with every RESULT line.
set -u
tag="$1"; list="$2"; N="${3:-8}"; arm="${4:-lane}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/drop1raw-runtest.R"
OUT="$ROOT/dev/drop1raw-out/$tag-files"
LOG="$ROOT/dev/drop1raw-out/$tag.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$ROOT" && "$RS" "$RUN" "$arm" "$p" "$f" > "$o" 2>&1 ) &
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
    echo "$p NO RESULT LINE for $f" >> "$LOG"; tail -3 "$o" >> "$LOG"
  fi
done < "$list"
echo "ran $ran of $want" >> "$LOG"
echo "ran $ran of $want"
