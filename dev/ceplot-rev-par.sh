#!/usr/bin/env bash
# Reviewer driver (lane ceplot): one R process per job line
# "<arm> <package> <test file> [mutant]", N at a time.
# Output: dev/ceplot-rev-log/<tag>/<arm>--<pkg>--<file>[--<mutant>].txt
# and dev/ceplot-rev-log/<tag>.log with every RESULT line in list order.
#   bash dev/ceplot-rev-par.sh <tag> <job list> [N] [gated]
set -u
tag="$1"; list="$2"; N="${3:-16}"; gated="${4:-}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-ceplot
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/ceplot-rev-log/$tag"
LOG="$ROOT/dev/ceplot-rev-log/$tag.log"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/ceplot-rev-stan-cache"
export NOT_CRAN=true
if [ "$gated" = "gated" ]; then
  export FRMTMB_BRMS_FIT_TESTS=true
  export FRMTMB_DRMTMB_FIT_TESTS=true
  export FRMTMB_FUZZ=true
fi
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r arm p f m; do
  [ -n "$arm" ] || continue
  [ -f "$f" ] || { echo "missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  sfx=""; [ -n "${m:-}" ] && sfx="--$m"
  o="$OUT/$arm--$p--$(basename "$f" .R)$sfx.txt"
  ( cd "$(dirname "$f")" && "$RS" "$ROOT/dev/ceplot-rev-runtest.R" "$arm" "$p" "$f" ${m:-} > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r arm p f m; do
  [ -n "$arm" ] || continue
  want=$((want + 1))
  sfx=""; [ -n "${m:-}" ] && sfx="--$m"
  o="$OUT/$arm--$p--$(basename "$f" .R)$sfx.txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1 | tr -d '\r')
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$hit" >> "$LOG"
    grep -a -e "^  fails:" -e "^  warns:" "$o" | tr -d '\r' >> "$LOG"
  else
    echo "NO RESULT LINE for $arm $p $f ${m:-}" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "$tag ran $ran of $want" >> "$LOG"
