#!/bin/sh
# Lane catmulti: run a job list ("<package> <test file>" lines), one R
# process per file, N at a time, through dev/catmulti-test1.R on the
# lane library. Gated tiers: the caller exports FRMTMB_BRMS_FIT_TESTS.
#   sh dev/catmulti-run.sh <tag> <job list> [N]
# Output: dev/catmulti-out/<tag>-files/*.txt and dev/catmulti-out/<tag>.log
tag="$1"; list="$2"; N="${3:-8}"
WT=/c/Users/adf44/source/r/frmtmb-wt-catmulti
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/catmulti-out/$tag-files"
LOG="$WT/dev/catmulti-out/$tag.log"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache"
export NOT_CRAN=true
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
cd "$WT"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -p | wc -l)" -ge "$N" ]; do sleep 2; jobs > /dev/null; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/catmulti-test1.R lane "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else echo "$p NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"; fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
