#!/usr/bin/env bash
# bash dev/review-splitsurv-par.sh <out dir> <job list> [N]
set -u
out="$1"; list="$2"; N="${3:-6}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-splitsurv
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/review-splitsurv-run.R"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
mkdir -p "$out"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$out/$p--$(basename "$f" .R)"
  ( "$RS" "$RUN" "$p" "$f" > "$b.txt" 2>&1 ) &
done < "$list"
wait
for t in "$out"/*.txt; do
  echo "$(basename "$t" .txt) $(grep -a '^lib:' "$t" | head -1 | sed 's/.*; frmtmb //') :: $(grep -a '^RESULT' "$t" | tail -1)"
done > "$out/summary.log"
