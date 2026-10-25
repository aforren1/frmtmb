#!/usr/bin/env bash
# Reviewer: run the listed test files, one R process each, N at a time.
# bash dev/review-sarfcor-tests.sh <lib> <list> <outdir> [N]
set -u
lib="$1"; list="$2"; out="$3"; N="${4:-4}"
WT=/c/Users/adf44/source/r/frmtmb-wt-sarfcor
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
mkdir -p "$out"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$out/$p--$(basename "$f" .R).txt"
  ( "$RS" "$WT/dev/review-sarfcor-run1.R" "$lib" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
for o in "$out"/*.txt; do
  r=$(grep -a "^RESULT" "$o" | tail -1)
  echo "$(basename "$o"): ${r:-NO RESULT}"
done
