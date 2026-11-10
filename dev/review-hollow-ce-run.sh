#!/usr/bin/env bash
# Reviewer: core test files that call make_conditions or
# conditional_effects, one R process per file, 3 at a time, lane library
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-hollow
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-hollow/dev/stan-cache
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-hollow-lib
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
OUT=dev/review-hollow-ce
rm -rf "$OUT"; mkdir -p "$OUT"
i=0
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  i=$((i+1))
  ( R_PARALLEL_PORT=$((11900 + i)) "$RS" dev/release/run-tests.R "$p" "$f" > "$OUT/$(basename "$f" .R).txt" 2>&1 ) &
done < dev/review-hollow-ce-jobs.lst
wait
for o in "$OUT"/*.txt; do
  r=$(grep -a "^RESULT" "$o" | tail -1)
  l=$(grep -a "^lib:" "$o" | head -1)
  echo "${r:-NO RESULT $o} | ${l}"
done > dev/review-hollow-ce.log
echo DONE >> dev/review-hollow-ce.log
