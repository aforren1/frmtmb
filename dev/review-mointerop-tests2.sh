#!/usr/bin/env bash
# Reviewer: short check, one R process per file, 4 at a time.
ROOT="/c/Users/adf44/source/r/frmtmb-wt-mointerop"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/review-mointerop-out/tests2"; mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
T="$ROOT/tests/testthat"
for f in test-mo-coding test-marginaleffects-args test-sandwich test-interop test-se-check test-confint-anova test-message-uniqueness; do
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 2; done
  ( cd "$T" && "$RS" "$ROOT/dev/review-mointerop-run-test.R" lane frmtmb "$T/$f.R" > "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -a -h "^RESULT\|^lib:" "$OUT"/*.txt
