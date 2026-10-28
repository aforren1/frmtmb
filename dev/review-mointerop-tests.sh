#!/usr/bin/env bash
# Reviewer: the requested test files, one R process each, 5 at a time.
set -u
ROOT="/c/Users/adf44/source/r/frmtmb-wt-mointerop"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/review-mointerop-run-test.R"
OUT="$ROOT/dev/review-mointerop-out/tests"
mkdir -p "$OUT"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
T="$ROOT/tests/testthat"
S="$ROOT/extensions/frmtmb.sample/tests/testthat"
jobs_list=(
  "lane frmtmb $T/test-mo-coding.R plain"
  "lane frmtmb $T/test-marginaleffects-args.R plain"
  "lane frmtmb $T/test-emmeans.R plain"
  "lane frmtmb $T/test-se-check.R plain"
  "lane frmtmb $T/test-message-uniqueness.R plain"
  "lane frmtmb $T/test-compat.R plain"
  "lane frmtmb.sample $S/test-mo-coding.R gated"
  "lane frmtmb.sample $S/test-loglik-autocor.R plain"
  "lane frmtmb $T/test-brms-suite-methods.R gated"
  "base frmtmb $T/test-marginaleffects-args.R plain"
  "base frmtmb $T/test-mo-coding.R plain"
)
for j in "${jobs_list[@]}"; do
  set -- $j
  arm=$1; p=$2; f=$3; gate=$4
  while [ "$(jobs -rp | wc -l)" -ge 5 ]; do sleep 2; done
  o="$OUT/$arm--$p--$(basename "$f" .R)--$gate.txt"
  (
    cd "$(dirname "$f")"
    if [ "$gate" = gated ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
    "$RS" "$RUN" "$arm" "$p" "$f" > "$o" 2>&1
  ) &
done
wait
grep -a -H "^RESULT\|^lib:" "$OUT"/*.txt
