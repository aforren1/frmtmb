#!/usr/bin/env bash
# Review of lane superassign: the named test files on the lane library,
# one file per R process, two processes at a time (the third slot is kept
# for the review's own probes). Log per file in dev/review-superassign-t/.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-superassign
OUT=$WT/dev/review-superassign-t
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-superassign-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
jobs=(
  "frmtmb tests/testthat/test-no-superassign.R"
  "frmtmb tests/testthat/test-priors-autocor-classes.R"
  "frmtmb tests/testthat/test-priors-bounds-grcov.R"
  "frmtmb tests/testthat/test-compat.R"
  "frmtmb tests/testthat/test-nl.R"
  "frmtmb tests/testthat/test-me-newdata.R"
  "frmtmb tests/testthat/test-importance.R"
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-loo.R"
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-loo-predict.R"
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-sample-direct.R"
  "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-sample-seed.R"
)
i=0
for j in "${jobs[@]}"; do
  i=$((i + 1))
  printf '%s %s\n' "$i" "$j"
done | xargs -P 2 -L 1 bash -c '
  n=$0; p=$1; f=$2
  cd '"$WT"'
  R_PARALLEL_PORT=$((11730 + n)) "'"$RS"'" dev/superassign-run-tests.R "$p" "$f" \
    > "'"$OUT"'/$(basename "$f" .R).log" 2>&1
'
grep -h "^lib:\|^RESULT" "$OUT"/*.log
