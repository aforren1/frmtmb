#!/usr/bin/env bash
# Release review of 0.81.0: a few test files on rellib-r20 and rellib-r19,
# 3 R processes at a time. Logs dev/review-rel081-tf-<lib>-<file>.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel081
. dev/rel081-env.sh
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
run() {
  set -- $1
  lib=C:/Users/adf44/source/r/$1
  R_PARALLEL_PORT=$((11500 + RANDOM % 400)) "$RS" dev/review-rel081-tf.R "$lib" "$2" "$3" \
    > "dev/review-rel081-tf-$1-$(basename "$3" .R).txt" 2>&1
}
export -f run
T=tests/testthat
S=extensions/frmtmb.sample/tests/testthat
{
  echo "rellib-r20 frmtmb.sample $S/test-loo-newdata-fitted.R"
  echo "rellib-r20 frmtmb.sample $S/test-loo.R"
  echo "rellib-r20 frmtmb $T/test-se-tier.R"
  echo "rellib-r19 frmtmb $T/test-se-tier.R"
  echo "rellib-r20 frmtmb $T/test-families-extreme.R"
  echo "rellib-r20 frmtmb $T/test-cens-trunc.R"
  echo "rellib-r20 frmtmb $T/test-se-check.R"
  echo "rellib-r20 frmtmb $T/test-brms-suite-methods.R"
} | xargs -P 3 -I{} bash -c 'run "{}"'
grep -h "RESULT\|^lib:\|BAD" dev/review-rel081-tf-*.txt
echo TF DONE
