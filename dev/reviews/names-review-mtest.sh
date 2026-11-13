#!/usr/bin/env bash
# Review of lane names: run the frmtmb.sample files either lane changed on
# the merged build, 3 processes at a time. Log dev/reviews/names-review-mtest.log
cd /c/Users/adf44/source/r/frmtmb-wt-names || exit 1
files="test-brms-output.R test-draws-methods.R test-latres.R
test-rescor-draws.R test-sample-direct.R test-sampling-ported.R
test-brms-draws-order.R test-draws-spellings.R test-generic-collision.R
test-gev-tmp-xi.R test-laplace-draws.R test-loglik-autocor.R
test-loglik-combine.R test-loo-newdata-fitted.R test-loo-newdata-subset.R
test-loo-predict.R test-loo.R test-postfit-draws.R test-ppcheck-loo.R
test-prior-draws.R test-stan-control.R"
out=dev/reviews/names-review-mtest
mkdir -p $out
run() { bash dev/reviews/names-review-run.sh dev/reviews/names-review-mtest.R "$out/$1.txt" "$1"; }
n=0
for f in $files; do
  run "$f" &
  n=$((n + 1))
  if [ $n -ge 3 ]; then wait -n; n=$((n - 1)); fi
done
wait
grep -h "^RESULT\|^FAILURE" $out/*.txt > dev/reviews/names-review-mtest.log
grep -L "^RESULT" $out/*.txt | sed 's/^/NO RESULT /' >> dev/reviews/names-review-mtest.log
echo "MTEST DONE" >> dev/reviews/names-review-mtest.log
