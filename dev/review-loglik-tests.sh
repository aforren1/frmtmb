#!/usr/bin/env bash
# Review of lane loglik, attacks 7 and 8: the named core test files and
# the whole frmtmb.latent suite on the lane library, one file per R
# process, three at a time.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-loglik
OUT=$WT/dev/review-loglik-out/tests
mkdir -p "$OUT"
RS="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-loglik-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
jobs_list=()
for f in test-loglik-prior.R test-map.R test-nl-rtmb-scope.R \
         test-prior-compat.R test-refit-data.R \
         test-beta-binomial-density.R; do
  jobs_list+=("frmtmb $WT/tests/testthat/$f")
done
for f in $WT/extensions/frmtmb.latent/tests/testthat/test-*.R; do
  jobs_list+=("frmtmb.latent $f")
done
port=31300
n=0
for j in "${jobs_list[@]}"; do
  set -- $j
  port=$((port + 1))
  R_PARALLEL_PORT=$port "$RS" "$WT/dev/loglik-run-tests.R" "$1" "$2" \
    > "$OUT/$1-$(basename "$2").txt" 2>&1 &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
for j in "${jobs_list[@]}"; do
  set -- $j
  t="$OUT/$1-$(basename "$2").txt"
  echo "$(grep -h '^lib:' "$t" | head -1)"
  echo "  $(grep -h '^RESULT' "$t")"
done
echo done
