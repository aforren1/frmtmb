#!/usr/bin/env bash
# Review of lane loglik: run the lane's prior tests against each mutant
# of the review build (rev-loglik-lib, REV_MUT picks the mutant; "none"
# is the lane code in that build). Three R processes at most.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-loglik
OUT=$WT/dev/review-loglik-out/mut
mkdir -p "$OUT"
RS="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/rev-loglik-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
jobs_list=()
for m in none nojac twice nooffset halfnorm remlnoref; do
  for f in test-loglik-prior.R test-map.R test-prior-compat.R; do
    jobs_list+=("$m $f")
  done
done
port=31100
run_one() {
  local m=$1 f=$2 p=$3
  REV_MUT=$m R_PARALLEL_PORT=$p "$RS" "$WT/dev/loglik-run-tests.R" frmtmb \
    "$WT/tests/testthat/$f" > "$OUT/$m-$f.txt" 2>&1
}
n=0
for j in "${jobs_list[@]}"; do
  set -- $j
  port=$((port + 1))
  run_one "$1" "$2" "$port" &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
for j in "${jobs_list[@]}"; do
  set -- $j
  echo "MUT $1 $(grep -h '^lib:' "$OUT/$1-$2.txt" | head -1 | cut -c1-80) :: $(grep -h '^RESULT' "$OUT/$1-$2.txt")"
done
echo done
