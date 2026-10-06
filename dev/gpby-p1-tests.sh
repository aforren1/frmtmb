#!/usr/bin/env bash
# Punch round 1: the test files the round touched or that reach
# predict()/extra_cov, one R process per file, N at a time.
#   bash dev/gpby-p1-tests.sh <out-dir> [gated]
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-gpby
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/$1"; N=12
export NOT_CRAN=true
if [ "${2:-}" = "gated" ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
rm -rf "$OUT"; mkdir -p "$OUT"
T=$ROOT/tests/testthat
SP=$ROOT/extensions/frmtmb.spline/tests/testthat
SA=$ROOT/extensions/frmtmb.sample/tests/testthat
cat > "$OUT/jobs.txt" <<EOF
frmtmb.spline $SP/test-deriv.R
frmtmb.spline $SP/test-curve.R
frmtmb.spline $SP/test-difference.R
frmtmb.spline $SP/test-new-levels.R
frmtmb.sample $SA/test-gp-by-draws.R
frmtmb $T/test-gp-by.R
frmtmb $T/test-fd-chain.R
frmtmb $T/test-emmeans.R
frmtmb $T/test-gp-multidim.R
frmtmb $T/test-predict-re-uncertainty.R
frmtmb $T/test-prior-compat.R
frmtmb $T/test-brms-likelihood.R
EOF
cd "$ROOT"
while read -r p f; do
  [ -f "$f" ] || { echo "MISSING $f"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/gpby-runtest.R "$p" "$f" > "$o" 2>&1 ) &
done < "$OUT/jobs.txt"
wait
for o in "$OUT"/*--*.txt; do
  echo "$(basename "$o" .txt): $(grep -a '^RESULT' "$o" | tail -1)"
  grep -a '^BAD' "$o"
done
