#!/usr/bin/env bash
# Run the core test suite ONE FILE PER PROCESS, three at a time.
#
# A whole-suite run in one process has repeatedly hidden state leakage
# in this project, and this lane's change is about what the search path
# looks like when a package loads, which is exactly the kind of state a
# shared process would smear.
#
#   bash dev/generics-suite.sh <outdir> [library to put FIRST]
set -u
OUT="$1"
EXTRA="${2:-}"
mkdir -p "$OUT"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-generics/dev/stan-cache"
cd "C:/Users/adf44/source/r/frmtmb-wt-generics" || exit 1

n=0
for f in tests/testthat/test-*.R; do
  b=$(basename "$f")
  "$R" dev/generics-runtests.R "$f" $EXTRA > "$OUT/$b.log" 2>&1 &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
echo "SUITE DONE"
