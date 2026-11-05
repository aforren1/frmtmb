#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=review-impsamp-out/fr1/suites
mkdir -p $L
n=0
for p in frmtmb.coupling frmtmb.ode frmtmb.survival; do
  for f in ../extensions/$p/tests/testthat/test-*.R; do
    b=$(basename $f)
    "$RS" review-impsamp-fr1-suite.R $p $b > $L/$p-$b.log 2>&1 &
    n=$((n + 1))
    if [ $((n % 3)) -eq 0 ]; then wait; fi
  done
done
wait
echo done
