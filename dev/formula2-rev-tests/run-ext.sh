#!/bin/bash
# Reviewer driver: every test file of the six other extensions, both arms
# (after = the lane's frmtmb under rellib-r3's extension builds), ungated.
cd /c/Users/adf44/source/r/frmtmb-wt-formula2
for p in frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.spline; do
  for f in $(ls extensions/$p/tests/testthat/test-*.R | xargs -n1 basename); do
    for arm in before after; do echo "$arm $p $f"; done
  done
done | xargs -P 22 -n 3 sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest-gated.R "$0" "$1" "$2" FALSE > "dev/formula2-rev-tests/x-$0-$1-$2.log" 2>&1'
echo ALLDONE
